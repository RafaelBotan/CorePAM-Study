# =============================================================================
# SCRIPT: 11_incremental_value_and_dca.R
# PURPOSE: Incremental value of CorePAM beyond CORE-A:
#          Delta C-index bootstrap (CORE-A vs CORE-A+score),
#          60m calibration (when available),
#          DCA as sensitivity only.
#          Follows Memorial v6.1 sec.2.1 + FIGURES sec.Fig5.
# PROJETO: Core-PAM (Memorial v6.1 / Freeze Core-PAM)
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "11_incremental_value_and_dca.R"

suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
})

set.seed(FREEZE$seed_folds)
message(sprintf("[%s] Starting CorePAM incremental value analysis", SCRIPT_NAME))

COHORTS <- c("SCANB", "TCGA_BRCA", "METABRIC", "GSE20685")

# Endpoint per cohort (from freeze)
ENDPOINT_MAP <- list(
  SCANB     = list(time = "os_time",  event = "os_event"),
  TCGA_BRCA = list(time = "os_time",  event = "os_event"),
  METABRIC  = list(time = "dss_time", event = "dss_event"),
  GSE20685  = list(time = "os_time",  event = "os_event")
)

# --------------------------------------------------------------------------
# Helper: Bootstrap Delta C-index (incremental value)
# --------------------------------------------------------------------------
bootstrap_delta_cindex <- function(time, event, score_z, corea_mat,
                                   n_boot = 1000, seed = 42) {
  set.seed(seed)
  n <- length(time)

  # Base C-index (CORE-A) and CORE-A+score
  old_warn <- getOption("warn"); options(warn = 0)
  c_base <- tryCatch(
    concordance(Surv(time, event) ~ corea_mat)$concordance,
    error = function(e) NA_real_
  )
  combined_mat <- cbind(corea_mat, score_z)
  c_full <- tryCatch(
    concordance(Surv(time, event) ~ combined_mat)$concordance,
    error = function(e) NA_real_
  )
  options(warn = old_warn)

  delta_obs <- c_full - c_base
  deltas    <- numeric(n_boot)

  old_warn <- getOption("warn"); options(warn = 0)
  for (i in seq_len(n_boot)) {
    idx <- sample(n, n, replace = TRUE)
    t_b <- time[idx]; e_b <- event[idx]; z_b <- score_z[idx]
    m_b <- corea_mat[idx, , drop = FALSE]
    cb  <- tryCatch(concordance(Surv(t_b, e_b) ~ m_b)$concordance, error = function(e) NA_real_)
    cf  <- tryCatch(concordance(Surv(t_b, e_b) ~ cbind(m_b, z_b))$concordance, error = function(e) NA_real_)
    deltas[i] <- cf - cb
  }
  options(warn = old_warn)

  list(
    c_base    = c_base,
    c_full    = c_full,
    delta     = delta_obs,
    ci_low    = quantile(deltas, 0.025, na.rm = TRUE),
    ci_high   = quantile(deltas, 0.975, na.rm = TRUE)
  )
}

# Helper: 60m calibration
calibrate_60m <- function(time, event, score_z, horizon = 60) {
  df_cal <- data.frame(time = time, event = event, score_z = score_z)
  df_cal <- df_cal[df_cal$time > 0 & !is.na(df_cal$score_z), ]
  if (nrow(df_cal) < 20) return(NULL)

  # Fit Cox and obtain predicted probability at 60m
  old_warn <- getOption("warn"); options(warn = 0)
  cox_cal <- tryCatch(coxph(Surv(time, event) ~ score_z, data = df_cal, x = TRUE),
                      error = function(e) NULL)
  options(warn = old_warn)
  if (is.null(cox_cal)) return(NULL)

  # Predictions at 60m
  old_warn <- getOption("warn"); options(warn = 0)
  sf  <- tryCatch(survfit(cox_cal, newdata = df_cal), error = function(e) NULL)
  options(warn = old_warn)
  if (is.null(sf)) return(NULL)

  # Get S(t=60m) for each observation
  times_sf <- sf$time
  if (max(times_sf) < horizon) return(NULL)

  t_idx <- max(which(times_sf <= horizon))
  pred_surv <- sf$surv[t_idx, ]

  # Observed: KM stratified by risk deciles
  df_cal$pred_surv <- pred_surv
  df_cal$decile    <- cut(df_cal$pred_surv,
                          breaks = quantile(df_cal$pred_surv, probs = seq(0, 1, 0.1), na.rm = TRUE),
                          include.lowest = TRUE, labels = FALSE)

  old_warn <- getOption("warn"); options(warn = 0)
  cal_obs <- lapply(split(df_cal, df_cal$decile), function(g) {
    if (nrow(g) < 3) return(NULL)
    km <- tryCatch(survfit(Surv(time, event) ~ 1, data = g), error = function(e) NULL)
    if (is.null(km)) return(NULL)
    t_km <- km$time; s_km <- km$surv
    obs_surv <- if (max(t_km) >= horizon) s_km[max(which(t_km <= horizon))] else NA_real_
    tibble(
      decile     = unique(g$decile),
      pred_mean  = mean(g$pred_surv, na.rm = TRUE),
      obs_surv   = obs_surv,
      n          = nrow(g)
    )
  })
  options(warn = old_warn)

  bind_rows(Filter(Negate(is.null), cal_obs))
}

# --------------------------------------------------------------------------
# 1) Loop over cohorts
# --------------------------------------------------------------------------
results_list <- list()

for (coh in COHORTS) {
  message(sprintf("[%s] Processing cohort: %s", SCRIPT_NAME, coh))

  ready_path <- file.path(proc_cohort(coh), "analysis_ready.parquet")
  if (!file.exists(ready_path)) {
    message(sprintf("[%s] WARNING: %s not found. Skipping.", SCRIPT_NAME, ready_path))
    next
  }

  df         <- strict_parquet(ready_path)
  ep         <- ENDPOINT_MAP[[coh]]
  time_col   <- ep$time
  event_col  <- ep$event

  # Fallback to OS if DSS not available
  if (!time_col %in% names(df)) {
    time_col  <- "os_time"
    event_col <- "os_event"
    message(sprintf("[%s] %s: primary endpoint unavailable; using OS", SCRIPT_NAME, coh))
  }

  if (!all(c(time_col, event_col, "score_z") %in% names(df))) {
    message(sprintf("[%s] %s: essential columns missing. Skipping.", SCRIPT_NAME, coh))
    next
  }

  df <- df[!is.na(df[[time_col]]) & df[[time_col]] > 0 &
             !is.na(df[[event_col]]) & !is.na(df$score_z), ]

  # CORE-A: age and/or ER
  corea_vars  <- c("age", "er_status")
  corea_avail <- intersect(corea_vars, names(df))

  c_base <- c_full <- delta <- ci_lo <- ci_hi <- NA_real_
  corea_used <- paste(corea_avail, collapse = "+")

  if (length(corea_avail) > 0) {
    df_m <- df[complete.cases(df[, c(corea_avail, "score_z", time_col, event_col)]), ]
    if (nrow(df_m) >= 30) {
      # Convert er_status to numeric if needed
      corea_mat <- df_m[, corea_avail, drop = FALSE]
      for (v in names(corea_mat)) {
        if (is.character(corea_mat[[v]]) || is.factor(corea_mat[[v]])) {
          corea_mat[[v]] <- as.numeric(factor(corea_mat[[v]]))
        }
      }
      corea_mat <- as.matrix(corea_mat)

      boot_delta <- bootstrap_delta_cindex(
        df_m[[time_col]], df_m[[event_col]], df_m$score_z, corea_mat,
        n_boot = as.integer(FREEZE$bootstrap_n),
        seed   = as.integer(FREEZE$seed_folds)
      )
      c_base <- boot_delta$c_base
      c_full <- boot_delta$c_full
      delta  <- boot_delta$delta
      ci_lo  <- boot_delta$ci_low
      ci_hi  <- boot_delta$ci_high
      message(sprintf("[%s] %s: C_base=%.4f | C_full=%.4f | DC=%.4f (%.4f, %.4f)",
                      SCRIPT_NAME, coh, c_base, c_full, delta, ci_lo, ci_hi))
    } else {
      message(sprintf("[%s] %s: n=%d insufficient for delta C-index", SCRIPT_NAME, coh, nrow(df_m)))
    }
  } else {
    message(sprintf("[%s] %s: CORE-A unavailable (no age/er_status)", SCRIPT_NAME, coh))
  }

  # 60m calibration
  cal_df    <- NULL
  cal_note  <- NA_character_
  has_cal60 <- max(df[[time_col]], na.rm = TRUE) >= 60

  if (has_cal60) {
    cal_df <- calibrate_60m(df[[time_col]], df[[event_col]], df$score_z, horizon = 60)
    if (!is.null(cal_df)) {
      cal_df$cohort <- coh
      cal_note <- "calibration_60m_ok"
      message(sprintf("[%s] %s: 60m calibration computed (%d deciles)", SCRIPT_NAME, coh, nrow(cal_df)))
    } else {
      cal_note <- "calibration_60m_failed"
    }
  } else {
    cal_note <- "follow_up_less_than_60m"
    message(sprintf("[%s] %s: maximum follow-up < 60m; calibration omitted", SCRIPT_NAME, coh))
  }

  results_list[[coh]] <- tibble(
    cohort          = coh,
    corea_vars_used = corea_used,
    n               = nrow(df),
    c_corea         = round(c_base, 4),
    c_corea_score   = round(c_full, 4),
    delta_cindex    = round(delta, 4),
    delta_ci_lo95   = round(ci_lo, 4),
    delta_ci_hi95   = round(ci_hi, 4),
    cal60_note      = cal_note
  )

  # Save calibration per cohort
  if (!is.null(cal_df)) {
    cal_path <- file.path(PATHS$results$supp, sprintf("calibration_60m_%s.csv", coh))
    readr::write_csv(cal_df, cal_path)
    h_cal <- sha256_file(cal_path)
    registry_append(coh, "calibration_60m", cal_path, h_cal, "ok", SCRIPT_NAME,
                    file.info(cal_path)$size / 1e6)
  }
}

# --------------------------------------------------------------------------
# 2) Save incremental value results
# --------------------------------------------------------------------------
if (length(results_list) == 0) {
  stop(sprintf("[%s] No incremental value results generated.", SCRIPT_NAME))
}

incr_df <- bind_rows(results_list)

main_dir <- PATHS$results$main
dir.create(main_dir, showWarnings = FALSE, recursive = TRUE)

incr_path <- file.path(main_dir, "incremental_value_by_cohort.csv")
readr::write_csv(incr_df, incr_path)
h_incr <- sha256_file(incr_path)
registry_append("ALL", "incremental_value", incr_path, h_incr, "ok", SCRIPT_NAME,
                file.info(incr_path)$size / 1e6)

# --------------------------------------------------------------------------
# 3) Figure Fig5 — Delta C-index per cohort
# --------------------------------------------------------------------------
fig_dir <- PATHS$figures$main
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

df_plot <- incr_df[!is.na(incr_df$delta_cindex), ]

if (nrow(df_plot) > 0) {
  p_delta <- ggplot(df_plot, aes(x = delta_cindex, y = reorder(cohort, delta_cindex))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = delta_ci_lo95, xmax = delta_ci_hi95),
                   height = 0.3, linewidth = 0.8, color = "#2980B9") +
    geom_point(size = 4, color = "#2980B9") +
    labs(
      title = "CorePAM Incremental Value: Delta C-index (CORE-A vs CORE-A + CorePAM)",
      x     = "Delta C-index (95% CI bootstrap)",
      y     = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 12))

  fig5_pdf <- file.path(fig_dir, "Fig5_DeltaCindex_COREA_vs_COREAplus_CorePAM.pdf")
  fig5_png <- file.path(fig_dir, "Fig5_DeltaCindex_COREA_vs_COREAplus_CorePAM.png")

  old_warn <- getOption("warn"); options(warn = 0)
  pdf(fig5_pdf, width = 8, height = 5); print(p_delta); dev.off()
  png(fig5_png, width = 800, height = 500, res = 100); print(p_delta); dev.off()
  options(warn = old_warn)

  registry_append("ALL", "figure_delta_cindex", fig5_pdf, sha256_file(fig5_pdf), "ok",
                  SCRIPT_NAME, file.info(fig5_pdf)$size / 1e6)
  registry_append("ALL", "figure_delta_cindex_png", fig5_png, sha256_file(fig5_png), "ok",
                  SCRIPT_NAME, file.info(fig5_png)$size / 1e6)
}

# --------------------------------------------------------------------------
# 4) Figure Fig5 — 60m calibration (panel, sensitivity)
# --------------------------------------------------------------------------
cal_files <- list.files(PATHS$results$supp, pattern = "^calibration_60m_.+\\.csv$",
                         full.names = TRUE)

if (length(cal_files) > 0) {
  old_warn <- getOption("warn"); options(warn = 0)
  cal_all <- bind_rows(lapply(cal_files, function(f) {
    tryCatch(readr::read_csv(f, show_col_types = FALSE), error = function(e) NULL)
  }))
  options(warn = old_warn)

  if (!is.null(cal_all) && nrow(cal_all) > 0 &&
      all(c("pred_mean", "obs_surv", "cohort") %in% names(cal_all))) {

    p_cal <- ggplot(cal_all, aes(x = pred_mean, y = obs_surv)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
      geom_point(aes(color = cohort), size = 2) +
      geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 0.7,
                  formula = y ~ x) +
      facet_wrap(~cohort) +
      labs(
        title = "CorePAM 60m Calibration (sensitivity)",
        x     = "Predicted survival at 60m",
        y     = "Observed survival at 60m (KM)"
      ) +
      theme_classic(base_size = 11) +
      theme(legend.position = "none")

    fig5b_pdf <- file.path(fig_dir, "Fig5_Calibration_60m_Panels_CorePAM.pdf")
    fig5b_png <- file.path(fig_dir, "Fig5_Calibration_60m_Panels_CorePAM.png")

    old_warn <- getOption("warn"); options(warn = 0)
    pdf(fig5b_pdf, width = 10, height = 6); print(p_cal); dev.off()
    png(fig5b_png, width = 1000, height = 600, res = 100); print(p_cal); dev.off()
    options(warn = old_warn)

    registry_append("ALL", "figure_calibration_60m", fig5b_pdf, sha256_file(fig5b_pdf), "ok",
                    SCRIPT_NAME, file.info(fig5b_pdf)$size / 1e6)
    registry_append("ALL", "figure_calibration_60m_png", fig5b_png, sha256_file(fig5b_png), "ok",
                    SCRIPT_NAME, file.info(fig5b_png)$size / 1e6)
  }
}

message(sprintf("[%s] COMPLETED | %d cohorts processed", SCRIPT_NAME, nrow(incr_df)))
