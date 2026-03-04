# =============================================================================
# SCRIPT: 34_subgroup_er_stratified.R
# PURPOSE: ER-stratified subgroup analysis in SCAN-B and METABRIC.
#          Computes Cox HR + C-index per ER subgroup (ER+ vs ER-),
#          tests interaction (score * ER) per cohort, and generates
#          a combined forest plot.
#          Addresses reviewer request for subgroup analysis by receptor status.
# OUTPUT:
#   results/supp/er_stratified_analysis.csv
#   figures/supp/en/pdf/FigS_ER_Stratified_Forest.pdf
#   figures/supp/en/png/FigS_ER_Stratified_Forest.png
# PROJECT: Core-PAM (Memorial CorePAM — Major Revision)
# =============================================================================

source("scripts/00_setup.R")
source("scripts/00_colors.R")
SCRIPT_NAME <- "34_subgroup_er_stratified.R"

suppressPackageStartupMessages({
  library(survival)
  library(arrow)
  library(ggplot2)
})

set.seed(FREEZE$seed_folds)
message(sprintf("[%s] ========== ER-STRATIFIED ANALYSIS ==========", SCRIPT_NAME))

# --------------------------------------------------------------------------
# 1) Bootstrap C-index helper
# --------------------------------------------------------------------------
bootstrap_cindex <- function(time, event, score_z, n_boot = 1000, seed = 42) {
  set.seed(seed)
  n <- length(time)
  cvals <- numeric(n_boot)
  old_warn <- getOption("warn"); options(warn = 0)
  for (i in seq_len(n_boot)) {
    idx <- sample(n, n, replace = TRUE)
    cx <- tryCatch(
      concordance(Surv(time[idx], event[idx]) ~ score_z[idx])$concordance,
      error = function(e) NA_real_
    )
    cvals[i] <- max(cx, 1 - cx, na.rm = TRUE)
  }
  options(warn = old_warn)
  c_raw <- concordance(Surv(time, event) ~ score_z)$concordance
  list(
    c_index = max(c_raw, 1 - c_raw),
    ci_low  = quantile(cvals, 0.025, na.rm = TRUE),
    ci_high = quantile(cvals, 0.975, na.rm = TRUE)
  )
}

# --------------------------------------------------------------------------
# 2) Per-cohort ER-stratified analysis function
# --------------------------------------------------------------------------
run_er_analysis <- function(cohort_name, parquet_path, time_col, event_col,
                            endpoint_label) {
  message(sprintf("[%s] --- %s (%s) ---", SCRIPT_NAME, cohort_name,
                  endpoint_label))

  df <- as.data.frame(strict_parquet(parquet_path))
  stopifnot(all(c(time_col, event_col, "score_z", "er_status") %in% names(df)))

  df <- df[!is.na(df[[time_col]]) & df[[time_col]] > 0 &
             !is.na(df[[event_col]]) & !is.na(df$score_z) &
             !is.na(df$er_status), ]

  # Standardize ER status to binary (handle different encodings)
  if (is.numeric(df$er_status)) {
    df$er_binary <- ifelse(df$er_status == 1, "ER+", "ER-")
  } else {
    df$er_binary <- ifelse(
      tolower(df$er_status) %in% c("positive", "er+", "1"),
      "ER+", "ER-"
    )
  }

  n_erp <- sum(df$er_binary == "ER+")
  n_ern <- sum(df$er_binary == "ER-")
  message(sprintf("[%s] %s N=%d with ER status | ER+: %d | ER-: %d",
                  SCRIPT_NAME, cohort_name, nrow(df), n_erp, n_ern))

  # --- Overall Cox ---
  old_warn <- getOption("warn"); options(warn = 0)
  cox_overall <- coxph(Surv(df[[time_col]], df[[event_col]]) ~ score_z,
                       data = df)
  options(warn = old_warn)
  sm_all <- summary(cox_overall)

  boot_all <- bootstrap_cindex(
    df[[time_col]], df[[event_col]], df$score_z,
    n_boot = as.integer(FREEZE$bootstrap_n),
    seed   = as.integer(FREEZE$seed_folds)
  )

  # --- Stratified analysis: ER+ and ER- ---
  results <- list()

  for (er_group in c("ER+", "ER-")) {
    df_sub <- df[df$er_binary == er_group, ]
    n_sub  <- nrow(df_sub)
    n_ev   <- sum(df_sub[[event_col]])

    message(sprintf("[%s] %s %s: N=%d, Events=%d",
                    SCRIPT_NAME, cohort_name, er_group, n_sub, n_ev))

    if (n_sub < 30 || n_ev < 10) {
      message(sprintf("[%s] %s %s: insufficient sample size, skipping",
                      SCRIPT_NAME, cohort_name, er_group))
      results[[er_group]] <- data.frame(
        cohort = cohort_name, endpoint = endpoint_label,
        subgroup = er_group, n = n_sub, n_events = n_ev,
        hr = NA, hr_lo95 = NA, hr_hi95 = NA, p_value = NA,
        c_index = NA, c_lo95 = NA, c_hi95 = NA,
        stringsAsFactors = FALSE
      )
      next
    }

    # Cox HR
    old_warn <- getOption("warn"); options(warn = 0)
    cox_sub <- tryCatch(
      coxph(Surv(df_sub[[time_col]], df_sub[[event_col]]) ~ score_z,
            data = df_sub),
      error = function(e) NULL
    )
    options(warn = old_warn)

    hr <- lo <- hi <- p_val <- NA_real_
    if (!is.null(cox_sub)) {
      sm <- summary(cox_sub)
      hr    <- sm$conf.int[1, "exp(coef)"]
      lo    <- sm$conf.int[1, "lower .95"]
      hi    <- sm$conf.int[1, "upper .95"]
      p_val <- sm$coefficients[1, "Pr(>|z|)"]
    }

    # Bootstrap C-index
    boot_sub <- bootstrap_cindex(
      df_sub[[time_col]], df_sub[[event_col]], df_sub$score_z,
      n_boot = as.integer(FREEZE$bootstrap_n),
      seed   = as.integer(FREEZE$seed_folds)
    )

    results[[er_group]] <- data.frame(
      cohort    = cohort_name,
      endpoint  = endpoint_label,
      subgroup  = er_group,
      n         = n_sub,
      n_events  = n_ev,
      hr        = round(hr, 4),
      hr_lo95   = round(lo, 4),
      hr_hi95   = round(hi, 4),
      p_value   = signif(p_val, 4),
      c_index   = round(boot_sub$c_index, 4),
      c_lo95    = round(boot_sub$ci_low, 4),
      c_hi95    = round(boot_sub$ci_high, 4),
      stringsAsFactors = FALSE
    )

    message(sprintf("[%s] %s %s: HR=%.3f (%.3f-%.3f) p=%.4g | C=%.4f",
                    SCRIPT_NAME, cohort_name, er_group,
                    hr, lo, hi, p_val, boot_sub$c_index))
  }

  # --- Interaction test: score_z * er_binary ---
  df$er_numeric <- ifelse(df$er_binary == "ER+", 1, 0)

  old_warn <- getOption("warn"); options(warn = 0)
  cox_interaction <- tryCatch(
    coxph(Surv(df[[time_col]], df[[event_col]]) ~ score_z * er_numeric,
          data = df),
    error = function(e) NULL
  )
  options(warn = old_warn)

  interaction_hr <- interaction_p <- NA_real_
  if (!is.null(cox_interaction)) {
    sm_int <- summary(cox_interaction)
    int_idx <- grep("score_z:er_numeric", rownames(sm_int$coefficients))
    if (length(int_idx) > 0) {
      interaction_hr <- round(exp(sm_int$coefficients[int_idx, "coef"]), 4)
      interaction_p  <- signif(sm_int$coefficients[int_idx, "Pr(>|z|)"], 4)
      message(sprintf("[%s] %s Interaction HR=%.4f, p=%.4g",
                      SCRIPT_NAME, cohort_name, interaction_hr, interaction_p))
    }
  }

  # --- Compile results ---
  res_df <- do.call(rbind, results)

  # Add overall row
  res_df <- rbind(
    data.frame(
      cohort   = cohort_name,
      endpoint = endpoint_label,
      subgroup = "Overall",
      n        = nrow(df),
      n_events = sum(df[[event_col]]),
      hr       = round(sm_all$conf.int[1, "exp(coef)"], 4),
      hr_lo95  = round(sm_all$conf.int[1, "lower .95"], 4),
      hr_hi95  = round(sm_all$conf.int[1, "upper .95"], 4),
      p_value  = signif(sm_all$coefficients[1, "Pr(>|z|)"], 4),
      c_index  = round(boot_all$c_index, 4),
      c_lo95   = round(boot_all$ci_low, 4),
      c_hi95   = round(boot_all$ci_high, 4),
      stringsAsFactors = FALSE
    ),
    res_df
  )

  # Attach interaction test as metadata columns
  res_df$interaction_hr <- interaction_hr
  res_df$interaction_p  <- interaction_p

  res_df
}

# --------------------------------------------------------------------------
# 3) Run analysis for both cohorts
# --------------------------------------------------------------------------
res_scanb <- run_er_analysis(
  cohort_name    = "SCAN-B",
  parquet_path   = file.path(proc_cohort("SCANB"), "analysis_ready.parquet"),
  time_col       = "os_time_months",
  event_col      = "os_event",
  endpoint_label = "OS"
)

res_metabric <- run_er_analysis(
  cohort_name    = "METABRIC",
  parquet_path   = file.path(proc_cohort("METABRIC"), "analysis_ready.parquet"),
  time_col       = "dss_time_months",
  event_col      = "dss_event",
  endpoint_label = "DSS"
)

# --------------------------------------------------------------------------
# 4) Combine and save results
# --------------------------------------------------------------------------
res_all <- rbind(res_scanb, res_metabric)

out_path <- file.path(PATHS$results$supp, "er_stratified_analysis.csv")
write.csv(res_all, out_path, row.names = FALSE)
h <- sha256_file(out_path)
registry_append("SCANB_METABRIC", "er_stratified_analysis", out_path, h,
                "ok", SCRIPT_NAME, file.info(out_path)$size / 1e6)
message(sprintf("[%s] Results saved: %s", SCRIPT_NAME, out_path))

# --------------------------------------------------------------------------
# 5) Combined forest plot: SCAN-B + METABRIC, ER-stratified
# --------------------------------------------------------------------------
plot_df <- res_all[!is.na(res_all$hr), ]
plot_df$label <- sprintf("%s %s (n=%d, ev=%d)",
                         plot_df$cohort, plot_df$subgroup,
                         plot_df$n, plot_df$n_events)
plot_df$hr_label <- sprintf("%.2f (%.2f\u2013%.2f)",
                            plot_df$hr, plot_df$hr_lo95, plot_df$hr_hi95)

# Row ordering: SCAN-B Overall / ER+ / ER- then METABRIC Overall / ER+ / ER-
row_order <- c(
  "SCAN-B Overall", "SCAN-B ER+", "SCAN-B ER-",
  "METABRIC Overall", "METABRIC ER+", "METABRIC ER-"
)
plot_df$row_id <- paste(plot_df$cohort, plot_df$subgroup)
plot_df$row_f  <- factor(plot_df$row_id, levels = rev(row_order))

# Assign cohort-specific colors
plot_df$pt_color <- ifelse(plot_df$cohort == "SCAN-B", COL$scanb, COL$metabric)

# Interaction annotation per cohort
int_scanb <- res_scanb[1, c("interaction_hr", "interaction_p")]
int_metab <- res_metabric[1, c("interaction_hr", "interaction_p")]

fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return(sprintf("%.2e", p))
  sprintf("%.3f", p)
}
subtitle_txt <- sprintf(
  "Interaction: SCAN-B HR=%.2f, p=%s | METABRIC HR=%.2f, p=%s",
  ifelse(is.na(int_scanb$interaction_hr), NA, int_scanb$interaction_hr),
  fmt_p(int_scanb$interaction_p),
  ifelse(is.na(int_metab$interaction_hr), NA, int_metab$interaction_hr),
  fmt_p(int_metab$interaction_p)
)

# Build y-axis labels with endpoint info
plot_df$y_label <- sprintf("%s %s\n(n=%d, ev=%d)",
                           plot_df$cohort, plot_df$subgroup,
                           plot_df$n, plot_df$n_events)
# Use same factor ordering for y_label
plot_df$y_label_f <- factor(plot_df$y_label,
                            levels = plot_df$y_label[match(rev(row_order),
                                                           plot_df$row_id)])

# Determine right margin x for HR labels
x_label_pos <- max(plot_df$hr_hi95, na.rm = TRUE) * 1.15

p_forest <- ggplot(plot_df, aes(x = hr, y = row_f, color = cohort)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  # Horizontal separator between cohorts
  geom_hline(yintercept = 3.5, linetype = "solid", color = "grey85",
             linewidth = 0.4) +
  geom_errorbar(aes(xmin = hr_lo95, xmax = hr_hi95),
                width = 0.2, linewidth = 0.8) +
  geom_point(aes(shape = ifelse(subgroup == "Overall", "diamond", "square"),
                 size  = ifelse(subgroup == "Overall", 5, 3.5))) +
  geom_text(aes(x = x_label_pos, label = hr_label),
            hjust = 0, size = 3.2, color = "grey30", show.legend = FALSE) +
  scale_shape_manual(values = c("diamond" = 18, "square" = 15), guide = "none") +
  scale_size_identity() +
  scale_color_manual(
    values = c("SCAN-B" = COL$scanb, "METABRIC" = COL$metabric),
    name   = "Cohort"
  ) +
  scale_y_discrete(labels = setNames(plot_df$y_label, plot_df$row_id)[rev(row_order)]) +
  scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3)) +
  coord_cartesian(xlim = c(0.3, max(plot_df$hr_hi95, na.rm = TRUE) * 2.5)) +
  labs(
    title    = "CorePAM HR per 1 SD by ER Status",
    subtitle = subtitle_txt,
    x        = "Hazard Ratio (log scale)",
    y        = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9.5, color = "grey40"),
    axis.text.y   = element_text(size = 10),
    legend.position = "bottom"
  )

fig_save(p_forest, "FigS_ER_Stratified_Forest", lang = "en", section = "supp",
         w = 10, h = 6)

message(sprintf("[%s] Forest plot saved", SCRIPT_NAME))
message(sprintf("[%s] ========== DONE ==========", SCRIPT_NAME))
