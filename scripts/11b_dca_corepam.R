# =============================================================================
# SCRIPT: 11b_dca_corepam.R
# PURPOSE: Decision Curve Analysis (DCA) for CorePAM score.
#          - OS/DSS: survival DCA at 60m (SCAN-B, METABRIC, GSE20685) and
#            24m (TCGA-BRCA) using dcurves::dca() with survival outcome.
#          - pCR: binary DCA using dcurves::dca() for the pCR cohorts.
#          Models compared: CorePAM score, age-only (CORE-A), age+score.
#
# Output:
#   results/supp/dca_os_summary.csv
#   results/supp/dca_pcr_summary.csv
#   figures/supp/en/pdf/FigS_DCA_OS_EN.pdf + figures/supp/en/png/FigS_DCA_OS_EN.png
#   figures/supp/pt/pdf/FigS_DCA_OS_PT.pdf + figures/supp/pt/png/FigS_DCA_OS_PT.png  — OS survival DCA
#   figures/supp/en/pdf/FigS_DCA_pCR_EN.pdf + figures/supp/en/png/FigS_DCA_pCR_EN.png
#   figures/supp/pt/pdf/FigS_DCA_pCR_PT.pdf + figures/supp/pt/png/FigS_DCA_pCR_PT.png — pCR binary DCA
#
# PROJECT: Core-PAM (Memorial v6.1 / Freeze Core-PAM)
# =============================================================================

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(survival)
  library(dcurves)
  library(ggplot2)
  library(patchwork)
})

source("scripts/00_setup.R")

SCRIPT_NAME <- "11b_dca_corepam.R"
FORCE       <- as.logical(Sys.getenv("FORCE_RERUN", "FALSE"))

out_os  <- file.path("results/supp", "dca_os_summary.csv")
out_pcr <- file.path("results/supp", "dca_pcr_summary.csv")
if (!FORCE && file.exists(out_os) && file.exists(out_pcr)) {
  message(sprintf("[%s] Outputs exist — skipping. Set FORCE_RERUN=TRUE.", SCRIPT_NAME))
  quit(save = "no", status = 0)
}

dir.create("results/supp", showWarnings = FALSE, recursive = TRUE)
# figure dirs created by 00_setup.R

message(sprintf("[%s] DCA analysis for OS and pCR", SCRIPT_NAME))

# --------------------------------------------------------------------------
# Helper: save PDF + PNG — bilingual (uses supp_en or supp_pt based on lang)
# --------------------------------------------------------------------------
save_fig <- function(p, name, w = 9, h = 6, lang = "en", dpi = 300) {
  pdf_path <- file.path(PATHS$figures[[paste0("supp_", lang, "_pdf")]], paste0(name, ".pdf"))
  png_path <- file.path(PATHS$figures[[paste0("supp_", lang, "_png")]], paste0(name, ".png"))
  old_warn <- getOption("warn"); options(warn = 0)
  tryCatch({
    cairo_pdf(pdf_path, width = w, height = h); print(p); dev.off()
    png(png_path, width = w, height = h, units = "in", res = dpi); print(p); dev.off()
    registry_append("ALL", name, pdf_path, sha256_file(pdf_path), "ok",
                    SCRIPT_NAME, file.info(pdf_path)$size / 1e6)
    message(sprintf("[%s] Saved: %s", SCRIPT_NAME, basename(pdf_path)))
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    message(sprintf("[%s] ERROR saving %s: %s", SCRIPT_NAME, name, e$message))
  })
  options(warn = old_warn)
  invisible(list(pdf = pdf_path, png = png_path))
}

# --------------------------------------------------------------------------
# PART 1: OS/DSS Decision Curve Analysis
# --------------------------------------------------------------------------
message(sprintf("[%s] === PART 1: OS/DSS DCA ===", SCRIPT_NAME))

os_cohort_cfg <- list(
  SCANB     = list(time_col = "os_time_months",  event_col = "os_event",
                   label = "SCAN-B (OS, 60m)",   horizon = 60),
  TCGA_BRCA = list(time_col = "os_time_months",  event_col = "os_event",
                   label = "TCGA-BRCA (OS, 24m)", horizon = 24),
  METABRIC  = list(time_col = "dss_time_months", event_col = "dss_event",
                   label = "METABRIC (DSS, 60m)", horizon = 60),
  GSE20685  = list(time_col = "os_time_months",  event_col = "os_event",
                   label = "GSE20685 (OS, 60m)",  horizon = 60)
)

dca_plots_os <- list()

for (cohort in names(os_cohort_cfg)) {
  cfg <- os_cohort_cfg[[cohort]]
  message(sprintf("[%s] %s: DCA at %dm", SCRIPT_NAME, cohort, cfg$horizon))

  pq_path <- file.path(proc_cohort(cohort), "analysis_ready.parquet")
  if (!file.exists(pq_path)) {
    message(sprintf("[%s] %s: parquet not found — skipping", SCRIPT_NAME, cohort))
    next
  }

  old_warn <- getOption("warn"); options(warn = 0)
  df <- arrow::read_parquet(pq_path)
  options(warn = old_warn)

  # Filter valid rows
  df <- df[!is.na(df[[cfg$time_col]]) & df[[cfg$time_col]] > 0 &
           !is.na(df[[cfg$event_col]]) & !is.na(df$score_z), ]

  message(sprintf("[%s] %s: N=%d events=%d", SCRIPT_NAME, cohort,
                  nrow(df), sum(df[[cfg$event_col]])))

  # Build CORE-A models
  has_age    <- "age" %in% names(df) && mean(!is.na(df$age)) > 0.5
  has_er     <- "er_status" %in% names(df) && mean(!is.na(df$er_status)) > 0.3

  time_vec  <- df[[cfg$time_col]]
  event_vec <- df[[cfg$event_col]]

  # Predicted probability using Cox model at horizon
  # Model 1: CorePAM score only
  old_warn2 <- getOption("warn"); options(warn = 0)
  cox1 <- tryCatch(coxph(Surv(time_vec, event_vec) ~ score_z, data = df, x = TRUE),
                   error = function(e) NULL)

  # Model 2: age only
  cox2 <- if (has_age) {
    tryCatch(coxph(Surv(time_vec, event_vec) ~ age, data = df, x = TRUE),
             error = function(e) NULL)
  } else NULL

  # Model 3: score_z + age (or score_z alone if no age)
  cox3 <- if (has_age) {
    tryCatch(coxph(Surv(time_vec, event_vec) ~ score_z + age, data = df, x = TRUE),
             error = function(e) NULL)
  } else cox1

  options(warn = old_warn2)

  if (is.null(cox1)) {
    message(sprintf("[%s] %s: Cox model failed — skipping DCA", SCRIPT_NAME, cohort))
    next
  }

  # Predicted event probability at horizon using survfit
  get_event_prob <- function(cox_fit, horizon_m) {
    tryCatch({
      sf <- survfit(cox_fit, newdata = df)
      # survfit returns survival probabilities per sample
      # Find time index closest to horizon
      times     <- summary(sf)$time
      surv_mat  <- t(summary(sf, times = horizon_m)$surv)
      event_prob <- 1 - as.vector(surv_mat)
      if (length(event_prob) == 1) {
        # Single summary time — expand to all patients
        # survfit stores individual predictions in matrix form
        event_prob <- rep(event_prob, nrow(df))
      }
      # Clamp to [0.001, 0.999]
      pmax(0.001, pmin(0.999, event_prob))
    }, error = function(e) NULL)
  }

  # Alternative: use Breslow estimator directly
  get_event_prob_v2 <- function(cox_fit, new_df, horizon_m) {
    tryCatch({
      lp <- predict(cox_fit, newdata = new_df, type = "lp")
      # Breslow baseline survival
      sf0 <- basehaz(cox_fit, centered = FALSE)
      # Find H0 at horizon
      idx <- max(which(sf0$time <= horizon_m), 1)
      h0  <- sf0$hazard[idx]
      S0  <- exp(-h0)
      # Individual survival at horizon
      S_i <- S0 ^ exp(lp)
      event_prob <- 1 - S_i
      pmax(0.001, pmin(0.999, event_prob))
    }, error = function(e) NULL)
  }

  p_score <- get_event_prob_v2(cox1, df, cfg$horizon)
  p_age   <- if (!is.null(cox2)) get_event_prob_v2(cox2, df, cfg$horizon) else NULL
  p_comb  <- if (!is.null(cox3)) get_event_prob_v2(cox3, df, cfg$horizon) else NULL

  if (is.null(p_score)) {
    message(sprintf("[%s] %s: prob estimation failed — skipping", SCRIPT_NAME, cohort))
    next
  }

  # Build data frame for dca
  dca_df <- data.frame(
    time      = time_vec,
    event     = event_vec,
    score_prb = p_score
  )
  if (!is.null(p_age))  dca_df$age_prb  <- p_age
  if (!is.null(p_comb)) dca_df$comb_prb <- p_comb

  # Run DCA
  old_warn3 <- getOption("warn"); options(warn = 0)
  dca_formula <- reformulate(
    c("score_prb",
      if (!is.null(p_age))  "age_prb"  else NULL,
      if (!is.null(p_comb)) "comb_prb" else NULL),
    response = sprintf("Surv(time, event)")
  )

  dca_res <- tryCatch(
    dca(dca_formula, data = dca_df, thresholds = seq(0, 0.5, by = 0.02),
        time = cfg$horizon),
    error = function(e) {
      message(sprintf("[%s] %s: dca() failed: %s", SCRIPT_NAME, cohort, e$message))
      NULL
    }
  )
  options(warn = old_warn3)

  if (is.null(dca_res)) next

  # Plot
  p_dca <- plot(dca_res, smooth = TRUE) +
    labs(
      title    = sprintf("DCA — %s", cfg$label),
      subtitle = sprintf("CorePAM, age (CORE-A), and combined at %d months", cfg$horizon),
      x        = "Decision threshold (probability)",
      y        = "Net benefit"
    ) +
    coord_cartesian(ylim = c(-0.05, NA)) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "bottom",
      plot.title      = element_text(face = "bold", size = 10),
      plot.subtitle   = element_text(size = 8),
      panel.grid.minor = element_blank()
    )

  dca_plots_os[[cohort]] <- p_dca

  message(sprintf("[%s] %s: DCA complete", SCRIPT_NAME, cohort))
}

# Combine OS DCA plots into a 2×2 panel
if (length(dca_plots_os) >= 2) {
  make_os_dca <- function(lang = "EN") {
    title_str <- if (lang == "EN") "Decision Curve Analysis — OS/DSS"
                 else "Análise de Curva de Decisão — OS/DSS"
    subtitle_str <- if (lang == "EN") "Net clinical benefit of CorePAM vs treat-all and treat-none"
                    else "Benefício clínico líquido do CorePAM vs tratar-todos e não-tratar"
    n_plots <- length(dca_plots_os)
    layout  <- if (n_plots == 4) c(2, 2) else c(2, ceiling(n_plots / 2))
    wrap_plots(dca_plots_os, ncol = layout[1]) +
      plot_annotation(
        title    = title_str,
        subtitle = subtitle_str,
        theme    = theme(
          plot.title    = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 10)
        )
      )
  }

  p_os_en <- make_os_dca("EN")
  p_os_pt <- make_os_dca("PT")

  save_fig(p_os_en, "FigS_DCA_OS_EN", w = 12, h = 9, lang = "en")
  save_fig(p_os_pt, "FigS_DCA_OS_PT", w = 12, h = 9, lang = "pt")
} else {
  message(sprintf("[%s] Not enough OS cohorts for combined DCA figure (%d)", SCRIPT_NAME,
                  length(dca_plots_os)))
}

# Save DCA net benefit at representative threshold (10%) per cohort for table
dca_os_summary <- tryCatch({
  nb_rows <- lapply(names(dca_plots_os), function(coh) {
    data.frame(
      cohort   = coh,
      horizon  = os_cohort_cfg[[coh]]$horizon,
      dca_done = TRUE,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, nb_rows)
}, error = function(e) NULL)

if (!is.null(dca_os_summary)) {
  write.csv(dca_os_summary, out_os, row.names = FALSE)
  message(sprintf("[%s] OS DCA summary: %s", SCRIPT_NAME, out_os))
}

# --------------------------------------------------------------------------
# PART 2: pCR Binary DCA
# --------------------------------------------------------------------------
message(sprintf("[%s] === PART 2: pCR Binary DCA ===", SCRIPT_NAME))

pcr_cohorts <- c("GSE25066", "GSE20194", "GSE32646", "ISPY1")
dca_plots_pcr <- list()

for (cohort in pcr_cohorts) {
  pq_path <- file.path(proc_pcr_cohort(cohort), "analysis_ready.parquet")
  if (!file.exists(pq_path)) next

  old_warn <- getOption("warn"); options(warn = 0)
  df <- arrow::read_parquet(pq_path)
  options(warn = old_warn)

  df <- df[!is.na(df$pcr) & !is.na(df$score_z), ]
  message(sprintf("[%s] %s: N=%d pCR=%d", SCRIPT_NAME, cohort, nrow(df), sum(df$pcr)))

  if (nrow(df) < 30 || sum(df$pcr) < 5) next

  # Predicted probability via logistic regression
  old_warn2 <- getOption("warn"); options(warn = 0)
  lr_score <- tryCatch(glm(pcr ~ score_z, data = df, family = binomial()), error = function(e) NULL)
  options(warn = old_warn2)

  if (is.null(lr_score)) next

  dca_df <- data.frame(
    pcr       = df$pcr,
    score_prb = predict(lr_score, type = "response")
  )

  old_warn3 <- getOption("warn"); options(warn = 0)
  dca_res <- tryCatch(
    dca(pcr ~ score_prb, data = dca_df, thresholds = seq(0, 0.6, by = 0.02)),
    error = function(e) {
      message(sprintf("[%s] %s pCR dca() failed: %s", SCRIPT_NAME, cohort, e$message))
      NULL
    }
  )
  options(warn = old_warn3)

  if (is.null(dca_res)) next

  pcr_rate_str <- sprintf("pCR rate: %.1f%%", 100 * mean(df$pcr))
  p_dca_pcr <- plot(dca_res, smooth = TRUE) +
    labs(
      title    = sprintf("%s", cohort),
      subtitle = pcr_rate_str,
      x        = "Decision threshold",
      y        = "Net benefit"
    ) +
    coord_cartesian(ylim = c(-0.05, NA)) +
    theme_bw(base_size = 10) +
    theme(
      legend.position  = "bottom",
      plot.title       = element_text(face = "bold", size = 10),
      plot.subtitle    = element_text(size = 8),
      panel.grid.minor = element_blank()
    )

  dca_plots_pcr[[cohort]] <- p_dca_pcr
  message(sprintf("[%s] %s: pCR DCA complete", SCRIPT_NAME, cohort))
}

if (length(dca_plots_pcr) >= 2) {
  make_pcr_dca <- function(lang = "EN") {
    title_str <- if (lang == "EN") "Decision Curve Analysis — pCR"
                 else "Análise de Curva de Decisão — pCR"
    subtitle_str <- if (lang == "EN")
      "Net clinical benefit of CorePAM score for pCR prediction"
      else "Benefício clínico líquido do CorePAM na predição de pCR"
    wrap_plots(dca_plots_pcr, ncol = 2) +
      plot_annotation(
        title    = title_str,
        subtitle = subtitle_str,
        theme    = theme(
          plot.title    = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 10)
        )
      )
  }

  p_pcr_en <- make_pcr_dca("EN")
  p_pcr_pt <- make_pcr_dca("PT")

  save_fig(p_pcr_en, "FigS_DCA_pCR_EN", w = 11, h = 8, lang = "en")
  save_fig(p_pcr_pt, "FigS_DCA_pCR_PT", w = 11, h = 8, lang = "pt")
}

# Save pCR DCA summary
dca_pcr_summary <- data.frame(
  cohort   = names(dca_plots_pcr),
  dca_done = TRUE,
  stringsAsFactors = FALSE
)
write.csv(dca_pcr_summary, out_pcr, row.names = FALSE)

h_os  <- sha256_file(out_os)
h_pcr <- sha256_file(out_pcr)
registry_append("ALL", "dca_os_summary",  out_os,  h_os,  "ok", SCRIPT_NAME,
                file.info(out_os)$size / 1e6)
registry_append("ALL", "dca_pcr_summary", out_pcr, h_pcr, "ok", SCRIPT_NAME,
                file.info(out_pcr)$size / 1e6)

message(sprintf("[%s] DONE | OS cohorts: %d | pCR cohorts: %d",
                SCRIPT_NAME, length(dca_plots_os), length(dca_plots_pcr)))
