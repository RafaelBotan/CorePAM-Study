# =============================================================================
# SCRIPT: 37_R2_km_multipanel.R
# PURPOSE: Generate consolidated multi-panel KM figure (new Figure 2) for R2
#          revision. Addresses Reviewer 1 comments:
#          - Single p-value per panel (log-rank for dichotomised comparison)
#          - Consolidated multi-panel (one cohort per panel)
#          - Consistent colour scheme, legend format, annotation style
# PROJECT: Core-PAM — Breast Cancer Research R2
# =============================================================================

source("scripts/00_setup.R")
source("scripts/00_colors.R")
SCRIPT_NAME <- "37_R2_km_multipanel.R"

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
  library(cowplot)
})

set.seed(FREEZE$seed_folds)
message(sprintf("[%s] Starting multi-panel KM figure generation", SCRIPT_NAME))

# --------------------------------------------------------------------------
# Cohort definitions (training first, then validation in manuscript order)
# --------------------------------------------------------------------------
make_cohort_defs <- function(lang = "EN") {
  tr <- if (lang == "EN") list(
    scanb_title = "SCAN-B (training)",
    ylab_os     = "Overall survival",
    ylab_dss    = "Disease-specific survival",
    gse1456     = "GSE1456 (Stockholm)",
    training_note = "Training cohort"
  ) else list(
    scanb_title = "SCAN-B (treinamento)",
    ylab_os     = "Sobrevida global",
    ylab_dss    = "Sobrevida doen\u00e7a-espec\u00edfica",
    gse1456     = "GSE1456 (Estocolmo)",
    training_note = "Coorte de treinamento"
  )
  list(
    list(id = "SCANB", label = "A", title = tr$scanb_title,
         time_col = "os_time_months", event_col = "os_event",
         ylab = tr$ylab_os, endpoint = "OS", note = tr$training_note),
    list(id = "TCGA_BRCA", label = "B", title = "TCGA-BRCA",
         time_col = "os_time_months", event_col = "os_event",
         ylab = tr$ylab_os, endpoint = "OS", note = NULL),
    list(id = "METABRIC", label = "C", title = "METABRIC",
         time_col = "dss_time_months", event_col = "dss_event",
         ylab = tr$ylab_dss, endpoint = "DSS", note = NULL),
    list(id = "GSE20685", label = "D", title = "GSE20685",
         time_col = "os_time_months", event_col = "os_event",
         ylab = tr$ylab_os, endpoint = "OS", note = NULL),
    list(id = "GSE1456", label = "E", title = tr$gse1456,
         time_col = "os_time_months", event_col = "os_event",
         ylab = tr$ylab_os, endpoint = "OS", note = NULL)
  )
}

# --------------------------------------------------------------------------
# Consistent styling
# --------------------------------------------------------------------------
PAL        <- c("High" = COL$km_high, "Low" = COL$km_low)
ANNO_SIZE  <- 2.8
TITLE_SIZE <- 11

lang_strings <- function(lang) {
  if (lang == "EN") list(
    legend = c("High risk", "Low risk"),
    xlab   = "Time (months)",
    events = "events",
    hr_prefix = "HR",
    lr_prefix = "Log-rank p",
    cutoff_note = "Cutoff: intra-cohort median\n(pre-specified, not optimised)"
  ) else list(
    legend = c("Alto risco", "Baixo risco"),
    xlab   = "Tempo (meses)",
    events = "eventos",
    hr_prefix = "HR",
    lr_prefix = "Log-rank p",
    cutoff_note = "Ponto de corte: mediana intracoorte\n(pr\u00e9-especificado, n\u00e3o otimizado)"
  )
}

# --------------------------------------------------------------------------
# Helper: x-axis cutoff (last time where min n-at-risk >= min_nrisk)
# --------------------------------------------------------------------------
x_cutoff <- function(km_fit, min_nrisk = 10) {
  s <- summary(km_fit)
  min_by_t <- tapply(s$n.risk, s$time, min)
  valid <- as.numeric(names(min_by_t[min_by_t >= min_nrisk]))
  if (length(valid) == 0) return(max(s$time))
  max(valid)
}

# --------------------------------------------------------------------------
# Generate one panel per cohort (wrapped in a language loop)
# --------------------------------------------------------------------------
build_all_panels <- function(lang) {
  cohort_defs <- make_cohort_defs(lang)
  tr <- lang_strings(lang)
  panels <- list()

  for (coh in cohort_defs) {
  message(sprintf("[%s] Processing %s (%s)...", SCRIPT_NAME, coh$id, coh$endpoint))

  # Load data
  ready_path <- file.path(proc_cohort(coh$id), "analysis_ready.parquet")
  df <- strict_parquet(ready_path)

  # Filter valid observations
  df <- df[!is.na(df[[coh$time_col]]) & df[[coh$time_col]] > 0 &
             !is.na(df[[coh$event_col]]) & !is.na(df$score_z), ]

  # Create standardised column names for survfit/ggsurvplot compatibility
  df$.time  <- df[[coh$time_col]]
  df$.event <- df[[coh$event_col]]

  n_total  <- nrow(df)
  n_events <- sum(df$.event)
  message(sprintf("[%s] %s: N=%d, Events=%d", SCRIPT_NAME, coh$id, n_total, n_events))

  # Median split
  km_cutpoint <- median(df$score_z, na.rm = TRUE)
  df$risk_group <- factor(
    ifelse(df$score_z >= km_cutpoint, "High", "Low"),
    levels = c("High", "Low")
  )

  # Dichotomised Cox HR (High vs Low)
  old_warn <- getOption("warn"); options(warn = 0)
  cox_km <- coxph(
    Surv(.time, .event) ~ I(risk_group == "High"),
    data = df
  )
  sm_km <- summary(cox_km)
  hr_km <- sm_km$conf.int[1, "exp(coef)"]
  lo_km <- sm_km$conf.int[1, "lower .95"]
  hi_km <- sm_km$conf.int[1, "upper .95"]

  # Log-rank test (single p-value for KM comparison)
  lr_test <- survdiff(
    Surv(.time, .event) ~ risk_group,
    data = df
  )
  lr_p <- 1 - pchisq(lr_test$chisq, df = 1)
  options(warn = old_warn)

  message(sprintf("[%s] %s: HR=%.2f (%.2f-%.2f), log-rank p=%.2e",
                  SCRIPT_NAME, coh$id, hr_km, lo_km, hi_km, lr_p))

  # KM fit
  old_warn <- getOption("warn"); options(warn = 0)
  km_fit <- survfit(
    Surv(.time, .event) ~ risk_group,
    data = df
  )
  options(warn = old_warn)

  # X-axis limit
  xlim_val  <- x_cutoff(km_fit, min_nrisk = 10)
  break_val <- max(12L, as.integer(round(xlim_val / 5 / 12) * 12))

  # Format p-value
  if (lr_p < 0.001) {
    p_label <- formatC(lr_p, format = "e", digits = 1)
  } else {
    p_label <- formatC(lr_p, format = "f", digits = 3)
  }

  # Annotation text — single annotation with HR + log-rank p
  anno_text <- sprintf(
    "%s = %.2f (%.2f\u2013%.2f)\n%s = %s",
    tr$hr_prefix, hr_km, lo_km, hi_km, tr$lr_prefix, p_label
  )

  # Panel title with cohort info
  panel_title <- sprintf(
    "%s.  %s (N = %s; %d %s; %s)",
    coh$label, coh$title,
    format(n_total, big.mark = ","), n_events, tr$events, coh$endpoint
  )

  # Create ggsurvplot (no risk table — too complex for multi-panel)
  old_warn <- getOption("warn"); options(warn = 0)
  km_p <- ggsurvplot(
    km_fit,
    data       = df,
    risk.table = FALSE,
    pval       = FALSE,
    conf.int   = TRUE,
    palette    = unname(PAL),
    xlab       = tr$xlab,
    ylab       = coh$ylab,
    legend.labs = tr$legend,
    legend.title = NULL,
    legend      = c(0.8, 0.95),
    xlim       = c(0, xlim_val),
    break.time.by = break_val,
    ggtheme    = theme_classic(base_size = 10)
  )
  options(warn = old_warn)

  # Add panel title and annotation
  km_p$plot <- km_p$plot +
    ggplot2::ggtitle(panel_title) +
    ggplot2::theme(
      plot.title = element_text(size = TITLE_SIZE, face = "bold", hjust = 0),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.4, "cm")
    ) +
    ggplot2::annotate(
      "text",
      x = xlim_val * 0.03, y = 0.12,
      label = anno_text,
      hjust = 0, size = ANNO_SIZE, colour = "grey20",
      lineheight = 0.9
    )

  panels[[coh$id]] <- km_p$plot
  }

  empty_panel <- ggplot() +
    theme_void() +
    annotate("text", x = 0.5, y = 0.5,
             label = tr$cutoff_note,
             size = 3.5, colour = "grey40", fontface = "italic")

  plot_grid(
    panels[["SCANB"]],
    panels[["TCGA_BRCA"]],
    panels[["METABRIC"]],
    panels[["GSE20685"]],
    panels[["GSE1456"]],
    empty_panel,
    ncol   = 2,
    nrow   = 3,
    align  = "hv",
    axis   = "tblr"
  )
}

# --------------------------------------------------------------------------
# Build + save for each language
# --------------------------------------------------------------------------
for (lang in c("EN", "PT")) {
  message(sprintf("[%s] Assembling multi-panel figure (%s)...", SCRIPT_NAME, lang))
  combined <- build_all_panels(lang)

  pdf_dir <- if (lang == "EN") PATHS$figures$main_en_pdf else PATHS$figures$main_pt_pdf
  png_dir <- if (lang == "EN") PATHS$figures$main_en_png else PATHS$figures$main_pt_png
  suffix <- lang
  out_pdf_main <- file.path(pdf_dir, sprintf("Fig2_KM_MultiPanel_%s.pdf", suffix))
  out_png_main <- file.path(png_dir, sprintf("Fig2_KM_MultiPanel_%s.png", suffix))

  cairo_pdf(out_pdf_main, width = 12, height = 16)
  print(combined)
  dev.off()
  message(sprintf("[%s] Saved: %s", SCRIPT_NAME, out_pdf_main))

  png(out_png_main, width = 12, height = 16, units = "in", res = 600)
  print(combined)
  dev.off()
  message(sprintf("[%s] Saved: %s", SCRIPT_NAME, out_png_main))

  if (lang == "EN") {
    out_png_sub <- file.path(ROOT_REPO, "submission", "figures", "Fig2_KM_MultiPanel_EN.png")
    file.copy(out_png_main, out_png_sub, overwrite = TRUE)
    message(sprintf("[%s] Copied to: %s", SCRIPT_NAME, out_png_sub))
    registry_append("ALL", "figure_km_multipanel_R2", out_pdf_main,
                    sha256_file(out_pdf_main), "ok", SCRIPT_NAME,
                    file.info(out_pdf_main)$size / 1e6)
  }
}

message(sprintf("[%s] COMPLETED — multi-panel KM figure generated (EN+PT)", SCRIPT_NAME))
