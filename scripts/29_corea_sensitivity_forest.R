# =============================================================================
# SCRIPT: 29_corea_sensitivity_forest.R
# PURPOSE: CORE-A sensitivity forest plot — compares univariate HR (CorePAM
#          score only) vs adjusted HR (CORE-A: age +/- ER) across all
#          validation cohorts. Shows robustness of CorePAM prognostic value
#          after adjustment for clinical covariates.
#
# INPUTS:
#   results/supp/survival_results_{COHORT}.csv  (contains hr_uni and hr_multi)
#
# OUTPUTS:
#   results/supp/corea_sensitivity_forest.csv
#   figures/supp/en/png/FigS_COREA_Sensitivity_EN.png
#   figures/supp/en/pdf/FigS_COREA_Sensitivity_EN.pdf
#
# PROJECT: Core-PAM (Memorial v6.1 / Freeze Core-PAM)
# =============================================================================

source("scripts/00_setup.R")

SCRIPT_NAME <- "29_corea_sensitivity_forest.R"

# Skip/force pattern
FORCE <- as.logical(Sys.getenv("FORCE_RERUN", "FALSE"))
out_csv <- file.path(PATHS$results$supp, "corea_sensitivity_forest.csv")
if (!FORCE && file.exists(out_csv)) {
  message(sprintf("[%s] Output exists — skipping. Set FORCE_RERUN=TRUE.", SCRIPT_NAME))
  quit(save = "no", status = 0)
}

old_warn_pkg <- getOption("warn"); options(warn = 0)
suppressPackageStartupMessages({
  library(ggplot2)
})
options(warn = old_warn_pkg)

# =============================================================================
# 1) LOAD SURVIVAL RESULTS FROM ALL COHORTS
# =============================================================================
COHORTS <- c("SCANB", "TCGA_BRCA", "METABRIC", "GSE20685", "GSE1456")

results_list <- list()

for (cohort in COHORTS) {
  path <- file.path(PATHS$results$supp, paste0("survival_results_", cohort, ".csv"))
  if (!file.exists(path)) {
    message(sprintf("[%s] WARNING: %s not found — skipping", SCRIPT_NAME, basename(path)))
    next
  }

  sr <- strict_csv(path)

  # Univariate HR
  results_list[[length(results_list) + 1]] <- tibble::tibble(
    cohort       = cohort,
    model        = "Univariate",
    hr           = sr$hr_uni,
    hr_lo        = sr$hr_uni_lo95,
    hr_hi        = sr$hr_uni_hi95,
    p_value      = sr$p_uni,
    endpoint     = sr$endpoint,
    n_samples    = sr$n_samples,
    n_events     = sr$n_events,
    corea_vars   = "score_z only"
  )

  # Adjusted HR (CORE-A model)
  results_list[[length(results_list) + 1]] <- tibble::tibble(
    cohort       = cohort,
    model        = "Adjusted (CORE-A)",
    hr           = sr$hr_multi,
    hr_lo        = sr$hr_multi_lo95,
    hr_hi        = sr$hr_multi_hi95,
    p_value      = sr$p_multi,
    endpoint     = sr$endpoint,
    n_samples    = sr$n_samples,
    n_events     = sr$n_events,
    corea_vars   = sr$corea_vars_used
  )
}

result_df <- dplyr::bind_rows(results_list)

# =============================================================================
# 2) SAVE RESULTS
# =============================================================================
readr::write_csv(result_df, out_csv)
message(sprintf("[%s] Saved: %s (%d rows)", SCRIPT_NAME, basename(out_csv), nrow(result_df)))

# =============================================================================
# 3) FIGURE: Forest plot comparing univariate vs adjusted HRs
# =============================================================================

# Clean cohort labels and set order
cohort_order <- c("SCANB", "TCGA_BRCA", "METABRIC", "GSE20685", "GSE1456")
cohort_labels <- c(
  SCANB      = "SCAN-B (Training)",
  TCGA_BRCA  = "TCGA-BRCA",
  METABRIC   = "METABRIC",
  GSE20685   = "GSE20685",
  GSE1456    = "GSE1456"
)

plot_df <- result_df |>
  dplyr::mutate(
    cohort_label = cohort_labels[cohort],
    cohort_label = factor(cohort_label, levels = rev(cohort_labels)),
    model = factor(model, levels = c("Univariate", "Adjusted (CORE-A)")),
    # Annotations
    hr_label = sprintf("%.2f (%.2f-%.2f)", hr, hr_lo, hr_hi),
    sig_label = dplyr::case_when(
      p_value < 0.001  ~ "***",
      p_value < 0.01   ~ "**",
      p_value < 0.05   ~ "*",
      TRUE             ~ "ns"
    ),
    # Annotation with CORE-A variables used
    corea_note = ifelse(model == "Adjusted (CORE-A)",
                        paste0("  [", corea_vars, "]"), "")
  )

# Dodge position for side-by-side
dodge_width <- 0.5

# Pre-compute axis limits to avoid ggplot scoping warnings
x_lo <- min(plot_df$hr_lo, na.rm = TRUE) * 0.85
x_hi <- max(plot_df$hr_hi, na.rm = TRUE) * 1.55

p <- ggplot(plot_df, aes(x = hr, y = cohort_label, color = model, shape = model)) +
  # Reference line at HR=1
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.4) +

  # Error bars (CI) — using geom_errorbar with orientation
  geom_errorbar(aes(xmin = hr_lo, xmax = hr_hi),
                width = 0.2, linewidth = 0.6,
                orientation = "y",
                position = position_dodge(width = dodge_width)) +

  # Point estimates
  geom_point(size = 3.5,
             position = position_dodge(width = dodge_width)) +

  # Scale
  scale_x_log10(
    breaks = c(0.5, 0.75, 1, 1.25, 1.5, 2, 2.5),
    limits = c(x_lo, x_hi)
  ) +

  scale_color_manual(
    values = c("Univariate" = "#2166AC", "Adjusted (CORE-A)" = "#D6604D"),
    name   = "Model"
  ) +
  scale_shape_manual(
    values = c("Univariate" = 16, "Adjusted (CORE-A)" = 17),
    name   = "Model"
  ) +

  labs(
    title    = "CorePAM HR: Univariate vs CORE-A Adjusted",
    subtitle = "Hazard Ratio per 1-SD increase in CorePAM score (95% CI)",
    x        = "Hazard Ratio (log scale)",
    y        = NULL
  ) +

  theme_minimal(base_size = 11) +
  theme(
    legend.position    = "bottom",
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 9.5, color = "grey40"),
    axis.text.y        = element_text(size = 10)
  )

# Save
out_png <- file.path(PATHS$figures$supp_en_png, "FigS_COREA_Sensitivity_EN.png")
out_pdf <- file.path(PATHS$figures$supp_en_pdf, "FigS_COREA_Sensitivity_EN.pdf")

old_warn_gs <- getOption("warn"); options(warn = 0)
ggsave(out_png, p, width = 200, height = 140, units = "mm", dpi = 300, bg = "white")
ggsave(out_pdf, p, width = 200, height = 140, units = "mm", device = cairo_pdf)
options(warn = old_warn_gs)

message(sprintf("[%s] Figure saved: %s", SCRIPT_NAME, basename(out_png)))

# =============================================================================
# 4) REGISTRY
# =============================================================================
for (f in c(out_csv, out_png)) {
  h  <- sha256_file(f)
  sz <- file.info(f)$size / 1e6
  registry_append("MULTI", "corea_sensitivity", f, h, "ok", SCRIPT_NAME, sz, list())
}

# Summary
for (i in seq_len(nrow(result_df))) {
  r <- result_df[i, ]
  message(sprintf("[%s]   %s | %s | HR=%.3f (%.3f-%.3f) | p=%s",
                  SCRIPT_NAME, r$cohort, r$model, r$hr, r$hr_lo, r$hr_hi,
                  formatC(r$p_value, format = "e", digits = 2)))
}
message(sprintf("[%s] DONE", SCRIPT_NAME))
