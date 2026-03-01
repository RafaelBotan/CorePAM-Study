# =============================================================================
# SCRIPT: 23b_pCR_ispy2_figures.R
# PURPOSE: Extended pCR forest plot including I-SPY2 as Exploratory External
#          Cohort. Generates bilingual figure (EN + PT).
#
# FIGURE:
#   Fig_pCR1_Extended_OR_{EN,PT}.pdf/png  — Forest: 4 primary + I-SPY2 + 2 pooled rows
#
# Visual design:
#   - Primary cohorts (4): solid circles, primary palette, sized by RE weight
#   - I-SPY2 (exploratory): triangle (pch=17), COL$ispy2, labelled "Exploratory"
#   - Separator line between primary and I-SPY2 sections
#   - Two pooled rows: RE without I-SPY2 (k=4) and RE with I-SPY2 (k=5)
#
# INPUTS:
#   results/pcr/meta_pCR_cohort_weights.csv          — primary cohort weights
#   results/pcr/meta_pCR_cohort_weights_with_ispy2.csv — all 5 cohorts
#   results/pcr/meta_pCR_without_ispy2.csv
#   results/pcr/meta_pCR_with_ispy2.csv
#   results/pcr/ispy2_results.csv
#
# PROJETO: Core-PAM (Memorial v6.1 / Freeze Core-PAM — pCR block)
# =============================================================================

source("scripts/00_setup.R")
source("scripts/00_colors.R")
SCRIPT_NAME <- "23b_pCR_ispy2_figures.R"

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

FORCE <- as.logical(Sys.getenv("FORCE_RERUN", "FALSE"))

fig_pcr <- PATHS$figures$pcr
dir.create(fig_pcr, showWarnings = FALSE, recursive = TRUE)

out_en <- file.path(fig_pcr, "Fig_pCR1_Extended_OR_EN.pdf")
out_pt <- file.path(fig_pcr, "Fig_pCR1_Extended_OR_PT.pdf")

if (!FORCE && file.exists(out_en) && file.exists(out_pt)) {
  message(sprintf("[%s] Outputs exist — skipping. Set FORCE_RERUN=TRUE to rerun.",
                  SCRIPT_NAME))
  quit(save = "no", status = 0)
}

# --------------------------------------------------------------------------
# 1) Load meta results
# --------------------------------------------------------------------------
req_files <- c(
  file.path(PATHS$results$pcr, "meta_pCR_cohort_weights.csv"),
  file.path(PATHS$results$pcr, "meta_pCR_cohort_weights_with_ispy2.csv"),
  file.path(PATHS$results$pcr, "meta_pCR_without_ispy2.csv"),
  file.path(PATHS$results$pcr, "meta_pCR_with_ispy2.csv"),
  file.path(PATHS$results$pcr, "ispy2_results.csv")
)
for (f in req_files) {
  if (!file.exists(f)) stop(sprintf("[%s] Missing: %s", SCRIPT_NAME, f))
}

# Primary cohort weights (from original 22_meta_pCR.R)
primary_wt <- strict_csv(file.path(PATHS$results$pcr, "meta_pCR_cohort_weights.csv"))

# I-SPY2 univariate result
ispy2_res <- strict_csv(file.path(PATHS$results$pcr, "ispy2_results.csv"))
ispy2_uni <- ispy2_res[ispy2_res$analysis == "univariate", ]

# Meta summaries
meta_no  <- strict_csv(file.path(PATHS$results$pcr, "meta_pCR_without_ispy2.csv"))
meta_yes <- strict_csv(file.path(PATHS$results$pcr, "meta_pCR_with_ispy2.csv"))

re_no   <- meta_no[meta_no$model == "RE", ]
re_yes  <- meta_yes[meta_yes$model == "RE", ]

message(sprintf("[%s] RE (k=4): OR=%.3f (%.3f–%.3f) I²=%.1f%%",
                SCRIPT_NAME, re_no$pooled_OR, re_no$pooled_OR_lo95, re_no$pooled_OR_hi95,
                re_no$I2_pct))
message(sprintf("[%s] RE (k=5): OR=%.3f (%.3f–%.3f) I²=%.1f%%",
                SCRIPT_NAME, re_yes$pooled_OR, re_yes$pooled_OR_lo95, re_yes$pooled_OR_hi95,
                re_yes$I2_pct))

# --------------------------------------------------------------------------
# 2) Build forest data
#    y positions:  primary cohorts top, I-SPY2 below separator, pooled rows bottom
# --------------------------------------------------------------------------
n_primary <- nrow(primary_wt)

# Primary cohort rows (y = n_primary down to 1)
primary_rows <- primary_wt |>
  dplyr::mutate(
    y_pos     = rev(seq_len(n_primary)) + 3,   # +3 to leave room for I-SPY2 + pooled
    row_type  = "primary",
    shape_val = "circle",
    size_val  = weight_re_pct / 10 + 2,
    label = sprintf("%s    N=%d, pCR=%d (%.0f%%)",
                    cohort, n_total, n_pcr1, 100 * n_pcr1 / n_total),
    color_cat = cohort
  ) |>
  dplyr::rename(or_val = or_uni, lo_val = or_lo95, hi_val = or_hi95)

# I-SPY2 row (y = 2)
ispy2_row <- tibble(
  cohort    = "ISPY2",
  label     = sprintf("I-SPY2 \u25b2    N=%d, pCR=%d (%.0f%%) [Exploratory]",
                      ispy2_uni$n_total, ispy2_uni$n_pcr1,
                      100 * ispy2_uni$pcr_rate),
  y_pos     = 2,
  or_val    = ispy2_uni$or_per_1SD,
  lo_val    = ispy2_uni$or_ci_lo95_wald,
  hi_val    = ispy2_uni$or_ci_hi95_wald,
  row_type  = "exploratory",
  shape_val = "triangle",
  size_val  = 3,
  color_cat = "ISPY2"
)

# Pooled rows (y = 1 and y = 0)
pool_no_lbl  <- sprintf("Pooled RE  (k=4, primary)   OR=%.3f (%.3f\u2013%.3f) I\u00b2=%.0f%%",
                        re_no$pooled_OR, re_no$pooled_OR_lo95, re_no$pooled_OR_hi95,
                        re_no$I2_pct)
pool_yes_lbl <- sprintf("Pooled RE  (k=5, +I-SPY2)  OR=%.3f (%.3f\u2013%.3f) I\u00b2=%.0f%%",
                        re_yes$pooled_OR, re_yes$pooled_OR_lo95, re_yes$pooled_OR_hi95,
                        re_yes$I2_pct)

pooled_rows <- tibble(
  cohort    = c("Pooled_k4", "Pooled_k5"),
  label     = c(pool_no_lbl, pool_yes_lbl),
  y_pos     = c(1, 0),
  or_val    = c(re_no$pooled_OR, re_yes$pooled_OR),
  lo_val    = c(re_no$pooled_OR_lo95, re_yes$pooled_OR_lo95),
  hi_val    = c(re_no$pooled_OR_hi95, re_yes$pooled_OR_hi95),
  row_type  = c("pooled_k4", "pooled_k5"),
  shape_val = c("diamond", "diamond"),
  size_val  = c(5, 5),
  color_cat = c("Pooled_k4", "Pooled_k5")
)

# Combine
forest_df <- dplyr::bind_rows(
  primary_rows[, c("cohort", "label", "y_pos", "or_val", "lo_val", "hi_val",
                    "row_type", "shape_val", "size_val", "color_cat")],
  ispy2_row,
  pooled_rows
)

# Color map
col_map <- c(
  "GSE25066"  = COL$gse25066,
  "GSE20194"  = COL$gse20194,
  "GSE32646"  = COL$gse32646,
  "ISPY1"     = COL$ispy1,
  "ISPY2"     = COL$ispy2,
  "Pooled_k4" = COL$summary,
  "Pooled_k5" = "#B22222"    # dark red to distinguish +ISPY2 pooled from k=4
)

# Y positions for separator lines
sep_y_primary  <- min(primary_rows$y_pos) - 0.5   # below primary cohorts
sep_y_pooled   <- 1.5                              # above pooled rows

# --------------------------------------------------------------------------
# 3) Forest plot function (bilingual)
# --------------------------------------------------------------------------
make_forest_extended <- function(lang = "EN") {
  xl  <- if (lang == "EN")
    "Odds Ratio (95% CI) per 1 SD CorePAM score_z"
  else
    "Odds Ratio (IC 95%) por 1 DP do escore_z CorePAM"

  ttl <- if (lang == "EN")
    "CorePAM score predicts pCR — primary + exploratory analysis"
  else
    "Escore CorePAM prediz pCR — an\u00e1lise prim\u00e1ria + explorat\u00f3ria"

  subt <- if (lang == "EN")
    "Random-effects (DerSimonian-Laird) | I-SPY2 = Exploratory External Cohort"
  else
    "Efeitos aleat\u00f3rios (DerSimonian-Laird) | I-SPY2 = Coorte Externa Explorat\u00f3ria"

  p <- ggplot(forest_df, aes(x = or_val, y = y_pos, colour = color_cat)) +
    # Reference line at OR=1
    geom_vline(xintercept = 1, linetype = "dashed",
               colour = COL$grey_mid, linewidth = 0.6) +
    # Separator lines
    geom_hline(yintercept = sep_y_primary, linetype = "solid",
               colour = COL$grey_mid, linewidth = 0.4) +
    geom_hline(yintercept = sep_y_pooled, linetype = "solid",
               colour = COL$grey_mid, linewidth = 0.4) +
    # Error bars
    geom_errorbarh(aes(xmin = lo_val, xmax = hi_val),
                   height = 0.25, linewidth = 0.65) +
    # Points (circles for primary, triangle for I-SPY2, diamond for pooled)
    geom_point(data = forest_df[forest_df$row_type == "primary", ],
               aes(size = size_val), shape = 16) +
    geom_point(data = forest_df[forest_df$row_type == "exploratory", ],
               aes(size = size_val), shape = 17) +
    geom_point(data = forest_df[grepl("^pooled", forest_df$row_type), ],
               aes(size = size_val), shape = 18) +
    scale_size_identity() +
    scale_colour_manual(values = col_map, guide = "none") +
    scale_x_log10(
      breaks = c(0.5, 0.75, 1, 1.25, 1.5, 2, 3),
      labels = c("0.5", "0.75", "1", "1.25", "1.5", "2", "3")
    ) +
    scale_y_continuous(
      breaks = forest_df$y_pos,
      labels = forest_df$label
    ) +
    labs(title = ttl, subtitle = subt, x = xl, y = NULL) +
    theme_classic(base_size = 10.5) +
    theme(
      axis.text.y  = element_text(size = 8.5),
      plot.subtitle = element_text(size = 8, colour = COL$grey_dark),
      panel.grid.major.x = element_line(colour = COL$grey_light, linewidth = 0.3)
    )

  p
}

# --------------------------------------------------------------------------
# 4) Save figures
# --------------------------------------------------------------------------
old_warn <- getOption("warn"); options(warn = 0)
for (lang in c("EN", "PT")) {
  p_ext  <- make_forest_extended(lang)
  pdf_f  <- file.path(fig_pcr, sprintf("Fig_pCR1_Extended_OR_%s.pdf", lang))
  png_f  <- file.path(fig_pcr, sprintf("Fig_pCR1_Extended_OR_%s.png", lang))
  cairo_pdf(pdf_f, width = 10, height = 5.5); print(p_ext); dev.off()
  png(png_f, width = 10, height = 5.5, units = "in", res = 600); print(p_ext); dev.off()
  h  <- sha256_file(pdf_f)
  sz <- file.info(pdf_f)$size / 1e6
  registry_append("META_PCR_ISPY2", sprintf("fig_pcr1_extended_%s", lang),
                  pdf_f, h, "ok", SCRIPT_NAME, sz)
  message(sprintf("[%s] [%s] Saved: %s", SCRIPT_NAME, lang, pdf_f))
}
options(warn = old_warn)

message(sprintf("[%s] COMPLETED.", SCRIPT_NAME))
