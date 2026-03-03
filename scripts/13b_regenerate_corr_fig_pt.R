# =============================================================================
# SCRIPT: 13b_regenerate_corr_fig_pt.R
# PURPOSE: Regenerate Portuguese version of FigS3 (Correlation Off-Diagonal)
#          with corrected labels: Spearman -> Pearson (matching EN version).
#          The underlying data uses Pearson correlation on quantile values,
#          NOT Spearman rank correlation.
#
# OUTPUT:
#   - figures/supp/pt/png/FigS3_Correlation_OffDiagonal_PT.png
#   - figures/supp/pt/pdf/FigS3_Correlation_OffDiagonal_PT.pdf
#
# PROJETO: Core-PAM (Memorial v6.1 / Freeze Core-PAM)
# =============================================================================

source("scripts/00_setup.R")
options(warn = 0)
SCRIPT_NAME <- "13b_regenerate_corr_fig_pt.R"

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
})

message(sprintf("[%s] Regenerating FigS3 PT (Pearson on quantile values)", SCRIPT_NAME))

# --------------------------------------------------------------------------
# 1) Load correlation matrix from CSV
# --------------------------------------------------------------------------
cor_csv <- file.path(PATHS$results$supp, "qc_score_correlations_offdiag.csv")
stopifnot(file.exists(cor_csv))

cor_df       <- read_csv(cor_csv, show_col_types = FALSE)
cohort_names <- setdiff(names(cor_df), "cohort_row")
cor_mat      <- as.matrix(cor_df[, cohort_names])
rownames(cor_mat) <- cor_df$cohort_row

message(sprintf("[%s] Correlation matrix loaded (%d cohorts)", SCRIPT_NAME, length(cohort_names)))

# --------------------------------------------------------------------------
# 2) Build long-format data (off-diagonal only)
# --------------------------------------------------------------------------
cor_long <- do.call(rbind, lapply(cohort_names, function(i) {
  do.call(rbind, lapply(cohort_names, function(j) {
    data.frame(
      Coorte_A   = i,
      Coorte_B   = j,
      Pearson_r  = cor_mat[i, j],
      stringsAsFactors = FALSE
    )
  }))
}))

# --------------------------------------------------------------------------
# 3) Plot — Portuguese labels, Pearson (corrected from Spearman)
# --------------------------------------------------------------------------
p_cor <- ggplot(cor_long, aes(x = Coorte_A, y = Coorte_B, fill = Pearson_r)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(
    aes(label = ifelse(is.na(Pearson_r), "", sprintf("%.3f", Pearson_r))),
    color = "black", size = 3.5
  ) +
  scale_fill_gradient2(
    low      = "#2166AC",
    mid      = "white",
    high     = "#D6604D",
    midpoint = 0,
    limits   = c(-1, 1),
    na.value = "gray90",
    name     = "Pearson r\n(valores dos quantis)"
  ) +
  labs(
    title    = "Correla\u00e7\u00f5es do escore CorePAM (fora da diagonal, via valores dos quantis)",
    subtitle = "Pearson sobre valores dos quantis | Diagonal exclu\u00edda",
    x        = NULL,
    y        = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x   = element_text(angle = 45, hjust = 1),
    plot.title     = element_text(face = "bold"),
    plot.subtitle  = element_text(color = "gray40")
  )

# --------------------------------------------------------------------------
# 4) Save PDF + PNG
# --------------------------------------------------------------------------
pdf_dir <- PATHS$figures$supp_pt_pdf
png_dir <- PATHS$figures$supp_pt_png
dir.create(pdf_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(png_dir, showWarnings = FALSE, recursive = TRUE)

fig_name <- "FigS3_Correlation_OffDiagonal_PT"
pdf_path <- file.path(pdf_dir, paste0(fig_name, ".pdf"))
png_path <- file.path(png_dir, paste0(fig_name, ".png"))

dpi <- 300
cairo_pdf(pdf_path, width = 7, height = 6)
print(p_cor)
dev.off()

png(png_path, width = round(7 * dpi), height = round(6 * dpi), res = dpi)
print(p_cor)
dev.off()

message(sprintf("[%s] Saved: %s", SCRIPT_NAME, pdf_path))
message(sprintf("[%s] Saved: %s", SCRIPT_NAME, png_path))

# --------------------------------------------------------------------------
# 5) Registry
# --------------------------------------------------------------------------
h_pdf <- sha256_file(pdf_path)
h_png <- sha256_file(png_path)
registry_append("ALL", "figure_qc_corr_offdiag_pt", pdf_path, h_pdf, "ok", SCRIPT_NAME,
                file.info(pdf_path)$size / 1e6)
registry_append("ALL", "figure_qc_corr_offdiag_pt_png", png_path, h_png, "ok", SCRIPT_NAME,
                file.info(png_path)$size / 1e6)

message(sprintf("[%s] COMPLETED | FigS3 PT regenerated with Pearson labels", SCRIPT_NAME))
