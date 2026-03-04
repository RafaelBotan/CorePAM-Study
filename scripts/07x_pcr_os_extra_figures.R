# =============================================================================
# SCRIPT: 07x_pcr_os_extra_figures.R
# PURPOSE: Generate extra OS-block supplementary figures:
#   Fig A: Raw CorePAM score distribution by OS cohort (violin + boxplot)
#   Fig B: CorePAM score_z vs PAM50-full score_z scatter (per cohort, LOESS)
#
# INPUTS:
#   results/corepam/CorePAM_weights.csv
#   results/corepam/CorePAM_model.rds
#   01_Base_Pura_CorePAM/PROCESSED/{COHORT}/expression_genelevel_preZ.parquet
#     (genes as rows, first col = "gene", remaining cols = sample IDs)
#
# OUTPUTS:
#   figures/supp/en/pdf/FigS_CorePAM_RawScore_Distribution_EN.pdf
#   figures/supp/en/png/FigS_CorePAM_RawScore_Distribution_EN.png
#   figures/supp/pt/pdf/FigS_CorePAM_RawScore_Distribution_PT.pdf
#   figures/supp/pt/png/FigS_CorePAM_RawScore_Distribution_PT.png
#   figures/supp/en/pdf/FigS_CorePAM_vs_PAM50full_EN.pdf
#   figures/supp/en/png/FigS_CorePAM_vs_PAM50full_EN.png
#   figures/supp/pt/pdf/FigS_CorePAM_vs_PAM50full_PT.pdf
#   figures/supp/pt/png/FigS_CorePAM_vs_PAM50full_PT.png
#
# PROJECT: Core-PAM (Memorial v6.1 / Freeze Core-PAM — OS block extra figs)
# =============================================================================

source("scripts/00_setup.R")
source("scripts/00_colors.R")
SCRIPT_NAME <- "07x_pcr_os_extra_figures.R"

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(glmnet)
})

FORCE <- as.logical(Sys.getenv("FORCE_RERUN", "FALSE"))

out_dist_en <- file.path(PATHS$figures$supp_en_pdf, "FigS_CorePAM_Score_Distribution_EN.pdf")
out_dist_pt <- file.path(PATHS$figures$supp_pt_pdf, "FigS_CorePAM_Score_Distribution_PT.pdf")
out_scat_en <- file.path(PATHS$figures$supp_en_pdf, "FigS_CorePAM_vs_PAM50full_EN.pdf")
out_scat_pt <- file.path(PATHS$figures$supp_pt_pdf, "FigS_CorePAM_vs_PAM50full_PT.pdf")

if (!FORCE && all(file.exists(c(out_dist_en, out_dist_pt, out_scat_en, out_scat_pt)))) {
  message(sprintf("[%s] All outputs exist — skipping. Set FORCE_RERUN=TRUE to rerun.",
                  SCRIPT_NAME))
  quit(save = "no", status = 0)
}

# --------------------------------------------------------------------------
# Load CorePAM weights
# --------------------------------------------------------------------------
weights_path <- file.path(PATHS$results$corepam, "CorePAM_weights.csv")
if (!file.exists(weights_path)) stop(sprintf("[%s] CorePAM_weights.csv not found: %s",
                                              SCRIPT_NAME, weights_path))
corepam_weights <- strict_csv(weights_path)
corepam_genes   <- corepam_weights$gene
w_vec           <- setNames(corepam_weights$weight, corepam_weights$gene)
message(sprintf("[%s] CorePAM weights loaded: %d genes", SCRIPT_NAME, length(corepam_genes)))

# --------------------------------------------------------------------------
# Extract PAM50-full weights from glmnet path at max df
# --------------------------------------------------------------------------
model_path <- file.path(PATHS$results$corepam, "CorePAM_model.rds")
if (!file.exists(model_path)) stop(sprintf("[%s] CorePAM_model.rds not found: %s",
                                            SCRIPT_NAME, model_path))
m <- strict_rds(model_path)
gfit <- m$glmnet_fit

# Find lambda at maximum df (use last occurrence = most regularized at max df)
df_path <- gfit$df
lambda_idx_max <- tail(which(df_path == max(df_path)), 1)
lambda_max <- gfit$lambda[lambda_idx_max]
message(sprintf("[%s] PAM50-full: using lambda index %d (df=%d, lambda=%.6f)",
                SCRIPT_NAME, lambda_idx_max, df_path[lambda_idx_max], lambda_max))

# Extract beta at max df — coef.glmnet returns sparse matrix
old_warn <- getOption("warn"); options(warn = 0)
beta_max_sparse <- glmnet::coef.glmnet(gfit, s = lambda_max)
options(warn = old_warn)

# Convert sparse matrix to named vector
beta_max_vec <- as.numeric(beta_max_sparse)
names(beta_max_vec) <- rownames(beta_max_sparse)
# Remove intercept and zero coefficients
beta_max_vec <- beta_max_vec[names(beta_max_vec) != "(Intercept)" & beta_max_vec != 0]
pam50full_genes   <- names(beta_max_vec)
pam50full_weights <- beta_max_vec
message(sprintf("[%s] PAM50-full weights extracted: %d genes", SCRIPT_NAME,
                length(pam50full_genes)))

# --------------------------------------------------------------------------
# OS cohorts
# --------------------------------------------------------------------------
OS_COHORTS <- c("SCANB", "TCGA_BRCA", "METABRIC", "GSE20685")

# --------------------------------------------------------------------------
# Helper: compute raw score from gene weights and expression matrix
# expr_mat: rows = samples, cols = genes (pre-Z or raw; we Z-score here)
# weights:  named numeric vector (gene -> weight)
# Returns: raw score vector (length = nrow(expr_mat))
# --------------------------------------------------------------------------
compute_score_from_expr <- function(expr_mat, weights) {
  genes_present <- intersect(names(weights), colnames(expr_mat))
  if (length(genes_present) == 0) return(rep(NA_real_, nrow(expr_mat)))
  w_sub <- weights[genes_present]
  # Z-score each gene intra-cohort
  e_sub <- scale(as.matrix(expr_mat[, genes_present, drop = FALSE]))
  # score_raw = sum(w_i * z_i) / sum(|w_i|)
  as.numeric((e_sub %*% w_sub) / sum(abs(w_sub)))
}

# --------------------------------------------------------------------------
# Load expression data, compute scores for each OS cohort
# expression_genelevel_preZ.parquet layout:
#   rows = genes (first col = "gene"), remaining cols = sample IDs
# We transpose to: rows = samples, cols = genes
# --------------------------------------------------------------------------
score_df_list <- list()

for (cohort in OS_COHORTS) {
  expr_path <- file.path(proc_cohort(cohort), "expression_genelevel_preZ.parquet")
  if (!file.exists(expr_path)) {
    message(sprintf("[%s] [%s] expression_genelevel_preZ.parquet not found — skipping.",
                    SCRIPT_NAME, cohort))
    next
  }

  message(sprintf("[%s] [%s] Loading expression matrix...", SCRIPT_NAME, cohort))
  old_warn <- getOption("warn"); options(warn = 0)
  expr_long <- arrow::read_parquet(expr_path)
  options(warn = old_warn)

  # expr_long: first col is "gene", rest are sample IDs
  gene_col_name <- names(expr_long)[1]  # should be "gene"
  gene_names    <- expr_long[[gene_col_name]]
  sample_ids    <- setdiff(names(expr_long), gene_col_name)

  message(sprintf("[%s] [%s] %d genes x %d samples loaded",
                  SCRIPT_NAME, cohort, length(gene_names), length(sample_ids)))

  # Subset to only genes needed (CorePAM + PAM50full) to save memory
  all_target_genes <- union(corepam_genes, pam50full_genes)
  keep_rows <- gene_names %in% all_target_genes
  if (sum(keep_rows) == 0) {
    message(sprintf("[%s] [%s] No target genes found in expression matrix — skipping.",
                    SCRIPT_NAME, cohort))
    next
  }
  expr_sub  <- expr_long[keep_rows, sample_ids, drop = FALSE]
  gene_sub  <- gene_names[keep_rows]
  message(sprintf("[%s] [%s] Using %d target genes", SCRIPT_NAME, cohort, length(gene_sub)))

  # Convert to numeric matrix: samples x genes
  expr_mat <- t(as.matrix(expr_sub))
  colnames(expr_mat) <- gene_sub
  rownames(expr_mat) <- sample_ids

  # Check CorePAM gene coverage
  corepam_present   <- intersect(corepam_genes, colnames(expr_mat))
  pam50full_present <- intersect(pam50full_genes, colnames(expr_mat))
  message(sprintf("[%s] [%s] CorePAM: %d/%d genes | PAM50-full: %d/%d genes",
                  SCRIPT_NAME, cohort,
                  length(corepam_present), length(corepam_genes),
                  length(pam50full_present), length(pam50full_genes)))

  if (length(corepam_present) < ceiling(0.8 * length(corepam_genes))) {
    message(sprintf("[%s] [%s] CorePAM gene coverage too low (%d/%d) — skipping.",
                    SCRIPT_NAME, cohort, length(corepam_present), length(corepam_genes)))
    next
  }

  # Compute raw scores (Z-scoring is done inside compute_score_from_expr)
  score_raw_corepam <- compute_score_from_expr(expr_mat, w_vec)
  # Overall Z-score of the raw score within cohort
  score_z_corepam   <- as.numeric(scale(score_raw_corepam))

  score_raw_pam50full <- if (length(pam50full_present) > 0) {
    compute_score_from_expr(expr_mat, pam50full_weights)
  } else {
    rep(NA_real_, nrow(expr_mat))
  }
  score_z_pam50full <- if (!all(is.na(score_raw_pam50full))) {
    as.numeric(scale(score_raw_pam50full))
  } else {
    rep(NA_real_, nrow(expr_mat))
  }

  score_df_list[[cohort]] <- data.frame(
    cohort            = cohort,
    score_raw         = score_raw_corepam,
    score_z_corepam   = score_z_corepam,
    score_z_pam50full = score_z_pam50full,
    stringsAsFactors  = FALSE
  )
  message(sprintf("[%s] [%s] Scores computed: %d samples", SCRIPT_NAME, cohort,
                  nrow(score_df_list[[cohort]])))

  # Free large matrix
  rm(expr_long, expr_sub, expr_mat); gc()
}

if (length(score_df_list) == 0) {
  stop(sprintf("[%s] No cohort data could be loaded — aborting.", SCRIPT_NAME))
}

all_scores <- bind_rows(score_df_list)
all_scores$cohort <- factor(all_scores$cohort, levels = OS_COHORTS)
message(sprintf("[%s] Total samples across all cohorts: %d", SCRIPT_NAME, nrow(all_scores)))

# --------------------------------------------------------------------------
# Colour map for OS cohorts
# --------------------------------------------------------------------------
COL_OS <- c(
  "SCANB"      = if (!is.null(COL$scanb))      COL$scanb      else "#2166AC",
  "TCGA_BRCA"  = if (!is.null(COL$tcga_brca))  COL$tcga_brca  else "#D6604D",
  "METABRIC"   = if (!is.null(COL$metabric))   COL$metabric   else "#4DAC26",
  "GSE20685"   = if (!is.null(COL$gse20685))   COL$gse20685   else "#762A83"
)

# --------------------------------------------------------------------------
# FIG A: Score distribution by cohort — violin + boxplot
# Plots score_z_corepam (= scale(Σwᵢzᵢ/Σ|wᵢ|), intra-cohort standardized)
# --------------------------------------------------------------------------
make_dist_fig <- function(lang = "EN") {
  xl  <- if (lang == "EN") "Cohort" else "Coorte"
  yl  <- if (lang == "EN") "CorePAM Score (z-standardized, SD units)" else
                            "Escore CorePAM (z-padronizado, unidades DP)"
  tt  <- if (lang == "EN") "CorePAM score distribution by OS cohort" else
                            "Distribuição do escore CorePAM por coorte OS"
  cap <- if (lang == "EN")
    "score_z = scale(\u03a3(w\u1d62\u00b7z\u1d62) / \u03a3|w\u1d62|), intra-cohort standardized; used in all Cox regressions"
  else
    "score_z = scale(\u03a3(w\u1d62\u00b7z\u1d62) / \u03a3|w\u1d62|), padronizado intracoorte; usado em todas as regressões Cox"

  ggplot(all_scores, aes(x = cohort, y = score_z_corepam, fill = cohort, colour = cohort)) +
    geom_violin(alpha = 0.35, linewidth = 0.5, trim = FALSE) +
    geom_boxplot(width = 0.18, alpha = 0.7, outlier.size = 0.5,
                 outlier.alpha = 0.4, colour = COL$black, fill = "white",
                 linewidth = 0.4) +
    scale_fill_manual(values = COL_OS, guide = "none") +
    scale_colour_manual(values = COL_OS, guide = "none") +
    labs(title = tt, x = xl, y = yl, caption = cap) +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(size = 10))
}

old_warn <- getOption("warn"); options(warn = 0)
for (lang in c("EN", "PT")) {
  lang_lc <- tolower(lang)
  pA    <- make_dist_fig(lang)
  pdf_f <- file.path(PATHS$figures[[paste0("supp_", lang_lc, "_pdf")]],
                     sprintf("FigS_CorePAM_Score_Distribution_%s.pdf", lang))
  png_f <- file.path(PATHS$figures[[paste0("supp_", lang_lc, "_png")]],
                     sprintf("FigS_CorePAM_Score_Distribution_%s.png", lang))
  cairo_pdf(pdf_f, width = 8, height = 5); print(pA); dev.off()
  png(png_f, width = 8, height = 5, units = "in", res = 600); print(pA); dev.off()
  registry_append("META_OS", sprintf("figS_corepam_score_dist_%s", lang),
                  pdf_f, sha256_file(pdf_f), "ok", SCRIPT_NAME,
                  file.info(pdf_f)$size / 1e6)
  message(sprintf("[%s] [%s] FigS_CorePAM_Score_Distribution saved: %s",
                  SCRIPT_NAME, lang, pdf_f))
}
gc(); options(warn = old_warn)

# --------------------------------------------------------------------------
# FIG B: CorePAM score_z vs PAM50-full score_z scatter — one panel per cohort
# --------------------------------------------------------------------------
scatter_data <- all_scores[!is.na(all_scores$score_z_pam50full), ]
message(sprintf("[%s] Scatter plot: %d samples with both CorePAM and PAM50-full scores",
                SCRIPT_NAME, nrow(scatter_data)))

make_scatter_fig <- function(lang = "EN") {
  xl  <- if (lang == "EN") "CorePAM score_z" else "Escore CorePAM score_z"
  yl  <- if (lang == "EN") "PAM50-full score_z" else "Escore PAM50-completo score_z"
  tt  <- if (lang == "EN") "CorePAM vs PAM50-full (within-cohort Z-scores)" else
                            "CorePAM vs PAM50-completo (Z-scores intracoorte)"
  cap <- if (lang == "EN")
    "Spearman \u03c1 shown per panel | LOESS smoother (span=0.75)"
  else
    "Spearman \u03c1 por painel | suavizador LOESS (span=0,75)"

  # Compute Spearman rho per cohort for annotation
  rho_df <- scatter_data |>
    dplyr::group_by(cohort) |>
    dplyr::summarise(
      rho   = cor(score_z_corepam, score_z_pam50full, method = "spearman",
                  use = "complete.obs"),
      x_pos = quantile(score_z_corepam, 0.05, na.rm = TRUE),
      y_pos = max(score_z_pam50full, na.rm = TRUE) * 0.92,
      .groups = "drop"
    ) |>
    dplyr::mutate(label = sprintf("\u03c1 = %.3f", rho))

  ggplot(scatter_data, aes(x = score_z_corepam, y = score_z_pam50full, colour = cohort)) +
    geom_point(alpha = 0.25, size = 0.8) +
    geom_smooth(method = "loess", span = 0.75, se = TRUE,
                colour = COL$black, fill = COL$grey_light,
                linewidth = 0.8) +
    geom_text(data = rho_df,
              aes(x = x_pos, y = y_pos, label = label),
              inherit.aes = FALSE,
              colour = COL$black, size = 3.5, hjust = 0, fontface = "bold") +
    facet_wrap(~ cohort, ncol = 2, scales = "free") +
    scale_colour_manual(values = COL_OS, guide = "none") +
    labs(title = tt, x = xl, y = yl, caption = cap) +
    theme_classic(base_size = 11) +
    theme(strip.text = element_text(size = 10, face = "bold"),
          strip.background = element_blank())
}

old_warn <- getOption("warn"); options(warn = 0)
if (nrow(scatter_data) > 0) {
  for (lang in c("EN", "PT")) {
    lang_lc <- tolower(lang)
    pB    <- make_scatter_fig(lang)
    pdf_f <- file.path(PATHS$figures[[paste0("supp_", lang_lc, "_pdf")]],
                       sprintf("FigS_CorePAM_vs_PAM50full_%s.pdf", lang))
    png_f <- file.path(PATHS$figures[[paste0("supp_", lang_lc, "_png")]],
                       sprintf("FigS_CorePAM_vs_PAM50full_%s.png", lang))
    cairo_pdf(pdf_f, width = 9, height = 8); print(pB); dev.off()
    png(png_f, width = 9, height = 8, units = "in", res = 600); print(pB); dev.off()
    registry_append("META_OS", sprintf("figS_corepam_vs_pam50full_%s", lang),
                    pdf_f, sha256_file(pdf_f), "ok", SCRIPT_NAME,
                    file.info(pdf_f)$size / 1e6)
    message(sprintf("[%s] [%s] FigS_CorePAM_vs_PAM50full saved: %s",
                    SCRIPT_NAME, lang, pdf_f))
  }
} else {
  message(sprintf("[%s] No valid scatter data — FigS_CorePAM_vs_PAM50full skipped.",
                  SCRIPT_NAME))
}
gc(); options(warn = old_warn)

message(sprintf("[%s] COMPLETED — extra OS supplementary figures saved to %s",
                SCRIPT_NAME, PATHS$figures$supp_en_pdf))
