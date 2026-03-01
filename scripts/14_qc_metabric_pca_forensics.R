# =============================================================================
# SCRIPT: 14_qc_metabric_pca_forensics.R
# PURPOSE: Forensic PCA on METABRIC expression data.
#          Detects batch effects, outlier samples, and platform artifacts
#          that could confound survival analysis.
#          Produces supplementary figure FigS4.
# PROJETO: Core-PAM (Memorial v6.1 §QC)
#
# INPUTS:
#   01_Base_Pura_CorePAM/PROCESSED/METABRIC/expression_genelevel_preZ.parquet
#   01_Base_Pura_CorePAM/PROCESSED/METABRIC/clinical_FINAL.parquet
#
# OUTPUTS:
#   figures/supp/en/pdf/FigS4_METABRIC_PCA_forensics.pdf
#   figures/supp/en/png/FigS4_METABRIC_PCA_forensics.png
#   results/supp/metabric_pca_forensics.csv
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "14_qc_metabric_pca_forensics.R"
COHORT      <- "METABRIC"

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

# Skip if output already exists (set FORCE_RERUN=TRUE to rerun)
FORCE <- as.logical(Sys.getenv("FORCE_RERUN", "FALSE"))
out_pdf <- file.path(PATHS$figures$supp_en_pdf, sprintf("FigS4_%s_PCA_forensics.pdf", COHORT))
if (!FORCE && file.exists(out_pdf)) {
  message(sprintf("[14] Output already exists: %s", out_pdf))
  message("[14] Skipping. Set FORCE_RERUN=TRUE to rerun.")
  quit(save = "no", status = 0)
}

message(sprintf("[14] Starting forensic PCA for %s", COHORT))

# ---------------------------------------------------------------------------
# 1) Load expression data
# ---------------------------------------------------------------------------
expr_path <- file.path(proc_cohort(COHORT), "expression_genelevel_preZ.parquet")
expr_df   <- strict_parquet(expr_path)

# Expect: first column = "gene", rest = sample columns
gene_col  <- names(expr_df)[1]
genes     <- expr_df[[gene_col]]
expr_mat  <- t(as.matrix(expr_df[, -1]))  # samples x genes
rownames(expr_mat) <- colnames(expr_df)[-1]
colnames(expr_mat) <- genes

message(sprintf("[14] Expression matrix: %d samples x %d genes", nrow(expr_mat), ncol(expr_mat)))

# ---------------------------------------------------------------------------
# 2) Select top 500 most variable genes for PCA
# ---------------------------------------------------------------------------
old_warn <- getOption("warn"); options(warn = 0)
gene_vars <- apply(expr_mat, 2, var, na.rm = TRUE)
options(warn = old_warn)

top_genes <- order(gene_vars, decreasing = TRUE)[seq_len(min(500L, ncol(expr_mat)))]
expr_top  <- expr_mat[, top_genes]

# Impute any remaining NAs with column means (rare after QC in step 03)
col_means <- colMeans(expr_top, na.rm = TRUE)
for (j in seq_len(ncol(expr_top))) {
  na_idx <- is.na(expr_top[, j])
  if (any(na_idx)) expr_top[na_idx, j] <- col_means[j]
}

message(sprintf("[14] PCA on top %d variable genes", ncol(expr_top)))

# ---------------------------------------------------------------------------
# 3) Run PCA
# ---------------------------------------------------------------------------
old_warn <- getOption("warn"); options(warn = 0)
pca_fit <- prcomp(expr_top, center = TRUE, scale. = TRUE)
options(warn = old_warn)

pct_var  <- 100 * pca_fit$sdev^2 / sum(pca_fit$sdev^2)
pc_scores <- as.data.frame(pca_fit$x[, 1:min(10L, ncol(pca_fit$x))])
pc_scores$sample_id <- rownames(pca_fit$x)

message(sprintf("[14] PC1: %.1f%% | PC2: %.1f%% | PC3: %.1f%% variance",
                pct_var[1], pct_var[2], pct_var[3]))

# ---------------------------------------------------------------------------
# 4) Load clinical data and annotate PCA
# ---------------------------------------------------------------------------
clin_path <- file.path(proc_cohort(COHORT), "clinical_FINAL.parquet")
clin_df   <- strict_parquet(clin_path)
clin_ids  <- normalize_id(clin_df$sample_id)
pc_scores$sample_norm <- normalize_id(pc_scores$sample_id)

annotated <- merge(pc_scores, clin_df, by.x = "sample_norm",
                   by.y = "sample_id", all.x = TRUE)
message(sprintf("[14] Annotated PCA scores: %d / %d matched to clinical",
                sum(!is.na(annotated$sample_id)), nrow(pc_scores)))

# ---------------------------------------------------------------------------
# 5) Detect outliers (Mahalanobis distance on PC1:PC3)
# ---------------------------------------------------------------------------
pc_matrix <- as.matrix(annotated[, c("PC1", "PC2", "PC3")])
cov_mat   <- cov(pc_matrix, use = "complete.obs")
center    <- colMeans(pc_matrix, na.rm = TRUE)

old_warn <- getOption("warn"); options(warn = 0)
maha_dist <- tryCatch(
  mahalanobis(pc_matrix, center = center, cov = cov_mat),
  error = function(e) rep(NA_real_, nrow(pc_matrix))
)
options(warn = old_warn)

# Chi-squared threshold for 3 df, p < 0.001
maha_threshold <- qchisq(0.999, df = 3)
annotated$mahal_dist  <- maha_dist
annotated$is_outlier  <- !is.na(maha_dist) & maha_dist > maha_threshold

n_outliers <- sum(annotated$is_outlier, na.rm = TRUE)
message(sprintf("[14] Mahalanobis outliers (p<0.001, df=3): %d / %d samples",
                n_outliers, nrow(annotated)))

if (n_outliers > 0) {
  out_ids <- annotated$sample_id[annotated$is_outlier]
  message(sprintf("[14] Outlier IDs: %s", paste(out_ids[seq_len(min(10L, length(out_ids)))],
                                                  collapse = ", ")))
}

# ---------------------------------------------------------------------------
# 6) PCA biplot (PC1 vs PC2), colored by ER status
# ---------------------------------------------------------------------------
# Note: er_status in clinical_FINAL.parquet has a bug — cBioPortal ER_IHC uses
# "Positve" (typo) not "Positive", so harmonize_clinical missed ER+ patients.
# Read ER directly from raw cBioPortal file here for correct coloring.
raw_clin_path <- "01_Base_Pura_CorePAM/RAW/METABRIC/brca_metabric/data_clinical_patient.txt"
er_col <- NULL
if (file.exists(raw_clin_path)) {
  raw_clin <- tryCatch(
    readr::read_tsv(raw_clin_path, skip = 4, show_col_types = FALSE),
    error = function(e) NULL
  )
  if (!is.null(raw_clin) && "ER_IHC" %in% names(raw_clin)) {
    er_map <- raw_clin[, c("PATIENT_ID", "ER_IHC")]
    names(er_map) <- c("sample_norm", "er_ihc_raw")
    er_map$sample_norm <- toupper(trimws(er_map$sample_norm))
    er_map$er_label_raw <- dplyr::case_when(
      er_map$er_ihc_raw == "Positve"  ~ "ER+",   # cBioPortal typo
      er_map$er_ihc_raw == "Negative" ~ "ER-",
      TRUE                            ~ NA_character_
    )
    annotated <- dplyr::left_join(annotated, er_map, by = "sample_norm")
    n_pos <- sum(annotated$er_label_raw == "ER+", na.rm = TRUE)
    n_neg <- sum(annotated$er_label_raw == "ER-", na.rm = TRUE)
    message(sprintf("[14] ER from raw cBioPortal: %d ER+ | %d ER-", n_pos, n_neg))
    er_col <- "er_label_raw"
  }
}

p1 <- ggplot(annotated, aes(x = PC1, y = PC2)) +
  {if (!is.null(er_col)) aes(color = as.factor(get(er_col))) else NULL} +
  geom_point(alpha = 0.5, size = 1.5) +
  {if (n_outliers > 0 && n_outliers <= 20) {
    geom_label_repel(
      data  = annotated[annotated$is_outlier, ],
      aes(label = sample_id),
      size  = 2.5,
      color = "red",
      max.overlaps = 20
    )
  }} +
  labs(
    title  = sprintf("METABRIC Forensic PCA — PC1 vs PC2 (top 500 variable genes)"),
    x      = sprintf("PC1 (%.1f%% variance)", pct_var[1]),
    y      = sprintf("PC2 (%.1f%% variance)", pct_var[2]),
    color  = "ER status",
    caption = sprintf("Outliers (Mahalanobis p<0.001): n=%d", n_outliers)
  ) +
  theme_classic(base_size = 11)

p2 <- ggplot(annotated, aes(x = PC1, y = PC3)) +
  {if (!is.null(er_col)) aes(color = as.factor(get(er_col))) else NULL} +
  geom_point(alpha = 0.5, size = 1.5) +
  labs(
    title  = "METABRIC Forensic PCA — PC1 vs PC3",
    x      = sprintf("PC1 (%.1f%%)", pct_var[1]),
    y      = sprintf("PC3 (%.1f%%)", pct_var[3]),
    color  = "ER status"
  ) +
  theme_classic(base_size = 11)

p3 <- ggplot(
  data.frame(pc = seq_len(10), pct = pct_var[seq_len(10)]),
  aes(x = pc, y = pct)
) +
  geom_col(fill = "#2980B9") +
  geom_line(color = "red", linewidth = 0.8) +
  geom_point(color = "red", size = 2) +
  labs(
    title = "Scree plot — PCA variance explained",
    x     = "Principal Component",
    y     = "Variance explained (%)"
  ) +
  theme_classic(base_size = 11)

# ---------------------------------------------------------------------------
# 7) Save figures
# ---------------------------------------------------------------------------
fig_dir <- PATHS$figures$supp
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

out_png <- file.path(fig_dir, sprintf("FigS4_%s_PCA_forensics.png", COHORT))

old_warn <- getOption("warn"); options(warn = 0)
pdf(out_pdf, width = 14, height = 5)
gridExtra::grid.arrange(p1, p2, p3, ncol = 3)
dev.off()

png(out_png, width = 1400, height = 500, res = 100)
gridExtra::grid.arrange(p1, p2, p3, ncol = 3)
dev.off()
options(warn = old_warn)

h_pdf <- sha256_file(out_pdf)
h_png <- sha256_file(out_png)
registry_append(COHORT, "figure_pca_forensics_pdf", out_pdf, h_pdf, "ok",
                SCRIPT_NAME, file.info(out_pdf)$size / 1e6)
registry_append(COHORT, "figure_pca_forensics_png", out_png, h_png, "ok",
                SCRIPT_NAME, file.info(out_png)$size / 1e6)

# ---------------------------------------------------------------------------
# 8) Save forensics report
# ---------------------------------------------------------------------------
forensics_df <- annotated[, c("sample_id", "PC1", "PC2", "PC3", "mahal_dist", "is_outlier")]
forensics_df$pct_var_pc1 <- round(pct_var[1], 2)
forensics_df$pct_var_pc2 <- round(pct_var[2], 2)
forensics_df$pct_var_pc3 <- round(pct_var[3], 2)
forensics_df$n_genes_pca  <- ncol(expr_top)

out_csv <- file.path(PATHS$results$supp, sprintf("metabric_pca_forensics.csv"))
readr::write_csv(forensics_df, out_csv)
h_csv <- sha256_file(out_csv)
registry_append(COHORT, "pca_forensics_report", out_csv, h_csv, "ok",
                SCRIPT_NAME, file.info(out_csv)$size / 1e6)

message(sprintf("[14] PCA forensics complete: %d outliers detected", n_outliers))
message(sprintf("[14] Saved: %s", out_pdf))
message("[14] Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
