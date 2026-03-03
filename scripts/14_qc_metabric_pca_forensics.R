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

# --- Helper: build PCA plots with language support ---
make_pca_plots <- function(lang = "EN") {
  if (lang == "PT") {
    t1 <- sprintf("PCA Forense METABRIC — PC1 vs PC2 (500 genes mais variáveis)")
    t2 <- "PCA Forense METABRIC — PC1 vs PC3"
    t3 <- "Scree plot — Variância explicada pela PCA"
    x1 <- sprintf("PC1 (%.1f%% variância)", pct_var[1])
    y1 <- sprintf("PC2 (%.1f%% variância)", pct_var[2])
    x2 <- sprintf("PC1 (%.1f%%)", pct_var[1])
    y2 <- sprintf("PC3 (%.1f%%)", pct_var[3])
    x3 <- "Componente Principal"
    y3 <- "Variância explicada (%)"
    col_lab <- "Status RE"
    cap1 <- sprintf("Outliers (Mahalanobis p<0,001): n=%d", n_outliers)
  } else {
    t1 <- sprintf("METABRIC Forensic PCA — PC1 vs PC2 (top 500 variable genes)")
    t2 <- "METABRIC Forensic PCA — PC1 vs PC3"
    t3 <- "Scree plot — PCA variance explained"
    x1 <- sprintf("PC1 (%.1f%% variance)", pct_var[1])
    y1 <- sprintf("PC2 (%.1f%% variance)", pct_var[2])
    x2 <- sprintf("PC1 (%.1f%%)", pct_var[1])
    y2 <- sprintf("PC3 (%.1f%%)", pct_var[3])
    x3 <- "Principal Component"
    y3 <- "Variance explained (%)"
    col_lab <- "ER status"
    cap1 <- sprintf("Outliers (Mahalanobis p<0.001): n=%d", n_outliers)
  }

  p1 <- ggplot(annotated, aes(x = PC1, y = PC2)) +
    {if (!is.null(er_col)) aes(color = as.factor(get(er_col))) else NULL} +
    geom_point(alpha = 0.5, size = 1.8) +
    {if (n_outliers > 0 && n_outliers <= 20) {
      geom_label_repel(
        data  = annotated[annotated$is_outlier, ],
        aes(label = sample_id),
        size  = 2.5,
        color = "red",
        max.overlaps = 20
      )
    }} +
    labs(title = t1, x = x1, y = y1, color = col_lab, caption = cap1) +
    theme_classic(base_size = 12)

  p2 <- ggplot(annotated, aes(x = PC1, y = PC3)) +
    {if (!is.null(er_col)) aes(color = as.factor(get(er_col))) else NULL} +
    geom_point(alpha = 0.5, size = 1.8) +
    labs(title = t2, x = x2, y = y2, color = col_lab) +
    theme_classic(base_size = 12)

  p3 <- ggplot(
    data.frame(pc = seq_len(10), pct = pct_var[seq_len(10)]),
    aes(x = pc, y = pct)
  ) +
    geom_col(fill = "#2980B9") +
    geom_line(color = "red", linewidth = 0.8) +
    geom_point(color = "red", size = 2) +
    labs(title = t3, x = x3, y = y3) +
    theme_classic(base_size = 12)

  list(scatter1 = p1, scatter2 = p2, scree = p3)
}

plots_en <- make_pca_plots("EN")
plots_pt <- make_pca_plots("PT")
p1 <- plots_en$scatter1; p2 <- plots_en$scatter2; p3 <- plots_en$scree

# ---------------------------------------------------------------------------
# 7) Save figures — split: scatter (2 panels) + scree (1 panel), EN + PT
# ---------------------------------------------------------------------------
for (d in c(PATHS$figures$supp_en_pdf, PATHS$figures$supp_en_png,
            PATHS$figures$supp_pt_pdf, PATHS$figures$supp_pt_png)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

old_warn <- getOption("warn"); options(warn = 0)

# --- Legacy 3-panel EN (backwards compat) ---
pdf(out_pdf, width = 14, height = 5)
gridExtra::grid.arrange(p1, p2, p3, ncol = 3)
dev.off()
out_png_legacy <- file.path(PATHS$figures$supp_en_png,
                            sprintf("FigS4_%s_PCA_forensics.png", COHORT))
png(out_png_legacy, width = 1400, height = 500, res = 100)
gridExtra::grid.arrange(p1, p2, p3, ncol = 3)
dev.off()

# --- EN: scatter (2-panel, larger) ---
en_scatter_pdf <- file.path(PATHS$figures$supp_en_pdf, "FigS4_METABRIC_PCA_scatter_EN.pdf")
en_scatter_png <- file.path(PATHS$figures$supp_en_png, "FigS4_METABRIC_PCA_scatter_EN.png")
cairo_pdf(en_scatter_pdf, width = 12, height = 6)
gridExtra::grid.arrange(plots_en$scatter1, plots_en$scatter2, ncol = 2)
dev.off()
png(en_scatter_png, width = 3600, height = 1800, res = 300)
gridExtra::grid.arrange(plots_en$scatter1, plots_en$scatter2, ncol = 2)
dev.off()

# --- EN: scree (standalone, larger) ---
en_scree_pdf <- file.path(PATHS$figures$supp_en_pdf, "FigS4_METABRIC_PCA_scree_EN.pdf")
en_scree_png <- file.path(PATHS$figures$supp_en_png, "FigS4_METABRIC_PCA_scree_EN.png")
cairo_pdf(en_scree_pdf, width = 7, height = 5)
print(plots_en$scree)
dev.off()
png(en_scree_png, width = 2100, height = 1500, res = 300)
print(plots_en$scree)
dev.off()

# --- PT: scatter (2-panel, larger) ---
pt_scatter_pdf <- file.path(PATHS$figures$supp_pt_pdf, "FigS4_METABRIC_PCA_scatter_PT.pdf")
pt_scatter_png <- file.path(PATHS$figures$supp_pt_png, "FigS4_METABRIC_PCA_scatter_PT.png")
cairo_pdf(pt_scatter_pdf, width = 12, height = 6)
gridExtra::grid.arrange(plots_pt$scatter1, plots_pt$scatter2, ncol = 2)
dev.off()
png(pt_scatter_png, width = 3600, height = 1800, res = 300)
gridExtra::grid.arrange(plots_pt$scatter1, plots_pt$scatter2, ncol = 2)
dev.off()

# --- PT: scree (standalone, larger) ---
pt_scree_pdf <- file.path(PATHS$figures$supp_pt_pdf, "FigS4_METABRIC_PCA_scree_PT.pdf")
pt_scree_png <- file.path(PATHS$figures$supp_pt_png, "FigS4_METABRIC_PCA_scree_PT.png")
cairo_pdf(pt_scree_pdf, width = 7, height = 5)
print(plots_pt$scree)
dev.off()
png(pt_scree_png, width = 2100, height = 1500, res = 300)
print(plots_pt$scree)
dev.off()

# --- PT: legacy combined (for backwards compat with existing QMD reference) ---
pt_combined_pdf <- file.path(PATHS$figures$supp_pt_pdf, "FigS4_METABRIC_PCA_forensics_PT.pdf")
pt_combined_png <- file.path(PATHS$figures$supp_pt_png, "FigS4_METABRIC_PCA_forensics_PT.png")
cairo_pdf(pt_combined_pdf, width = 14, height = 5)
gridExtra::grid.arrange(plots_pt$scatter1, plots_pt$scatter2, plots_pt$scree, ncol = 3)
dev.off()
png(pt_combined_png, width = 4200, height = 1500, res = 300)
gridExtra::grid.arrange(plots_pt$scatter1, plots_pt$scatter2, plots_pt$scree, ncol = 3)
dev.off()

options(warn = old_warn)

# Registry
for (f in c(out_pdf, en_scatter_pdf, pt_scatter_pdf, pt_scree_pdf)) {
  if (file.exists(f)) {
    h <- sha256_file(f)
    registry_append(COHORT, "figure_pca_forensics", f, h, "ok",
                    SCRIPT_NAME, file.info(f)$size / 1e6)
  }
}

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
