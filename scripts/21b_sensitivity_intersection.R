# =============================================================================
# SCRIPT: 21b_sensitivity_intersection.R
# PURPOSE: Sensitivity analysis — recompute pCR OR using only the 21-gene
#          intersection present in ALL 4 primary pCR cohorts.
#
# Context: 3 genes (GPR160, CXXC5, KRT17) are absent in >= 1 cohort due to
# platform-level probe absence (structural missingness). Together they
# contribute 7.2% of Σ|w|.
#
# For GSE25066 and GSE20194: already use 21-gene subset (scores unchanged).
# For GSE32646 (23/24 → 21) and ISPY1 (24/24 → 21): recompute from expression.
#
# Output:
#   results/pcr/pcr_sensitivity_intersection.csv  — OR, AUC, CI per cohort
#   results/pcr/pcr_sensitivity_meta_intersection.csv — pooled meta OR
#
# PROJECT: Core-PAM (Memorial v6.1 / Freeze Core-PAM — pCR block)
# =============================================================================

source("scripts/00_setup.R")
source("scripts/03_utils_gene_mapping.R")

suppressPackageStartupMessages({
  library(arrow)
  library(GEOquery)
  library(Biobase)
  library(pROC)
  library(metafor)
})

SCRIPT_NAME <- "21b_sensitivity_intersection.R"

FORCE <- as.logical(Sys.getenv("FORCE_RERUN", "FALSE"))
out_path <- file.path("results/pcr", "pcr_sensitivity_intersection.csv")
if (!FORCE && file.exists(out_path)) {
  message(sprintf("[%s] Output exists — skipping. Set FORCE_RERUN=TRUE.", SCRIPT_NAME))
  quit(save = "no", status = 0)
}

message(sprintf("[%s] Sensitivity analysis: 21-gene intersection score", SCRIPT_NAME))

# --------------------------------------------------------------------------
# CorePAM weights — restrict to intersection (21 genes)
# --------------------------------------------------------------------------
weights_df    <- strict_csv(file.path(PATHS$results$corepam, "CorePAM_weights.csv"))
weights_df    <- weights_df[weights_df$weight != 0, ]
panel_genes   <- weights_df$gene
panel_weights <- setNames(weights_df$weight, weights_df$gene)

# 3 genes absent in >= 1 cohort (identified by check_pcr_gene_coverage.R)
GENES_EXCLUDED <- c("GPR160", "CXXC5", "KRT17")
inter_genes    <- setdiff(panel_genes, GENES_EXCLUDED)
inter_weights  <- panel_weights[inter_genes]
message(sprintf("[%s] Intersection: %d genes | excluded: %s | Σ|w|=%.4f (full: %.4f)",
                SCRIPT_NAME, length(inter_genes), paste(GENES_EXCLUDED, collapse=", "),
                sum(abs(inter_weights)), sum(abs(panel_weights))))

# --------------------------------------------------------------------------
# Helper: compute intersection score from series matrix
# --------------------------------------------------------------------------
compute_inter_score_from_geo <- function(cohort, geo_acc) {
  raw_dir <- raw_pcr_cohort(cohort)
  sm_path <- file.path(raw_dir, sprintf("%s_series_matrix.txt.gz", geo_acc))
  if (!file.exists(sm_path)) {
    stop(sprintf("[%s] Series matrix not found: %s", SCRIPT_NAME, sm_path))
  }

  old_warn <- getOption("warn"); options(warn = 0)
  gse <- GEOquery::getGEO(filename = sm_path)
  options(warn = old_warn)

  if (is.list(gse)) {
    sizes <- sapply(gse, function(x) ncol(Biobase::exprs(x)))
    gse   <- gse[[which.max(sizes)]]
  }

  expr_mat  <- Biobase::exprs(gse)
  fdata     <- Biobase::fData(gse)
  plat      <- annotation(gse)
  message(sprintf("[%s] %s: platform=%s, probes=%d, samples=%d",
                  SCRIPT_NAME, cohort, plat, nrow(expr_mat), ncol(expr_mat)))

  # Gene symbol mapping
  gene_symbols <- NULL
  for (cand in c("Gene Symbol", "GENE_SYMBOL", "gene_symbol", "GeneName",
                 "GENE_NAME", "SystematicName")) {
    if (cand %in% names(fdata)) {
      gs_raw <- as.character(fdata[[cand]])
      valid  <- gs_raw[nchar(trimws(gs_raw)) > 0 & !is.na(gs_raw) &
                       gs_raw != "---" & gs_raw != ""]
      if (length(valid) > 100) {
        gene_symbols <- trimws(gs_raw)
        break
      }
    }
  }
  if (is.null(gene_symbols)) {
    stop(sprintf("[%s] %s: no gene symbol column in fData", SCRIPT_NAME, cohort))
  }

  valid_idx    <- which(!is.na(gene_symbols) & nchar(gene_symbols) > 0 &
                        gene_symbols != "---" & gene_symbols != "")
  mapped_probes <- rownames(expr_mat)[valid_idx]
  gene_labels   <- gene_symbols[valid_idx]

  # Remove duplicate probe IDs
  dup_probes    <- duplicated(mapped_probes)
  mapped_probes <- mapped_probes[!dup_probes]
  gene_labels   <- gene_labels[!dup_probes]

  expr_sub       <- expr_mat[mapped_probes, , drop = FALSE]
  expr_collapsed <- collapse_probes_by_variance(expr_sub, gene_labels)

  # Z-score
  z_t        <- t(expr_collapsed)
  gene_means <- colMeans(z_t, na.rm = TRUE)
  old_warn   <- getOption("warn"); options(warn = 0)
  gene_sds   <- apply(z_t, 2, sd, na.rm = TRUE)
  options(warn = old_warn)
  zero_sd    <- gene_sds == 0 | is.na(gene_sds)
  gene_sds[zero_sd]  <- 1
  z_t[, zero_sd]     <- 0
  z_mat <- sweep(sweep(z_t, 2, gene_means, "-"), 2, gene_sds, "/")

  # Intersection score
  genes_use    <- intersect(inter_genes, colnames(z_mat))
  w_use        <- inter_weights[genes_use]
  denom        <- sum(abs(w_use))
  score_inter  <- as.vector(z_mat[, genes_use] %*% w_use) / denom
  old_warn     <- getOption("warn"); options(warn = 0)
  score_z_inter <- as.vector(scale(score_inter))
  options(warn = old_warn)
  message(sprintf("[%s] %s: %d/%d intersection genes used",
                  SCRIPT_NAME, cohort, length(genes_use), length(inter_genes)))

  data.frame(
    sample_id     = colnames(expr_collapsed),
    score_inter   = score_inter,
    score_z_inter = score_z_inter,
    stringsAsFactors = FALSE
  )
}

# --------------------------------------------------------------------------
# Main loop: for each cohort, get intersection score and run logistic
# --------------------------------------------------------------------------
PCR_COHORTS <- list(
  list(cohort = "GSE25066", geo_acc = "GSE25066"),
  list(cohort = "GSE20194", geo_acc = "GSE20194"),
  list(cohort = "GSE32646", geo_acc = "GSE32646"),
  list(cohort = "ISPY1",    geo_acc = "GSE22226")
)

sens_results <- list()

for (ch in PCR_COHORTS) {
  cohort  <- ch$cohort
  geo_acc <- ch$geo_acc

  message(sprintf("\n[%s] === %s ===", SCRIPT_NAME, cohort))

  # Load pCR data from analysis_ready.parquet
  pq_path <- file.path(proc_pcr_cohort(cohort), "analysis_ready.parquet")
  old_warn <- getOption("warn"); options(warn = 0)
  pdata    <- arrow::read_parquet(pq_path)
  options(warn = old_warn)
  message(sprintf("[%s] %s: N=%d, pCR=%d", SCRIPT_NAME, cohort,
                  nrow(pdata), sum(pdata$pcr == 1L, na.rm = TRUE)))

  # Compute intersection score
  inter_scores <- tryCatch(
    compute_inter_score_from_geo(cohort, geo_acc),
    error = function(e) {
      message(sprintf("[%s] %s: ERROR computing intersection score: %s",
                      SCRIPT_NAME, cohort, e$message))
      NULL
    }
  )
  if (is.null(inter_scores)) next

  # Merge
  merged <- merge(pdata, inter_scores, by = "sample_id", sort = FALSE)
  merged <- merged[!is.na(merged$pcr), ]

  # Logistic regression: pcr ~ score_z_inter
  fit  <- glm(pcr ~ score_z_inter, data = merged, family = binomial())
  coef_fit <- summary(fit)$coefficients
  or_inter  <- exp(coef_fit["score_z_inter", "Estimate"])
  se_log    <- coef_fit["score_z_inter", "Std. Error"]
  p_inter   <- coef_fit["score_z_inter", "Pr(>|z|)"]
  ci_lo     <- exp(coef_fit["score_z_inter", "Estimate"] - 1.96 * se_log)
  ci_hi     <- exp(coef_fit["score_z_inter", "Estimate"] + 1.96 * se_log)

  # AUC
  old_warn <- getOption("warn"); options(warn = 0)
  roc_obj  <- pROC::roc(merged$pcr, merged$score_z_inter, quiet = TRUE)
  auc_val  <- as.numeric(pROC::auc(roc_obj))
  options(warn = old_warn)

  sens_results[[cohort]] <- data.frame(
    cohort      = cohort,
    n_samples   = nrow(merged),
    pcr_rate    = mean(merged$pcr, na.rm = TRUE),
    n_genes     = length(inter_genes),
    or_inter    = or_inter,
    ci_lo       = ci_lo,
    ci_hi       = ci_hi,
    p_value     = p_inter,
    log_or      = log(or_inter),
    se_log_or   = se_log,
    auc_inter   = auc_val,
    stringsAsFactors = FALSE
  )

  message(sprintf("[%s] %s: OR=%.3f (%.3f-%.3f) p=%.4f AUC=%.3f",
                  SCRIPT_NAME, cohort, or_inter, ci_lo, ci_hi, p_inter, auc_val))
}

if (length(sens_results) == 0) stop(sprintf("[%s] No sensitivity results", SCRIPT_NAME))

sens_df <- do.call(rbind, sens_results)
rownames(sens_df) <- NULL

# --------------------------------------------------------------------------
# Pooled meta-analysis (intersection score)
# --------------------------------------------------------------------------
old_warn <- getOption("warn"); options(warn = 0)
meta_inter <- metafor::rma(yi = log_or, sei = se_log_or, data = sens_df,
                           method = "DL", test = "z")
options(warn = old_warn)
pooled_or_inter  <- exp(meta_inter$beta[1])
pooled_ci_lo     <- exp(meta_inter$ci.lb)
pooled_ci_hi     <- exp(meta_inter$ci.ub)
pooled_p         <- meta_inter$pval
i2_inter         <- max(0, meta_inter$I2)

message(sprintf("\n[%s] POOLED (intersection) OR=%.3f (%.3f-%.3f) p=%.2e I²=%.0f%%",
                SCRIPT_NAME, pooled_or_inter, pooled_ci_lo, pooled_ci_hi,
                pooled_p, i2_inter))

# --------------------------------------------------------------------------
# Save (before optional comparison to avoid aborting on column name issues)
# --------------------------------------------------------------------------
dir.create("results/pcr", showWarnings = FALSE, recursive = TRUE)
write.csv(sens_df, out_path, row.names = FALSE)

meta_sum_df <- data.frame(
  analysis   = "intersection_21genes",
  or_pooled  = pooled_or_inter,
  ci_lo      = pooled_ci_lo,
  ci_hi      = pooled_ci_hi,
  p_value    = pooled_p,
  i2         = i2_inter,
  n_cohorts  = nrow(sens_df),
  genes_used = length(inter_genes),
  stringsAsFactors = FALSE
)
meta_path <- file.path("results/pcr", "pcr_sensitivity_meta_intersection.csv")
write.csv(meta_sum_df, meta_path, row.names = FALSE)

h <- sha256_file(out_path)
registry_append("PCR_SENSITIVITY_INTER", "pcr_sensitivity_intersection",
                out_path, h, "ok", SCRIPT_NAME,
                file.info(out_path)$size / 1e6,
                list(genes = length(inter_genes), cohorts = nrow(sens_df)))

message(sprintf("\n[%s] DONE | Sensitivity table: %s | Meta: %s",
                SCRIPT_NAME, out_path, meta_path))

# --------------------------------------------------------------------------
# Optional comparison with original results
# --------------------------------------------------------------------------
orig_df <- tryCatch({
  old_warn <- getOption("warn"); options(warn = 0)
  df <- read.csv(file.path("results/pcr", "pcr_validation_table.csv"),
                 stringsAsFactors = FALSE)
  options(warn = old_warn)
  df
}, error = function(e) NULL)

message("\n========== COMPARISON: original vs intersection ==========")
if (!is.null(orig_df)) {
  or_col    <- grep("^or_uni$|^or$", names(orig_df), value = TRUE)[1]
  ci_lo_col <- grep("lo95|ci_lo|lo_wald", names(orig_df), value = TRUE)[1]
  ci_hi_col <- grep("hi95|ci_hi|hi_wald", names(orig_df), value = TRUE)[1]
  for (i in seq_len(nrow(sens_df))) {
    coh      <- sens_df$cohort[i]
    orig_row <- orig_df[orig_df$cohort == coh, ]
    if (nrow(orig_row) > 0 && !is.na(or_col)) {
      message(sprintf("  %s: orig OR=%.3f (%.3f-%.3f) | inter OR=%.3f (%.3f-%.3f)",
                      coh,
                      orig_row[[or_col]][1],
                      if (!is.na(ci_lo_col)) orig_row[[ci_lo_col]][1] else NA,
                      if (!is.na(ci_hi_col)) orig_row[[ci_hi_col]][1] else NA,
                      sens_df$or_inter[i], sens_df$ci_lo[i], sens_df$ci_hi[i]))
    }
  }
}
