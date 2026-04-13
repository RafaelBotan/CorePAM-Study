# =============================================================================
# SCRIPT: 30_frozen_zscore_sensitivity.R
# PURPOSE: Frozen z-score sensitivity analysis — freezes mean/SD from SCAN-B
#          training cohort and applies to all validation cohorts (TRUE single-
#          sample scoring, NO final scale()). Compares HR + C-index frozen vs
#          intra-cohort z-score on the same raw-score scale.
#          Addresses reviewer concern: "z-score should use training reference."
# NOTE:    No scale() is applied after weighted sum — this ensures each
#          patient's score depends ONLY on their expression + frozen SCAN-B
#          reference. C-index comparison is rank-invariant; HR units are
#          per-unit of raw weighted score (not per-SD).
#          TABLE S4 reports C-index and correlation (r) as primary metrics;
#          raw-scale HRs are retained in CSV but omitted from the manuscript
#          table to avoid confusion with per-1-SD HRs in the main Table 3.
# OUTPUT:
#   results/supp/frozen_zscore_sensitivity.csv
#   results/corepam/SCANB_reference_meanSD.csv
# PROJECT: Core-PAM (Memorial CorePAM — Major Revision)
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "30_frozen_zscore_sensitivity.R"

suppressPackageStartupMessages({
  library(survival)
  library(arrow)
})

set.seed(FREEZE$seed_folds)
message(sprintf("[%s] ========== FROZEN Z-SCORE SENSITIVITY ==========", SCRIPT_NAME))

# --------------------------------------------------------------------------
# 1) Load frozen CorePAM weights
# --------------------------------------------------------------------------
weights_path <- file.path(PATHS$results$corepam, "CorePAM_weights.csv")
weights_df   <- strict_csv(weights_path)
weights_df   <- weights_df[weights_df$weight != 0, ]
panel_genes  <- weights_df$gene
panel_weights <- setNames(weights_df$weight, weights_df$gene)
message(sprintf("[%s] CorePAM panel: %d genes", SCRIPT_NAME, length(panel_genes)))

# --------------------------------------------------------------------------
# 2) Compute reference mean/SD from SCAN-B pre-Z expression
# --------------------------------------------------------------------------
scanb_expr_path <- file.path(proc_cohort("SCANB"), "expression_genelevel_preZ.parquet")
scanb_expr      <- strict_parquet(scanb_expr_path)

gene_names <- scanb_expr[["gene"]]
sample_ids <- setdiff(names(scanb_expr), "gene")

expr_vals <- t(as.matrix(scanb_expr[, sample_ids]))
colnames(expr_vals) <- gene_names

old_warn <- getOption("warn"); options(warn = 0)
ref_means <- colMeans(expr_vals, na.rm = TRUE)
ref_sds   <- apply(expr_vals, 2, sd, na.rm = TRUE)
options(warn = old_warn)

# Save reference statistics (only for CorePAM panel genes present in SCAN-B)
genes_in_scanb <- intersect(panel_genes, gene_names)
ref_df <- data.frame(
  gene    = genes_in_scanb,
  mean_scanb = ref_means[genes_in_scanb],
  sd_scanb   = ref_sds[genes_in_scanb],
  stringsAsFactors = FALSE
)
ref_path <- file.path(PATHS$results$corepam, "SCANB_reference_meanSD.csv")
write.csv(ref_df, ref_path, row.names = FALSE)
message(sprintf("[%s] Reference mean/SD saved: %s (%d genes)", SCRIPT_NAME, ref_path, nrow(ref_df)))

# --------------------------------------------------------------------------
# 3) Define cohorts + endpoints
# --------------------------------------------------------------------------
COHORTS <- list(
  SCANB     = list(time = "os_time_months",  event = "os_event"),
  TCGA_BRCA = list(time = "os_time_months",  event = "os_event"),
  METABRIC  = list(time = "dss_time_months", event = "dss_event"),
  GSE20685  = list(time = "os_time_months",  event = "os_event"),
  GSE1456   = list(time = "os_time_months",  event = "os_event")
)

# --------------------------------------------------------------------------
# 4) For each cohort: compute frozen z-score, compare vs intra-cohort
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

results <- list()

for (coh_name in names(COHORTS)) {
  message(sprintf("[%s] Processing %s...", SCRIPT_NAME, coh_name))
  ep <- COHORTS[[coh_name]]

  # Load pre-Z expression
  expr_path <- file.path(proc_cohort(coh_name), "expression_genelevel_preZ.parquet")
  if (!file.exists(expr_path)) {
    message(sprintf("[%s] %s: expression_genelevel_preZ.parquet not found, skipping", SCRIPT_NAME, coh_name))
    next
  }
  expr_raw <- strict_parquet(expr_path)
  g_names  <- expr_raw[["gene"]]
  s_ids    <- setdiff(names(expr_raw), "gene")
  e_vals   <- t(as.matrix(expr_raw[, s_ids]))
  colnames(e_vals) <- g_names

  # Genes present in both this cohort and reference
  genes_present <- intersect(panel_genes, intersect(g_names, ref_df$gene))
  genes_missing <- setdiff(panel_genes, genes_present)
  n_present <- length(genes_present)
  frac_present <- n_present / length(panel_genes)

  if (frac_present < FREEZE$min_genes_fraction) {
    message(sprintf("[%s] %s: insufficient gene coverage (%.1f%%), skipping",
                    SCRIPT_NAME, coh_name, frac_present * 100))
    next
  }

  # --- Frozen z-score: use SCAN-B reference mean/SD ---
  ref_m <- setNames(ref_df$mean_scanb, ref_df$gene)
  ref_s <- setNames(ref_df$sd_scanb,   ref_df$gene)
  # Protect against zero SD
  ref_s[ref_s == 0 | is.na(ref_s)] <- 1

  z_frozen <- sweep(e_vals[, genes_present, drop = FALSE], 2, ref_m[genes_present], "-")
  z_frozen <- sweep(z_frozen, 2, ref_s[genes_present], "/")

  w_present     <- panel_weights[genes_present]
  denom_absw    <- sum(abs(w_present))
  score_frozen  <- as.vector(z_frozen %*% w_present) / denom_absw

  # TRUE single-sample: NO scale() — each patient's score depends only on
  # their expression + frozen SCAN-B reference parameters. In clinical
  # deployment there is no cohort to scale against.
  score_z_frozen <- score_frozen

  # --- Intra-cohort z-score (original approach) ---
  old_warn <- getOption("warn"); options(warn = 0)
  intra_means <- colMeans(e_vals, na.rm = TRUE)
  intra_sds   <- apply(e_vals, 2, sd, na.rm = TRUE)
  options(warn = old_warn)
  intra_sds[intra_sds == 0 | is.na(intra_sds)] <- 1

  z_intra <- sweep(e_vals[, genes_present, drop = FALSE], 2, intra_means[genes_present], "-")
  z_intra <- sweep(z_intra, 2, intra_sds[genes_present], "/")
  score_intra <- as.vector(z_intra %*% w_present) / denom_absw
  # No scale() here either — fair comparison on same raw-score scale
  score_z_intra <- score_intra

  # --- Load clinical data ---
  ar_path <- file.path(proc_cohort(coh_name), "analysis_ready.parquet")
  ar_df   <- strict_parquet(ar_path)

  # Build comparison data frame
  comp_df <- data.frame(
    sample_id      = s_ids,
    score_z_frozen = score_z_frozen,
    score_z_intra  = score_z_intra,
    stringsAsFactors = FALSE
  )

  # Try sample_id first, then patient_id. Pick whichever yields more matches.
  merged <- data.frame()
  for (jk in c("sample_id", "patient_id")) {
    if (jk %in% names(ar_df)) {
      tmp <- comp_df
      if (jk != "sample_id") tmp[[jk]] <- s_ids
      candidate <- merge(ar_df, tmp, by = jk)
      if (nrow(candidate) > nrow(merged)) merged <- candidate
    }
  }
  if (nrow(merged) == 0) {
    message(sprintf("[%s] %s: no matching join key, skipping", SCRIPT_NAME, coh_name))
    next
  }

  time_col  <- ep$time
  event_col <- ep$event
  if (!all(c(time_col, event_col) %in% names(merged))) {
    message(sprintf("[%s] %s: endpoint columns missing, skipping", SCRIPT_NAME, coh_name))
    next
  }

  merged <- merged[!is.na(merged[[time_col]]) & merged[[time_col]] > 0 &
                     !is.na(merged[[event_col]]), ]
  n_analysis <- nrow(merged)

  if (n_analysis < 30) {
    message(sprintf("[%s] %s: n=%d too small, skipping", SCRIPT_NAME, coh_name, n_analysis))
    next
  }

  # --- Cox HR + C-index for frozen ---
  old_warn <- getOption("warn"); options(warn = 0)
  cox_frozen <- tryCatch(
    coxph(Surv(merged[[time_col]], merged[[event_col]]) ~ score_z_frozen, data = merged),
    error = function(e) NULL
  )
  options(warn = old_warn)

  hr_frozen <- lo_frozen <- hi_frozen <- p_frozen <- NA_real_
  if (!is.null(cox_frozen)) {
    sm <- summary(cox_frozen)
    hr_frozen <- sm$conf.int[1, "exp(coef)"]
    lo_frozen <- sm$conf.int[1, "lower .95"]
    hi_frozen <- sm$conf.int[1, "upper .95"]
    p_frozen  <- sm$coefficients[1, "Pr(>|z|)"]
  }

  boot_frozen <- bootstrap_cindex(
    merged[[time_col]], merged[[event_col]], merged$score_z_frozen,
    n_boot = as.integer(FREEZE$bootstrap_n),
    seed   = as.integer(FREEZE$seed_folds)
  )

  # --- Cox HR + C-index for intra-cohort ---
  old_warn <- getOption("warn"); options(warn = 0)
  cox_intra <- tryCatch(
    coxph(Surv(merged[[time_col]], merged[[event_col]]) ~ score_z_intra, data = merged),
    error = function(e) NULL
  )
  options(warn = old_warn)

  hr_intra <- lo_intra <- hi_intra <- p_intra <- NA_real_
  if (!is.null(cox_intra)) {
    sm2 <- summary(cox_intra)
    hr_intra <- sm2$conf.int[1, "exp(coef)"]
    lo_intra <- sm2$conf.int[1, "lower .95"]
    hi_intra <- sm2$conf.int[1, "upper .95"]
    p_intra  <- sm2$coefficients[1, "Pr(>|z|)"]
  }

  boot_intra <- bootstrap_cindex(
    merged[[time_col]], merged[[event_col]], merged$score_z_intra,
    n_boot = as.integer(FREEZE$bootstrap_n),
    seed   = as.integer(FREEZE$seed_folds)
  )

  # --- Correlation between frozen and intra-cohort scores ---
  cor_scores <- cor(merged$score_z_frozen, merged$score_z_intra, use = "complete.obs")

  # --- Store results ---
  results[[coh_name]] <- data.frame(
    cohort            = coh_name,
    platform          = ifelse(coh_name %in% c("SCANB", "TCGA_BRCA"), "RNA-seq", "Microarray"),
    n                 = n_analysis,
    n_events          = sum(merged[[event_col]]),
    genes_present     = n_present,
    # Frozen z-score results
    hr_frozen         = round(hr_frozen, 4),
    hr_frozen_lo95    = round(lo_frozen, 4),
    hr_frozen_hi95    = round(hi_frozen, 4),
    p_frozen          = signif(p_frozen, 4),
    c_frozen          = round(boot_frozen$c_index, 4),
    c_frozen_lo95     = round(boot_frozen$ci_low, 4),
    c_frozen_hi95     = round(boot_frozen$ci_high, 4),
    # Intra-cohort z-score results (reference)
    hr_intra          = round(hr_intra, 4),
    hr_intra_lo95     = round(lo_intra, 4),
    hr_intra_hi95     = round(hi_intra, 4),
    p_intra           = signif(p_intra, 4),
    c_intra           = round(boot_intra$c_index, 4),
    c_intra_lo95      = round(boot_intra$ci_low, 4),
    c_intra_hi95      = round(boot_intra$ci_high, 4),
    # Comparison
    delta_c_frozen_vs_intra = round(boot_frozen$c_index - boot_intra$c_index, 4),
    cor_frozen_intra  = round(cor_scores, 4),
    stringsAsFactors  = FALSE
  )

  message(sprintf("[%s] %s: HR_frozen=%.3f (%.3f-%.3f) C=%.4f | HR_intra=%.3f C=%.4f | dC=%.4f | r=%.3f",
                  SCRIPT_NAME, coh_name,
                  hr_frozen, lo_frozen, hi_frozen, boot_frozen$c_index,
                  hr_intra, boot_intra$c_index,
                  boot_frozen$c_index - boot_intra$c_index,
                  cor_scores))
}

# --------------------------------------------------------------------------
# 5) Save results
# --------------------------------------------------------------------------
if (length(results) == 0) {
  stop(sprintf("[%s] No results generated.", SCRIPT_NAME))
}

results_df <- do.call(rbind, results)
out_path   <- file.path(PATHS$results$supp, "frozen_zscore_sensitivity.csv")
write.csv(results_df, out_path, row.names = FALSE)
h <- sha256_file(out_path)
registry_append("ALL", "frozen_zscore_sensitivity", out_path, h, "ok", SCRIPT_NAME,
                file.info(out_path)$size / 1e6)

message(sprintf("[%s] Results saved: %s", SCRIPT_NAME, out_path))
message(sprintf("[%s] ========== DONE ==========", SCRIPT_NAME))
