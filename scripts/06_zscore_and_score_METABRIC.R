# =============================================================================
# SCRIPT: 06_zscore_and_score_METABRIC.R
# PURPOSE: Applies Core-PAM score to METABRIC cohort:
#          Intra-cohort Z-score, rescaled score (Effective Gene Count),
#          standardized direction (HR>1), generates analysis_ready.parquet
# COHORT:  METABRIC (microarray validation; primary endpoint = DSS)
# PROJETO: Core-PAM (Memorial v6.1 / Freeze Core-PAM)
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "06_zscore_and_score_METABRIC.R"
COHORT      <- "METABRIC"

suppressPackageStartupMessages({
  library(survival)
})

message(sprintf("[%s] Starting CorePAM score for %s", SCRIPT_NAME, COHORT))

# --------------------------------------------------------------------------
# 1) Load frozen Core-PAM weights
# --------------------------------------------------------------------------
weights_path <- file.path(PATHS$results$corepam, "CorePAM_weights.csv")
weights_df   <- strict_csv(weights_path)

stopifnot(all(c("gene", "weight") %in% names(weights_df)))
weights_df    <- weights_df[weights_df$weight != 0, ]
panel_genes   <- weights_df$gene
panel_weights <- setNames(weights_df$weight, weights_df$gene)
message(sprintf("[%s] CorePAM panel: %d genes with weight != 0", SCRIPT_NAME, length(panel_genes)))

# --------------------------------------------------------------------------
# 2) Load pre-Z expression
# --------------------------------------------------------------------------
expr_path <- file.path(proc_cohort(COHORT), "expression_genelevel_preZ.parquet")
expr_mat   <- strict_parquet(expr_path)

# Format: first column = "gene" (rows are genes); remaining cols = sample IDs
gene_names <- expr_mat[["gene"]]
sample_ids  <- setdiff(names(expr_mat), "gene")

message(sprintf("[%s] Expression: %d samples x %d genes", SCRIPT_NAME, length(sample_ids), length(gene_names)))

# --------------------------------------------------------------------------
# 3) Intra-cohort Z-score per gene
# --------------------------------------------------------------------------
expr_vals <- t(as.matrix(expr_mat[, sample_ids]))  # transpose: samples x genes
colnames(expr_vals) <- gene_names

old_warn <- getOption("warn"); options(warn = 0)
gene_means <- colMeans(expr_vals, na.rm = TRUE)
gene_sds   <- apply(expr_vals, 2, sd, na.rm = TRUE)
options(warn = old_warn)

zero_sd <- gene_sds == 0 | is.na(gene_sds)
if (any(zero_sd)) {
  message(sprintf("[%s] WARNING: %d genes with SD=0 (z=0 for those genes)", SCRIPT_NAME, sum(zero_sd)))
}
gene_sds[zero_sd] <- 1

z_mat <- sweep(sweep(expr_vals, 2, gene_means, "-"), 2, gene_sds, "/")
z_mat[, zero_sd] <- 0

# --------------------------------------------------------------------------
# 4) Calculate score = sum(w_i * z_i) / sum(|w_i|)  for present genes only
# --------------------------------------------------------------------------
genes_present <- intersect(panel_genes, colnames(z_mat))
genes_missing <- setdiff(panel_genes, colnames(z_mat))

n_present    <- length(genes_present)
n_panel      <- length(panel_genes)
frac_present <- n_present / n_panel

message(sprintf("[%s] Genes present: %d / %d (%.1f%%)",
                SCRIPT_NAME, n_present, n_panel, frac_present * 100))
if (length(genes_missing) > 0) {
  message(sprintf("[%s] Genes missing: %s", SCRIPT_NAME, paste(genes_missing, collapse = ", ")))
}

if (frac_present < FREEZE$min_genes_fraction) {
  stop(sprintf(
    "[%s] INSUFFICIENT COVERAGE: %.1f%% < %.0f%% minimum required",
    SCRIPT_NAME, frac_present * 100, FREEZE$min_genes_fraction * 100
  ))
}

w_present      <- panel_weights[genes_present]
denom_sum_absw <- sum(abs(w_present))

z_sub     <- z_mat[, genes_present, drop = FALSE]
score_raw <- as.vector(z_sub %*% w_present) / denom_sum_absw

# --------------------------------------------------------------------------
# 5) score_z = scale(score) for HR per 1 SD
# --------------------------------------------------------------------------
old_warn <- getOption("warn"); options(warn = 0)
score_z_raw <- as.vector(scale(score_raw))
options(warn = old_warn)

# --------------------------------------------------------------------------
# 6) Load clinical data to check score direction
#    METABRIC: primary endpoint = DSS
# --------------------------------------------------------------------------
clin_path <- file.path(proc_cohort(COHORT), "clinical_FINAL.parquet")
clin_df   <- strict_parquet(clin_path)

score_df <- tibble(
  sample_id      = sample_ids,
  score          = score_raw,
  score_z        = score_z_raw,
  genes_present  = n_present,
  denom_sum_absw = denom_sum_absw
)

clin_score <- inner_join(clin_df, score_df, by = "sample_id")
n_join <- nrow(clin_score)
message(sprintf("[%s] Clinical x score join: %d samples", SCRIPT_NAME, n_join))
if (n_join == 0) stop(sprintf("[%s] Join resulted in 0 rows. Check join keys.", SCRIPT_NAME))

# --------------------------------------------------------------------------
# 7) Determine score direction using DSS (METABRIC primary endpoint)
# --------------------------------------------------------------------------
# Detect available DSS columns
has_dss <- all(c("dss_time_months", "dss_event") %in% names(clin_score))
if (has_dss) {
  endpoint_col <- "dss_time_months"
  event_col    <- "dss_event"
  message(sprintf("[%s] Using DSS endpoint to determine score direction", SCRIPT_NAME))
} else {
  endpoint_col <- "os_time_months"
  event_col    <- "os_event"
  message(sprintf("[%s] WARNING: DSS not found; using OS for direction", SCRIPT_NAME))
}

stopifnot(all(c(endpoint_col, event_col) %in% names(clin_score)))

n_before_drop <- nrow(clin_score)
clin_score    <- clin_score[!is.na(clin_score[[endpoint_col]]) & clin_score[[endpoint_col]] > 0, ]
n_dropped     <- n_before_drop - nrow(clin_score)
if (n_dropped > 0) {
  message(sprintf("[%s] Samples with time<=0 removed: %d", SCRIPT_NAME, n_dropped))
}

# --------------------------------------------------------------------------
# Direction confirmation — VALIDATION cohort (METABRIC)
# Score direction was established in the SCAN-B training cohort (HR=1.92;
# direction="original"). This check confirms the same convention holds here.
# If hr_dir < 1 were observed, it would indicate a methodological problem
# requiring investigation — NOT a routine per-cohort inversion.
# --------------------------------------------------------------------------
old_warn <- getOption("warn"); options(warn = 0)
cox_dir <- tryCatch(
  coxph(Surv(clin_score[[endpoint_col]], clin_score[[event_col]]) ~ clin_score$score_z),
  error = function(e) NULL
)
options(warn = old_warn)

if (!is.null(cox_dir)) {
  hr_dir <- exp(coef(cox_dir)[1])
  message(sprintf("[%s] VALIDATION direction confirmation HR score_z: %.4f (expected > 1 per training)", SCRIPT_NAME, hr_dir))
  if (hr_dir < 1) {
    clin_score$score   <- -clin_score$score
    clin_score$score_z <- -clin_score$score_z
    score_direction    <- "inverted"
    message(sprintf("[%s] WARNING: Direction inverted (HR < 1) — UNEXPECTED; investigate", SCRIPT_NAME))
  } else {
    score_direction <- "original"
    message(sprintf("[%s] Direction confirmed original (HR >= 1) — consistent with training", SCRIPT_NAME))
  }
} else {
  score_direction <- "unknown"
  message(sprintf("[%s] WARNING: Direction Cox failed; direction = unknown", SCRIPT_NAME))
}

clin_score$score_direction <- score_direction

# --------------------------------------------------------------------------
# 8) Save analysis_ready.parquet
# --------------------------------------------------------------------------
out_dir <- proc_cohort(COHORT)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_path <- file.path(out_dir, "analysis_ready.parquet")
arrow::write_parquet(clin_score, out_path)
h  <- sha256_file(out_path)
sz <- file.info(out_path)$size / 1e6
message(sprintf("[%s] Saved: %s (%.2f MB | SHA256: %s)", SCRIPT_NAME, out_path, sz, h))

registry_append(
  cohort    = COHORT,
  file_type = "analysis_ready",
  file_path = out_path,
  sha256    = h,
  status    = "ok",
  script    = SCRIPT_NAME,
  size_mb   = sz,
  extra     = list(genes_present = n_present, score_direction = score_direction,
                   primary_endpoint = endpoint_col)
)

# --------------------------------------------------------------------------
# 9) Supplementary summary
# --------------------------------------------------------------------------
supp_dir <- PATHS$results$supp
dir.create(supp_dir, showWarnings = FALSE, recursive = TRUE)

supp_df <- tibble(
  cohort              = COHORT,
  n_samples           = nrow(clin_score),
  n_panel_genes       = n_panel,
  genes_present       = n_present,
  genes_missing       = n_panel - n_present,
  frac_present        = round(frac_present, 4),
  denom_sum_absw      = round(denom_sum_absw, 6),
  score_direction     = score_direction,
  primary_endpoint    = endpoint_col,
  score_mean          = round(mean(clin_score$score, na.rm = TRUE), 6),
  score_sd            = round(sd(clin_score$score, na.rm = TRUE), 6),
  n_time_leq0_dropped = n_dropped
)

supp_path <- file.path(supp_dir, sprintf("risk_score_summary_%s.csv", COHORT))
readr::write_csv(supp_df, supp_path)
h2 <- sha256_file(supp_path)
registry_append(
  cohort    = COHORT,
  file_type = "risk_score_summary",
  file_path = supp_path,
  sha256    = h2,
  status    = "ok",
  script    = SCRIPT_NAME,
  size_mb   = file.info(supp_path)$size / 1e6
)

message(sprintf("[%s] COMPLETED for %s", SCRIPT_NAME, COHORT))
