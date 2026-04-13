# =============================================================================
# SCRIPT: 30b_add_GSE1456_frozen.R
# PURPOSE: Add GSE1456 to frozen_zscore_sensitivity.csv using pre-existing
#          SCAN-B reference (avoids reloading 544 MB SCAN-B expression).
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "30b_add_GSE1456_frozen.R"

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
})

set.seed(FREEZE$seed_folds)
message(sprintf("[%s] Adding GSE1456 to frozen z-score table", SCRIPT_NAME))

# 1) Load CorePAM weights
weights_df    <- strict_csv(file.path(PATHS$results$corepam, "CorePAM_weights.csv"))
panel_genes   <- weights_df$gene
panel_weights <- setNames(weights_df$weight, weights_df$gene)
message(sprintf("[%s] CorePAM panel: %d genes", SCRIPT_NAME, length(panel_genes)))

# 2) Load SCAN-B reference (24 genes, pre-computed)
ref_df <- read.csv(file.path(PATHS$results$corepam, "SCANB_reference_meanSD.csv"),
                   stringsAsFactors = FALSE)
ref_means <- setNames(ref_df$mean_scanb, ref_df$gene)
ref_sds   <- setNames(ref_df$sd_scanb, ref_df$gene)
message(sprintf("[%s] SCAN-B reference loaded: %d genes", SCRIPT_NAME, nrow(ref_df)))

# 3) Load GSE1456 expression
coh_name <- "GSE1456"
expr_path <- file.path(proc_cohort(coh_name), "expression_genelevel_preZ.parquet")
expr_df   <- strict_parquet(expr_path)
gene_names <- expr_df[["gene"]]
s_ids      <- setdiff(names(expr_df), "gene")
expr_vals  <- t(as.matrix(expr_df[, s_ids]))
colnames(expr_vals) <- gene_names

# Panel genes present
genes_here <- intersect(panel_genes, gene_names)
n_present  <- length(genes_here)
message(sprintf("[%s] %s: %d/%d panel genes present", SCRIPT_NAME, coh_name, n_present, length(panel_genes)))

# 4) Compute frozen z-score
w_here  <- panel_weights[genes_here]
mu_ref  <- ref_means[genes_here]
sd_ref  <- ref_sds[genes_here]

frozen_z <- sweep(expr_vals[, genes_here, drop = FALSE], 2, mu_ref, "-")
frozen_z <- sweep(frozen_z, 2, sd_ref, "/")
frozen_z[is.na(frozen_z)] <- 0
score_z_frozen <- as.vector(frozen_z %*% w_here / sum(abs(w_here)))

# 5) Compute intra-cohort z-score
intra_z <- scale(expr_vals[, genes_here, drop = FALSE])
intra_z[is.na(intra_z)] <- 0
score_z_intra <- as.vector(intra_z %*% w_here / sum(abs(w_here)))

# 6) Load clinical data
ar_path <- file.path(proc_cohort(coh_name), "analysis_ready.parquet")
ar_df   <- strict_parquet(ar_path)

# Match samples
comp_df <- data.frame(
  sample_id      = s_ids,
  score_z_frozen = score_z_frozen,
  score_z_intra  = score_z_intra,
  stringsAsFactors = FALSE
)

merged <- data.frame()
for (jk in c("sample_id", "patient_id")) {
  if (jk %in% names(ar_df)) {
    tmp <- comp_df
    if (jk != "sample_id") tmp[[jk]] <- s_ids
    candidate <- merge(ar_df, tmp, by = jk)
    if (nrow(candidate) > nrow(merged)) merged <- candidate
  }
}
stopifnot(nrow(merged) > 0)

# Filter valid
time_col  <- "os_time_months"
event_col <- "os_event"
merged <- merged[!is.na(merged[[time_col]]) & merged[[time_col]] > 0 &
                   !is.na(merged[[event_col]]), ]
N <- nrow(merged)
n_events <- sum(merged[[event_col]])
message(sprintf("[%s] %s: N=%d, Events=%d", SCRIPT_NAME, coh_name, N, n_events))

# 7) Cox + C-index for frozen
cox_frozen <- coxph(Surv(merged[[time_col]], merged[[event_col]]) ~ merged$score_z_frozen)
sm_frozen  <- summary(cox_frozen)
hr_frozen     <- sm_frozen$conf.int[1, "exp(coef)"]
hr_frozen_lo  <- sm_frozen$conf.int[1, "lower .95"]
hr_frozen_hi  <- sm_frozen$conf.int[1, "upper .95"]
p_frozen      <- sm_frozen$coefficients[1, "Pr(>|z|)"]
c_frozen      <- sm_frozen$concordance[1]
c_frozen_se   <- sm_frozen$concordance[2]
c_frozen_lo   <- c_frozen - 1.96 * c_frozen_se
c_frozen_hi   <- c_frozen + 1.96 * c_frozen_se

# 8) Cox + C-index for intra
cox_intra <- coxph(Surv(merged[[time_col]], merged[[event_col]]) ~ merged$score_z_intra)
sm_intra  <- summary(cox_intra)
hr_intra     <- sm_intra$conf.int[1, "exp(coef)"]
hr_intra_lo  <- sm_intra$conf.int[1, "lower .95"]
hr_intra_hi  <- sm_intra$conf.int[1, "upper .95"]
p_intra      <- sm_intra$coefficients[1, "Pr(>|z|)"]
c_intra      <- sm_intra$concordance[1]
c_intra_se   <- sm_intra$concordance[2]
c_intra_lo   <- c_intra - 1.96 * c_intra_se
c_intra_hi   <- c_intra + 1.96 * c_intra_se

# 9) Correlation and delta-C
cor_val <- cor(merged$score_z_frozen, merged$score_z_intra, method = "pearson")
delta_c <- c_frozen - c_intra

message(sprintf("[%s] %s FROZEN: HR=%.4f (%.4f-%.4f) C=%.4f | INTRA: HR=%.4f (%.4f-%.4f) C=%.4f | dC=%.4f r=%.3f",
                SCRIPT_NAME, coh_name,
                hr_frozen, hr_frozen_lo, hr_frozen_hi, c_frozen,
                hr_intra, hr_intra_lo, hr_intra_hi, c_intra,
                delta_c, cor_val))

# 10) Append to CSV
new_row <- data.frame(
  cohort           = coh_name,
  platform         = "Microarray",
  n                = N,
  n_events         = n_events,
  genes_present    = n_present,
  hr_frozen        = hr_frozen,
  hr_frozen_lo95   = hr_frozen_lo,
  hr_frozen_hi95   = hr_frozen_hi,
  p_frozen         = p_frozen,
  c_frozen         = c_frozen,
  c_frozen_lo95    = c_frozen_lo,
  c_frozen_hi95    = c_frozen_hi,
  hr_intra         = hr_intra,
  hr_intra_lo95    = hr_intra_lo,
  hr_intra_hi95    = hr_intra_hi,
  p_intra          = p_intra,
  c_intra          = c_intra,
  c_intra_lo95     = c_intra_lo,
  c_intra_hi95     = c_intra_hi,
  delta_c_frozen_vs_intra = delta_c,
  cor_frozen_intra = cor_val,
  stringsAsFactors = FALSE
)

frozen_path <- file.path(PATHS$results$supp, "frozen_zscore_sensitivity.csv")
existing <- read.csv(frozen_path, stringsAsFactors = FALSE)
if (!"GSE1456" %in% existing$cohort) {
  updated <- rbind(existing, new_row)
  write.csv(updated, frozen_path, row.names = FALSE)
  message(sprintf("[%s] Appended GSE1456 to: %s", SCRIPT_NAME, frozen_path))
} else {
  message(sprintf("[%s] GSE1456 already in table, skipping", SCRIPT_NAME))
}

message(sprintf("[%s] COMPLETED", SCRIPT_NAME))
