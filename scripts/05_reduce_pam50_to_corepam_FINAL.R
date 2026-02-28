# =============================================================================
# SCRIPT: 05_reduce_pam50_to_corepam_FINAL.R
# PURPOSE: Core-PAM derivation — smallest PAM50 subset non-inferior to the
#          full PAM50 model by OOF C-index (Harrell) on SCAN-B.
# PROJETO: Core-PAM (Memorial v6.1 §5)
#
# ALGORITHM (frozen):
#   1. Prepare X (50 PAM50 genes, intra-SCANB Z-score) + y (Surv OS)
#   2. Fixed stratified K=10 folds (frozen seed)
#   3. cv.glmnet (alpha=0.5, keep=TRUE) → OOF linear predictors per lambda
#   4. For each point on path (lambda, df):
#        Cadj = max(Craw, 1-Craw)  — orientation-free
#   5. Cmax = max(Cadj); threshold = Cmax - delta_c
#   6. Select smallest df with Cadj >= threshold
#      Tiebreak: largest lambda (most parsimonious) within same df
#   7. Refit on full data: glmnet.fit at chosen lambda
#   8. Freeze: CorePAM_weights.csv, CorePAM_model.rds,
#              CorePAM_training_card.json, pareto_df_cindex_oof.csv
#
# INPUTS:
#   01_Base_Pura_CorePAM/PROCESSED/SCANB/expression_genelevel_preZ.parquet
#   01_Base_Pura_CorePAM/PROCESSED/SCANB/clinical_FINAL.parquet
#
# OUTPUTS (all in results/corepam/):
#   pareto_df_cindex_oof.csv
#   CorePAM_weights.csv
#   CorePAM_model.rds
#   CorePAM_training_card.json
#   artifact_hashes.csv
# =============================================================================

source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(glmnet)
  library(survival)
  library(jsonlite)
})

SCRIPT_NAME <- "05_reduce_pam50_to_corepam_FINAL.R"

PAM50_GENES <- c(
  "ACTR3B","ANLN","BAG1","BCL2","BIRC5","BLVRA","CCNB1","CCNE1","CDC20",
  "CDC6","CDH3","CENPF","CEP55","CXXC5","EGFR","ERBB2","ESR1","EXO1",
  "FGFR4","FOXA1","FOXC1","GPR160","GRB7","KIF2C","KRT14","KRT17","KRT5",
  "MAPT","MDM2","MELK","MIA","MKI67","MLPH","MMP11","MYBL2","MYC","NAT1",
  "NDC80","NUF2","ORC6","PGR","PHGDH","PTTG1","RRM2","SFRP1","SLC39A6",
  "TMEM45B","TYMS","UBE2C","UBE2T"
)

# =============================================================================
# 1) LOAD SCANB DATA
# =============================================================================
message("[05] Loading SCANB expression and clinical data...")

expr <- strict_parquet(file.path(proc_cohort("SCANB"),
                                 "expression_genelevel_preZ.parquet"))
clin <- strict_parquet(file.path(proc_cohort("SCANB"),
                                 "clinical_FINAL.parquet"))

message(sprintf("[05] Expression: %d genes x %d samples", nrow(expr), ncol(expr) - 1))
message(sprintf("[05] Clinical: %d samples | OS events: %d (%.1f%%)",
                nrow(clin),
                sum(clin$os_event, na.rm = TRUE),
                100 * mean(clin$os_event, na.rm = TRUE)))

# =============================================================================
# 2) ALIGN SAMPLES (expression x clinical)
# =============================================================================
expr_mat  <- as.matrix(expr[, -1])
rownames(expr_mat) <- expr$gene
expr_mat  <- t(expr_mat)           # samples x genes
rownames(expr_mat) <- normalize_id(rownames(expr_mat))

clin_ids  <- normalize_id(clin$sample_id)
common_ids <- intersect(rownames(expr_mat), clin_ids)

if (length(common_ids) == 0) {
  stop("[05] No samples in common between expression and SCANB clinical data.")
}
message(sprintf("[05] Common samples: %d", length(common_ids)))

expr_mat <- expr_mat[common_ids, , drop = FALSE]
clin_sub <- clin[match(common_ids, clin_ids), ]

# =============================================================================
# 3) FILTER PAM50 GENES AND APPLY INTRA-SCANB Z-SCORE
#    Z-score here is only for derivation (coefficients in SD units).
#    Validation Z-score will be intra-cohort in 06_zscore_and_score_<COHORT>.R
# =============================================================================
pam50_available <- intersect(PAM50_GENES, colnames(expr_mat))
pam50_missing   <- setdiff(PAM50_GENES, colnames(expr_mat))

if (length(pam50_missing) > 0) {
  message(sprintf("[05] WARNING: %d PAM50 genes missing in SCANB: %s",
                  length(pam50_missing), paste(pam50_missing, collapse = ",")))
}
if (length(pam50_available) < 40) {
  stop(sprintf("[05] Only %d PAM50 genes available — minimum 40 required for derivation.",
               length(pam50_available)))
}
message(sprintf("[05] PAM50 genes available for derivation: %d/50",
                length(pam50_available)))

X_raw <- expr_mat[, pam50_available, drop = FALSE]

# Z-score per gene intra-SCANB (for derivation)
X <- scale(X_raw)
# Remove genes with no variance (sd=0 → NA column after scale)
zero_var <- apply(X, 2, function(col) all(is.na(col)))
if (any(zero_var)) {
  message(sprintf("[05] Removing %d genes with zero variance.", sum(zero_var)))
  X <- X[, !zero_var, drop = FALSE]
}
message(sprintf("[05] Final X matrix: %d samples x %d genes (Z-scored)",
                nrow(X), ncol(X)))

# Survival response
y <- survival::Surv(time  = clin_sub$os_time_months,
                    event = clin_sub$os_event)

# =============================================================================
# 4) FIXED STRATIFIED FOLDS (K=10, frozen seed)
#    Stratification by event to balance events/censored across folds.
# =============================================================================
set.seed(FREEZE$seed_folds)

# Manual stratification by event status
event_idx   <- which(clin_sub$os_event == 1L)
censor_idx  <- which(clin_sub$os_event == 0L)

foldid <- integer(nrow(X))
foldid[event_idx]  <- sample(rep(1:FREEZE$k_folds,
                                  length.out = length(event_idx)))
foldid[censor_idx] <- sample(rep(1:FREEZE$k_folds,
                                  length.out = length(censor_idx)))

message(sprintf("[05] Folds created: K=%d | seed=%d | events per fold: ~%d",
                FREEZE$k_folds, FREEZE$seed_folds,
                round(sum(clin_sub$os_event) / FREEZE$k_folds)))

# =============================================================================
# 5) CV.GLMNET WITH keep=TRUE (stores OOF linear predictors)
# =============================================================================
message("[05] Running cv.glmnet (alpha=", FREEZE$alpha,
        ", K=", FREEZE$k_folds, ")...")

old_warn <- getOption("warn"); options(warn = 0)
cv_fit <- glmnet::cv.glmnet(
  x       = X,
  y       = y,
  family  = "cox",
  alpha   = FREEZE$alpha,
  foldid  = foldid,
  keep    = TRUE,   # stores fit.preval (OOF linear predictors)
  grouped = TRUE
)
options(warn = old_warn)

message(sprintf("[05] Lambda path: %d points | lambda.min=%.6f | lambda.1se=%.6f",
                length(cv_fit$lambda), cv_fit$lambda.min, cv_fit$lambda.1se))

# =============================================================================
# 6) CALCULATE OOF C-INDEX FOR EACH POINT ON PATH
#    fit.preval: n_samples x n_lambdas matrix with OOF linear predictors
# =============================================================================
message("[05] Calculating OOF C-index per lambda...")

oof_lp  <- cv_fit$fit.preval   # n × n_lambda
lambdas <- cv_fit$lambda
dfs     <- cv_fit$glmnet.fit$df  # número de coeficientes não-zero por lambda

pareto_list <- vector("list", length(lambdas))

for (j in seq_along(lambdas)) {
  lp <- oof_lp[, j]

  # Samples with valid predictions (some may be NA at first lambdas)
  valid_idx <- !is.na(lp)
  if (sum(valid_idx) < 20) {
    pareto_list[[j]] <- tibble(
      lambda = lambdas[j], df = dfs[j],
      craw = NA_real_, cadj = NA_real_, n_valid = sum(valid_idx)
    )
    next
  }

  old_warn <- getOption("warn"); options(warn = 0)
  craw <- tryCatch({
    conc <- survival::concordance(y[valid_idx] ~ lp[valid_idx])
    as.numeric(conc$concordance)
  }, error = function(e) NA_real_)
  options(warn = old_warn)

  cadj <- if (!is.na(craw)) max(craw, 1 - craw) else NA_real_

  pareto_list[[j]] <- tibble(
    lambda  = lambdas[j],
    df      = dfs[j],
    craw    = craw,
    cadj    = cadj,
    n_valid = sum(valid_idx)
  )
}

df_pareto <- bind_rows(pareto_list) |>
  filter(!is.na(cadj)) |>
  # For each df: keep only the largest lambda (most parsimonious)
  group_by(df) |>
  slice_max(lambda, n = 1, with_ties = FALSE) |>
  ungroup() |>
  arrange(df)

message(sprintf("[05] Pareto computed: %d unique df points", nrow(df_pareto)))

# =============================================================================
# 7) SELECTION BY NON-INFERIORITY
# =============================================================================
cmax      <- max(df_pareto$cadj, na.rm = TRUE)
threshold <- cmax - FREEZE$delta_c

message(sprintf("[05] Cmax=%.4f | threshold=%.4f (delta_c=%.3f)",
                cmax, threshold, FREEZE$delta_c))

# Smallest df with Cadj >= threshold; tiebreak: largest lambda
selected <- df_pareto |>
  filter(cadj >= threshold) |>
  arrange(df, desc(lambda)) |>
  slice(1)

if (nrow(selected) == 0) {
  stop("[05] No model satisfies the non-inferiority criterion. ",
       "Check delta_c in analysis_freeze.csv.")
}

message(sprintf(
  "[05] Core-PAM SELECTED: df=%d genes | Cadj=%.4f | lambda=%.6f",
  selected$df, selected$cadj, selected$lambda
))

# =============================================================================
# 8) REFIT ON FULL DATA (glmnet.fit at chosen lambda — Memorial §5.4)
#    Use glmnet.fit from cv_fit (full-data path), not single-lambda refit
# =============================================================================
message("[05] Extracting coefficients from full-data path at chosen lambda...")

old_warn <- getOption("warn"); options(warn = 0)
coefs_raw <- coef(cv_fit$glmnet.fit, s = selected$lambda)
options(warn = old_warn)

# Genes with non-zero coefficient
nonzero_idx <- which(coefs_raw != 0)
corepam_genes  <- rownames(coefs_raw)[nonzero_idx]
corepam_weights <- as.numeric(coefs_raw[nonzero_idx])
n_genes <- length(corepam_genes)

if (n_genes != selected$df) {
  message(sprintf("[05] WARNING: path df=%d but non-zero coefs=%d (normal rounding discrepancy).",
                  selected$df, n_genes))
}

message(sprintf("[05] Core-PAM: %d genes | Weights (first 10):", n_genes))
weights_df <- tibble(gene = corepam_genes, weight = corepam_weights) |>
  arrange(desc(abs(weight)))
print(head(weights_df, 10))

# =============================================================================
# 9) SCORE DIRECTION (Memorial §5.5)
#    If HR < 1 in training cohort → invert sign of weights
# =============================================================================
score_train <- as.numeric(X[, corepam_genes, drop = FALSE] %*% corepam_weights)
score_train_z <- scale(score_train)[, 1]

old_warn <- getOption("warn"); options(warn = 0)
cox_dir <- survival::coxph(y ~ score_train_z)
options(warn = old_warn)

hr_train <- exp(coef(cox_dir))
message(sprintf("[05] Training HR (per 1 SD): %.3f", hr_train))

if (hr_train < 1) {
  message("[05] HR < 1 → inverting sign of weights (score_direction = -1)")
  corepam_weights <- -corepam_weights
  score_direction <- -1L
} else {
  message("[05] HR >= 1 → original direction (score_direction = +1)")
  score_direction <- 1L
}

weights_final <- tibble(gene = corepam_genes, weight = corepam_weights) |>
  arrange(desc(abs(weight)))

# =============================================================================
# 10) SAVE ARTIFACTS (results/corepam/)
# =============================================================================
dir.create(PATHS$results$corepam, showWarnings = FALSE, recursive = TRUE)

artifact_hashes <- list()

# --- 10a) pareto_df_cindex_oof.csv ---
pareto_path <- file.path(PATHS$results$corepam, "pareto_df_cindex_oof.csv")
write_csv(df_pareto, pareto_path)
h <- sha256_file(pareto_path)
registry_append("SCANB", "Pareto_OOF", pareto_path, h, "INTEGRO", SCRIPT_NAME,
                file.info(pareto_path)$size / 1024^2)
artifact_hashes$pareto_df_cindex_oof <- h
message("[05] Saved: ", pareto_path)

# --- 10b) CorePAM_weights.csv ---
weights_path <- file.path(PATHS$results$corepam, "CorePAM_weights.csv")
write_csv(weights_final, weights_path)
h <- sha256_file(weights_path)
registry_append("SCANB", "CorePAM_Weights", weights_path, h, "INTEGRO", SCRIPT_NAME,
                file.info(weights_path)$size / 1024^2)
artifact_hashes$CorePAM_weights <- h
message("[05] Saved: ", weights_path)

# --- 10c) CorePAM_model.rds (complete glmnet.fit) ---
model_path <- file.path(PATHS$results$corepam, "CorePAM_model.rds")
saveRDS(cv_fit$glmnet.fit, model_path)
h <- sha256_file(model_path)
registry_append("SCANB", "CorePAM_Model", model_path, h, "INTEGRO", SCRIPT_NAME,
                file.info(model_path)$size / 1024^2)
artifact_hashes$CorePAM_model <- h
message("[05] Saved: ", model_path)

# --- 10d) CorePAM_training_card.json ---
training_card <- list(
  project          = "Core-PAM",
  memorial_version = "v6.1",
  script           = SCRIPT_NAME,
  timestamp        = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
  cohort_train     = "SCANB",
  n_samples_train  = nrow(X),
  n_events_train   = sum(clin_sub$os_event, na.rm = TRUE),
  freeze = list(
    alpha              = FREEZE$alpha,
    k_folds            = FREEZE$k_folds,
    seed_folds         = FREEZE$seed_folds,
    delta_c            = FREEZE$delta_c,
    min_genes_fraction = FREEZE$min_genes_fraction
  ),
  derivation = list(
    pam50_genes_available = length(pam50_available),
    pam50_genes_missing   = pam50_missing,
    lambda_chosen  = selected$lambda,
    df_chosen      = selected$df,
    cadj_chosen    = selected$cadj,
    cmax           = cmax,
    threshold      = threshold,
    score_direction = score_direction
  ),
  corepam = list(
    n_genes  = n_genes,
    genes    = corepam_genes,
    weights  = corepam_weights
  ),
  artifact_hashes = artifact_hashes
)

card_path <- file.path(PATHS$results$corepam, "CorePAM_training_card.json")
jsonlite::write_json(training_card, card_path, pretty = TRUE, auto_unbox = TRUE)
h <- sha256_file(card_path)
registry_append("SCANB", "CorePAM_Training_Card", card_path, h,
                "INTEGRO", SCRIPT_NAME,
                file.info(card_path)$size / 1024^2)
artifact_hashes$CorePAM_training_card <- h
message("[05] Saved: ", card_path)

# --- 10e) artifact_hashes.csv ---
hashes_df   <- tibble(artifact = names(artifact_hashes),
                      sha256   = unlist(artifact_hashes))
hashes_path <- file.path(PATHS$results$corepam, "artifact_hashes.csv")
write_csv(hashes_df, hashes_path)
message("[05] Saved: ", hashes_path)

# =============================================================================
# 11) FINAL REPORT (printed to console)
# =============================================================================
message("\n", strrep("=", 60))
message("[05] CORE-PAM — FINAL RESULT (FREEZE)")
message(strrep("=", 60))
message(sprintf("  Core-PAM genes    : %d", n_genes))
message(sprintf("  Cadj OOF          : %.4f", selected$cadj))
message(sprintf("  Cmax PAM50        : %.4f", cmax))
message(sprintf("  Margin used       : %.3f (frozen delta_c)", FREEZE$delta_c))
message(sprintf("  Lambda chosen     : %.6f", selected$lambda))
message(sprintf("  Score direction   : %+d", score_direction))
message(sprintf("  Training HR (1SD) : %.3f", exp(coef(cox_dir)) * score_direction^2))
message("\n  Genes and weights:")
print(weights_final, n = Inf)
message(strrep("=", 60))
message("[05] Artifacts in: ", PATHS$results$corepam)
message("[05] Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message("Next step: scripts/06_zscore_and_score_<COHORT>.R")
