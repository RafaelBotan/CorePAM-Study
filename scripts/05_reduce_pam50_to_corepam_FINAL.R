# =============================================================================
# SCRIPT: 05_reduce_pam50_to_corepam_FINAL.R
# PURPOSE: Reduce PAM50 → CorePAM (minimal gene set) via OOF non-inferiority.
# PROJETO: Core-PAM (Memorial v6.1 §5)
#
# DESIGN (A1-proof, audit-friendly):
#   • Train only: SCAN-B (all N=3069; Option A — no internal split)
#   • Elastic-Net Cox (alpha = FREEZE$alpha = 0.5)
#   • K=10 folds — deterministic via SHA-256(patient_id) + event stratification
#     → Reproducible WITHOUT set.seed() as primary driver
#     → Same patient always assigned to same fold, regardless of data order
#   • OOF Harrell C-index computed per unique df level in lambda path
#   • Core-PAM = smallest df with C_adj ≥ C_max − FREEZE$delta_c
#   • Gene count (XX) is DERIVED, not pre-specified
#
# INPUTS:
#   01_Base_Pura_CorePAM/PROCESSED/SCANB/expression_genelevel_preZ.parquet
#   01_Base_Pura_CorePAM/PROCESSED/SCANB/clinical_FINAL.parquet
#
# OUTPUTS:
#   results/corepam/pareto_df_cindex_oof.csv
#   results/corepam/selected_CorePAM_summary.json
#   results/corepam/CorePAM_weights.csv
#   results/corepam/CorePAM_model.rds
#   results/corepam/artifact_hashes.csv
# =============================================================================

source("scripts/00_setup.R")

SCRIPT_NAME <- "05_reduce_pam50_to_corepam_FINAL.R"
COHORT      <- "SCANB"

# Load modeling packages
old_warn_pkg <- getOption("warn"); options(warn = 0)
suppressPackageStartupMessages({
  library(glmnet)
  library(survival)
})
options(warn = old_warn_pkg)

# =============================================================================
# 0) FROZEN PARAMETERS
# =============================================================================
ALPHA_EN <- FREEZE$alpha                # 0.5 (elastic-net mixing)
K_FOLDS  <- as.integer(FREEZE$k_folds) # 10
DELTA_C  <- FREEZE$delta_c             # 0.010 (non-inferiority margin)

message(sprintf("[05] Frozen: alpha=%.1f | K=%d | delta_c=%.3f", ALPHA_EN, K_FOLDS, DELTA_C))
message("[05] Fold assignment: deterministic (SHA-256 of patient_id + event stratification)")
message("[05] NOTE: FREEZE$seed_folds is NOT the primary selection driver in this design.")

# =============================================================================
# 1) PAM50 CANDIDATE GENE LIST
#    50 genes — canonical PAM50 (Parker et al. 2009, J Clin Oncol)
#    Reference: Table 1 of Parker JS et al. Supervised Risk Predictor of Breast
#    Cancer Based on Intrinsic Subtypes. J Clin Oncol 2009;27:1160-1167.
#    KRT8 is NOT included (it is not part of the 50-gene canonical list).
# =============================================================================
PAM50_GENES <- c(
  "ACTR3B","ANLN","BAG1","BCL2","BIRC5","BLVRA","CCNB1","CCNE1","CDC20","CDC6",
  "CDH3","CENPF","CEP55","CXXC5","EGFR","ERBB2","ESR1","EXO1","FGFR4","FOXA1",
  "FOXC1","GPR160","GRB7","KIF2C","KRT14","KRT17","KRT5","MAPT","MDM2",
  "MELK","MIA","MKI67","MLPH","MMP11","MYBL2","MYC","NAT1","NDC80","NUF2",
  "ORC6","PGR","PHGDH","PTTG1","RRM2","SFRP1","SLC39A6","TMEM45B","TYMS","UBE2C","UBE2T"
)
message(sprintf("[05] PAM50 candidate pool: %d genes", length(PAM50_GENES)))

# =============================================================================
# 2) LOAD EXPRESSION (PAM50 genes only — memory-efficient lazy filter)
# =============================================================================
path_expr <- file.path(proc_cohort(COHORT), "expression_genelevel_preZ.parquet")
message("[05] Loading PAM50 genes from expression parquet...")

old_warn_arrow <- getOption("warn"); options(warn = 0)
expr_full <- arrow::read_parquet(path_expr)
pam50_df  <- expr_full[expr_full$gene %in% PAM50_GENES, ]
rm(expr_full)
options(warn = old_warn_arrow)

gene_names_found   <- pam50_df$gene
missing_genes      <- setdiff(PAM50_GENES, gene_names_found)
X_genes_x_samples <- as.matrix(pam50_df[, !names(pam50_df) %in% "gene"])
rownames(X_genes_x_samples) <- gene_names_found

if (length(missing_genes) > 0) {
  message(sprintf("[05] PAM50 genes absent in SCANB: %s", paste(missing_genes, collapse=",")))
}
message(sprintf("[05] PAM50 genes loaded: %d/%d | absent: %d",
                length(gene_names_found), length(PAM50_GENES), length(missing_genes)))

# Transpose: samples × genes (glmnet convention: X is n_samples × n_features)
X_raw <- t(X_genes_x_samples)   # rownames = patient_ids, colnames = genes
message(sprintf("[05] X transposed: %d samples × %d genes", nrow(X_raw), ncol(X_raw)))

# =============================================================================
# 3) LOAD CLINICAL + JOIN ON patient_id
# =============================================================================
path_clin <- file.path(proc_cohort(COHORT), "clinical_FINAL.parquet")
clin      <- strict_parquet(path_clin)

common_ids <- intersect(clin$patient_id, rownames(X_raw))
message(sprintf("[05] Sample overlap: %d / %d clinical | %d expression",
                length(common_ids), nrow(clin), nrow(X_raw)))

clin_sub <- clin[match(common_ids, clin$patient_id), ]
X_sub    <- X_raw[common_ids, , drop = FALSE]

df_train <- data.frame(
  patient_id = common_ids,
  OS_Time    = clin_sub$os_time_months,
  OS_Status  = as.integer(clin_sub$os_event),
  stringsAsFactors = FALSE
)

# =============================================================================
# 4) QC GUARDS
# =============================================================================
n_samples <- nrow(df_train)
n_events  <- sum(df_train$OS_Status)
n_genes   <- ncol(X_sub)

message(sprintf("[05] QC TRAIN | N=%d | events=%d (%.1f%%) | genes=%d",
                n_samples, n_events, 100 * n_events / n_samples, n_genes))

if (n_samples < 500)
  stop("[05] ABORT: too few samples (< 500).")
if (n_events < 100)
  stop("[05] ABORT: too few events (< 100) for penalized Cox.")
if (n_genes < floor(0.90 * length(PAM50_GENES)))
  stop("[05] ABORT: > 10% of PAM50 candidate genes missing in training data.")
if (any(df_train$OS_Time <= 0, na.rm = TRUE))
  stop("[05] ABORT: OS_Time <= 0 detected.")
if (any(!df_train$OS_Status %in% c(0L, 1L)))
  stop("[05] ABORT: OS_Status values outside {0, 1}.")

# =============================================================================
# 5) SCALE X (once upstream; standardize=FALSE in all glmnet calls)
# =============================================================================
old_warn_sc <- getOption("warn"); options(warn = 0)
X_scaled <- scale(X_sub)
options(warn = old_warn_sc)

sd_check <- attr(X_scaled, "scaled:scale")
if (any(!is.finite(sd_check)) || any(sd_check == 0))
  stop("[05] ABORT: zero or non-finite SD in a PAM50 gene — check expression preprocessing.")

y_surv <- survival::Surv(df_train$OS_Time, df_train$OS_Status)
message("[05] X scaled (center + scale). standardize=FALSE used in all glmnet calls.")

# =============================================================================
# 6) DETERMINISTIC FOLDID — SHA-256(patient_id) + event stratification
# =============================================================================
hash_to_int <- function(x) {
  h <- substr(digest::digest(x, algo = "sha256"), 1, 8)
  strtoi(h, base = 16L)
}

make_foldid_stratified <- function(patient_id, event, K) {
  event  <- as.integer(event > 0)
  id_int <- vapply(patient_id, hash_to_int, integer(1L))
  foldid <- integer(length(patient_id))
  for (e in c(0L, 1L)) {
    idx <- which(event == e)
    ord <- idx[order(id_int[idx])]
    foldid[ord] <- rep(seq_len(K), length.out = length(ord))
  }
  foldid
}

foldid <- make_foldid_stratified(df_train$patient_id, df_train$OS_Status, K_FOLDS)

fold_events <- vapply(seq_len(K_FOLDS),
                      function(k) sum(df_train$OS_Status[foldid == k]),
                      integer(1L))
message(sprintf("[05] Foldid assigned: K=%d | events/fold: %s",
                K_FOLDS, paste(fold_events, collapse = "/")))

# =============================================================================
# 7) CV.GLMNET — full lambda path (foldid fixed; used for path definition)
# =============================================================================
message("[05] Running cv.glmnet (defines lambda path)...")
old_warn_cv <- getOption("warn"); options(warn = 0)
cv_full <- glmnet::cv.glmnet(
  x           = X_scaled,
  y           = y_surv,
  family      = "cox",
  alpha       = ALPHA_EN,
  foldid      = foldid,
  standardize = FALSE
)
options(warn = old_warn_cv)

lambda_path <- cv_full$glmnet.fit$lambda
df_path     <- as.integer(cv_full$glmnet.fit$df)
message(sprintf("[05] Lambda path: %d steps | df range: %d – %d",
                length(lambda_path), min(df_path), max(df_path)))

# One lambda per unique df (most parsimonious = largest lambda at that df)
path_tbl <- tibble::tibble(lambda = lambda_path, df = df_path) |>
  dplyr::group_by(df) |>
  dplyr::summarise(lambda = max(lambda), .groups = "drop") |>
  dplyr::arrange(df)

message(sprintf("[05] Unique df levels to evaluate: %d", nrow(path_tbl)))

# =============================================================================
# 8) OOF C-INDEX PER df LEVEL
# =============================================================================
harrell_c <- function(time, status, lp) {
  old_w <- getOption("warn"); options(warn = 0)
  val   <- as.numeric(survival::concordance(survival::Surv(time, status) ~ lp)$concordance)
  options(warn = old_w)
  val
}
c_adj_fn <- function(c_raw) max(c_raw, 1 - c_raw, na.rm = TRUE)

oof_lp <- function(X, time, status, foldid, lam) {
  K  <- max(foldid)
  lp <- rep(NA_real_, nrow(X))
  for (k in seq_len(K)) {
    tr <- which(foldid != k)
    te <- which(foldid == k)
    old_w <- getOption("warn"); options(warn = 0)
    fit_k <- glmnet::glmnet(
      x           = X[tr, , drop = FALSE],
      y           = survival::Surv(time[tr], status[tr]),
      family      = "cox",
      alpha       = ALPHA_EN,
      lambda      = lam,
      standardize = FALSE
    )
    lp[te] <- as.numeric(stats::predict(fit_k,
                                         newx = X[te, , drop = FALSE],
                                         s    = lam,
                                         type = "link"))
    options(warn = old_w)
  }
  lp
}

message(sprintf("[05] OOF C-index loop: %d levels × K=%d folds...",
                nrow(path_tbl), K_FOLDS))
oof_list <- vector("list", nrow(path_tbl))

for (i in seq_len(nrow(path_tbl))) {
  lam <- path_tbl$lambda[i]
  dfi <- path_tbl$df[i]

  lp_i <- oof_lp(X_scaled, df_train$OS_Time, df_train$OS_Status, foldid, lam)

  if (any(is.na(lp_i))) {
    oof_list[[i]] <- tibble::tibble(df = dfi, lambda = lam,
                                     c_raw = NA_real_, c_adj = NA_real_)
  } else {
    c_raw <- harrell_c(df_train$OS_Time, df_train$OS_Status, lp_i)
    oof_list[[i]] <- tibble::tibble(df = dfi, lambda = lam,
                                     c_raw = c_raw, c_adj = c_adj_fn(c_raw))
  }

  if (i %% 5 == 0 || i == nrow(path_tbl)) {
    message(sprintf("[05]   %d/%d | df=%d | C_adj=%s",
                    i, nrow(path_tbl), dfi,
                    if (is.na(oof_list[[i]]$c_adj)) "NA"
                    else sprintf("%.4f", oof_list[[i]]$c_adj)))
  }
}

res_oof <- dplyr::bind_rows(oof_list) |> dplyr::arrange(df)

out_pareto <- file.path(PATHS$results$corepam, "pareto_df_cindex_oof.csv")
readr::write_csv(res_oof, out_pareto)
message("[05] Pareto table saved: ", basename(out_pareto))

# =============================================================================
# 9) SELECT MINIMAL NON-INFERIOR df
# =============================================================================
res_valid <- dplyr::filter(res_oof, !is.na(c_adj))
if (nrow(res_valid) == 0) stop("[05] ABORT: all OOF evaluations returned NA.")

C_max    <- max(res_valid$c_adj)
C_thresh <- C_max - DELTA_C

candidates <- res_valid |>
  dplyr::filter(c_adj >= C_thresh) |>
  dplyr::arrange(df, dplyr::desc(lambda))

if (nrow(candidates) == 0)
  stop(sprintf("[05] ABORT: no model meets non-inferiority (delta_c=%.3f | C_max=%.4f).",
               DELTA_C, C_max))

chosen        <- candidates[1, ]
lambda_chosen <- chosen$lambda
df_chosen     <- as.integer(chosen$df)

message(strrep("=", 60))
message(sprintf("[05] SELECTED Core-PAM | df = %d genes", df_chosen))
message(sprintf("[05] lambda_chosen = %.6g", lambda_chosen))
message(sprintf("[05] C_adj = %.4f | C_max = %.4f | gap = %.4f",
                chosen$c_adj, C_max, C_max - chosen$c_adj))
message(strrep("=", 60))

# =============================================================================
# 10) FREEZE WEIGHTS (from full-path fit at chosen lambda — avoids df drift)
# =============================================================================
old_warn_coef <- getOption("warn"); options(warn = 0)
b <- as.matrix(coef(cv_full$glmnet.fit, s = lambda_chosen))
options(warn = old_warn_coef)

nonzero_genes   <- rownames(b)[as.numeric(b) != 0]
nonzero_weights <- as.numeric(b[nonzero_genes, 1])

if (length(nonzero_genes) != df_chosen) {
  message(sprintf("[05] NOTE: coef() returned %d nonzero genes; adjusting df to match.",
                  length(nonzero_genes)))
  df_chosen <- length(nonzero_genes)
}

message(sprintf("[05] Core-PAM genes (%d): %s",
                length(nonzero_genes), paste(sort(nonzero_genes), collapse = ", ")))

weights_tbl <- tibble::tibble(
  gene   = nonzero_genes,
  weight = nonzero_weights
) |> dplyr::arrange(dplyr::desc(abs(weight)))

out_weights <- file.path(PATHS$results$corepam, "CorePAM_weights.csv")
readr::write_csv(weights_tbl, out_weights)
message("[05] Weights saved: ", basename(out_weights))

# =============================================================================
# 11) MODEL OBJECT (for scoring in scripts 06–07)
# =============================================================================
model_obj <- list(
  glmnet_fit     = cv_full$glmnet.fit,
  lambda_chosen  = lambda_chosen,
  df_chosen      = length(nonzero_genes),
  alpha          = ALPHA_EN,
  genes_selected = nonzero_genes,
  weights        = weights_tbl,
  pam50_candidate = PAM50_GENES,
  pam50_present  = gene_names_found,
  foldid         = foldid,
  scale_center   = attr(X_scaled, "scaled:center"),
  scale_sd       = attr(X_scaled, "scaled:scale"),
  c_adj_chosen   = chosen$c_adj,
  c_max_oof      = C_max,
  delta_c        = DELTA_C,
  n_train        = n_samples,
  n_events_train = n_events,
  script         = SCRIPT_NAME,
  timestamp      = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)

out_model <- file.path(PATHS$results$corepam, "CorePAM_model.rds")
old_warn_rds <- getOption("warn"); options(warn = 0)
saveRDS(model_obj, out_model)
options(warn = old_warn_rds)
message("[05] Model saved: ", basename(out_model))

# =============================================================================
# 12) TRAINING CARD (JSON)
# =============================================================================
summary_json <- list(
  method         = "Core-PAM: PAM50 reduction via OOF non-inferiority (Harrell C-index)",
  train_cohort   = COHORT,
  alpha_en       = ALPHA_EN,
  k_folds        = K_FOLDS,
  delta_c        = DELTA_C,
  fold_method    = "Deterministic SHA-256(patient_id) + event stratification (no seed as driver)",
  n_train        = n_samples,
  n_events_train = n_events,
  pam50_candidate_n  = length(PAM50_GENES),
  pam50_present_n    = length(gene_names_found),
  pam50_missing_list = missing_genes,
  selected = list(
    df     = length(nonzero_genes),
    lambda = lambda_chosen,
    c_adj  = chosen$c_adj,
    c_max  = C_max,
    gap    = C_max - chosen$c_adj,
    genes  = sort(nonzero_genes)
  ),
  input_sha256 = list(
    expr = sha256_file(path_expr),
    clin = sha256_file(path_clin)
  ),
  timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)

out_json <- file.path(PATHS$results$corepam, "selected_CorePAM_summary.json")
jsonlite::write_json(summary_json, out_json, pretty = TRUE, auto_unbox = TRUE)
message("[05] Training card saved: ", basename(out_json))

# =============================================================================
# 13) ARTIFACT HASHES + REGISTRY
# =============================================================================
artifacts  <- c(out_pareto, out_weights, out_model, out_json)
hashes_tbl <- tibble::tibble(
  artifact = basename(artifacts),
  path     = artifacts,
  sha256   = vapply(artifacts, sha256_file, character(1L))
)

out_hashes <- file.path(PATHS$results$corepam, "artifact_hashes.csv")
readr::write_csv(hashes_tbl, out_hashes)

for (i in seq_along(artifacts)) {
  registry_append(
    cohort    = COHORT,
    file_type = sub("\\..*$", "", basename(artifacts[i])),
    file_path = artifacts[i],
    sha256    = hashes_tbl$sha256[i],
    status    = "INTEGRO",
    script    = SCRIPT_NAME,
    size_mb   = file.info(artifacts[i])$size / 1024^2
  )
}
message("[05] Artifact hashes and registry updated.")

# =============================================================================
# FINAL SUMMARY
# =============================================================================
message("\n", strrep("=", 60))
message("[05] Core-PAM derivation complete")
message(strrep("=", 60))
message(sprintf("  Selected: %d genes", length(nonzero_genes)))
message(sprintf("  OOF C_adj: %.4f", chosen$c_adj))
message(sprintf("  C_max:     %.4f", C_max))
message(sprintf("  Gap:       %.4f (threshold: delta_c=%.3f)", C_max - chosen$c_adj, DELTA_C))
message(sprintf("  Genes: %s", paste(sort(nonzero_genes), collapse = ", ")))
message(strrep("=", 60))
message("\nNext step: scripts/06_zscore_and_score_SCANB.R")
message("[05] Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
