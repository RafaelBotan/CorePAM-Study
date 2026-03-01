# =============================================================================
# SCRIPT: 21_pCR_logistic_analysis.R
# PURPOSE: CorePAM score association with pCR (pathologic complete response
#          to NACT) across all pCR-block cohorts. Implements Memorial §11 protocol.
#
#   Step 2: Univariate logistic + AUC (primary)
#   Step 3: Incremental value vs CORE-NACT baseline (ER + HER2)
#   QA:     audit_pcr_definition.csv, pcr_validation_table.csv
#
# OUTPUTS (all in results/pcr/):
#   audit_pcr_definition.csv      — pCR endpoint mapping audit (reviewer-critical)
#   pcr_validation_table.csv      — minimum fields per cohort (Table T_pCR_2)
#   pCR_results_by_cohort.csv     — full per-cohort results
#   pcr_incremental_table.csv     — ΔAUC baseline vs baseline+score (T_pCR_3)
#   T_pCR_1_Cohort_Characteristics.csv — cohort characteristics (T_pCR_1)
#
# PROJETO: Core-PAM (Memorial v6.1 §11 / Freeze Core-PAM — pCR block)
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "21_pCR_logistic_analysis.R"

suppressPackageStartupMessages({
  library(pROC)
})

set.seed(FREEZE$seed_folds)
FORCE <- as.logical(Sys.getenv("FORCE_RERUN", "FALSE"))

out_path      <- file.path(PATHS$results$pcr, "pCR_results_by_cohort.csv")
out_val_table <- file.path(PATHS$results$pcr, "pcr_validation_table.csv")
out_incr      <- file.path(PATHS$results$pcr, "pcr_incremental_table.csv")
out_audit     <- file.path(PATHS$results$pcr, "audit_pcr_definition.csv")
out_char      <- file.path(PATHS$results$pcr, "T_pCR_1_Cohort_Characteristics.csv")

if (!FORCE && all(file.exists(c(out_path, out_val_table, out_incr, out_audit, out_char)))) {
  message(sprintf("[%s] All outputs exist — skipping. Set FORCE_RERUN=TRUE to rerun.",
                  SCRIPT_NAME))
  quit(save = "no", status = 0)
}

# --------------------------------------------------------------------------
# Load cohort manifest
# --------------------------------------------------------------------------
manifest_path <- file.path(PATHS$config, "pcr_cohort_manifest.csv")
if (file.exists(manifest_path)) {
  manifest <- readr::read_csv(manifest_path, show_col_types = FALSE)
  PCR_COHORTS <- manifest$cohort[manifest$status != "excluded"]
} else {
  PCR_COHORTS <- c("GSE25066", "GSE20194", "GSE32646", "ISPY1")
  warning(sprintf("[%s] Manifest not found at %s; using default cohort list.",
                  SCRIPT_NAME, manifest_path))
}
message(sprintf("[%s] Cohorts: %s", SCRIPT_NAME, paste(PCR_COHORTS, collapse = ", ")))

# --------------------------------------------------------------------------
# Helper: bootstrap OR CI
# --------------------------------------------------------------------------
bootstrap_or_ci <- function(pcr, score_z, n_boot = 1000, seed = 42) {
  set.seed(seed)
  n <- length(pcr); or_boot <- numeric(n_boot)
  old_warn <- getOption("warn"); options(warn = 0)
  for (i in seq_len(n_boot)) {
    idx   <- sample(n, n, replace = TRUE)
    fit_b <- tryCatch(
      glm(pcr[idx] ~ score_z[idx], family = binomial),
      error   = function(e) NULL, warning = function(w) glm(pcr[idx] ~ score_z[idx], family = binomial)
    )
    or_boot[i] <- if (!is.null(fit_b)) exp(coef(fit_b)[2]) else NA_real_
  }
  options(warn = old_warn)
  list(ci_low = quantile(or_boot, 0.025, na.rm = TRUE),
       ci_high = quantile(or_boot, 0.975, na.rm = TRUE))
}

# --------------------------------------------------------------------------
# Helper: logistic calibration (Hosmer-Lemeshow style + intercept/slope)
# --------------------------------------------------------------------------
logistic_calibration <- function(obs, prob) {
  old_warn <- getOption("warn"); options(warn = 0)
  cal_fit <- tryCatch(
    glm(obs ~ log(prob / (1 - prob)), family = binomial),
    error = function(e) NULL, warning = function(w) glm(obs ~ log(prob / (1 - prob)), family = binomial)
  )
  options(warn = old_warn)
  if (is.null(cal_fit)) return(list(intercept = NA, slope = NA))
  list(intercept = round(coef(cal_fit)[1], 4),
       slope     = round(coef(cal_fit)[2], 4))
}

# --------------------------------------------------------------------------
# Helper: Brier score
# --------------------------------------------------------------------------
brier_score <- function(obs, prob) round(mean((obs - prob)^2, na.rm = TRUE), 5)

# --------------------------------------------------------------------------
# Main analysis loop
# --------------------------------------------------------------------------
results_list    <- list()
audit_list      <- list()
incr_list       <- list()
char_list       <- list()
val_table_list  <- list()
er_strat_list_all <- list()

for (cohort in PCR_COHORTS) {
  message(sprintf("\n[%s] === %s ===", SCRIPT_NAME, cohort))

  ready_path <- file.path(proc_pcr_cohort(cohort), "analysis_ready.parquet")
  if (!file.exists(ready_path)) {
    message(sprintf("[%s] [%s] analysis_ready.parquet not found — skipping.",
                    SCRIPT_NAME, cohort))
    next
  }

  df <- strict_parquet(ready_path)
  df <- df[!is.na(df$pcr) & !is.na(df$score_z), ]

  n_total  <- nrow(df)
  n_pcr1   <- sum(df$pcr == 1L)
  n_rd     <- sum(df$pcr == 0L)
  pcr_rate <- n_pcr1 / n_total

  # Gene coverage (from manifest or analysis_ready)
  n_genes_present <- if ("genes_present" %in% names(df)) df$genes_present[1] else NA_integer_

  message(sprintf("[%s] [%s] N=%d | pCR=%d | RD=%d (%.1f%%)",
                  SCRIPT_NAME, cohort, n_total, n_pcr1, n_rd, 100 * pcr_rate))

  # --- Cohort characteristics row (T_pCR_1) ---
  char_list[[cohort]] <- tibble(
    cohort        = cohort,
    endpoint      = "pCR",
    n_total       = n_total,
    n_pcr1        = n_pcr1,
    n_rd          = n_rd,
    pcr_rate_pct  = round(100 * pcr_rate, 1),
    genes_present = n_genes_present,
    has_er        = "er_status" %in% names(df) && mean(!is.na(df$er_status)) >= 0.5,
    has_her2      = "her2_status" %in% names(df) && mean(!is.na(df$her2_status)) >= 0.5
  )

  if (n_total < 30 || n_pcr1 < 5 || n_rd < 5) {
    message(sprintf("[%s] [%s] Insufficient events — skipping inference.", SCRIPT_NAME, cohort))
    next
  }

  # ------- Univariate logistic -----------------------------------------------
  old_warn <- getOption("warn"); options(warn = 0)
  fit_uni <- tryCatch(
    glm(pcr ~ score_z, data = df, family = binomial),
    error   = function(e) { message("[", SCRIPT_NAME, "] [", cohort, "] glm error: ",
                                    conditionMessage(e)); NULL },
    warning = function(w) glm(pcr ~ score_z, data = df, family = binomial)
  )
  options(warn = old_warn)
  if (is.null(fit_uni)) next

  sm_uni     <- summary(fit_uni)
  or_uni     <- exp(coef(fit_uni)[["score_z"]])
  se_logOR   <- sm_uni$coefficients["score_z", "Std. Error"]
  p_uni      <- sm_uni$coefficients["score_z", "Pr(>|z|)"]
  or_lo_wald <- exp(coef(fit_uni)[["score_z"]] - 1.96 * se_logOR)
  or_hi_wald <- exp(coef(fit_uni)[["score_z"]] + 1.96 * se_logOR)

  message(sprintf("[%s] [%s] OR=%.3f (%.3f-%.3f) p=%.4g",
                  SCRIPT_NAME, cohort, or_uni, or_lo_wald, or_hi_wald, p_uni))

  # Bootstrap CI
  message(sprintf("[%s] [%s] Bootstrap OR CI (B=%d)...", SCRIPT_NAME, cohort, FREEZE$bootstrap_n))
  boot_ci <- bootstrap_or_ci(df$pcr, df$score_z,
                              n_boot = as.integer(FREEZE$bootstrap_n),
                              seed   = as.integer(FREEZE$seed_folds))

  # ------- AUC (DeLong) -------------------------------------------------------
  old_warn <- getOption("warn"); options(warn = 0)
  roc_obj <- tryCatch(pROC::roc(df$pcr, df$score_z, direction = "<", quiet = TRUE),
                      error = function(e) NULL)
  options(warn = old_warn)

  if (is.null(roc_obj)) {
    auc_val <- NA_real_; auc_lo <- NA_real_; auc_hi <- NA_real_
  } else {
    auc_val <- as.numeric(pROC::auc(roc_obj))
    old_warn <- getOption("warn"); options(warn = 0)
    ci_obj  <- tryCatch(pROC::ci.auc(roc_obj, method = "delong"), error = function(e) NULL)
    options(warn = old_warn)
    auc_lo  <- if (!is.null(ci_obj)) ci_obj[1] else NA_real_
    auc_hi  <- if (!is.null(ci_obj)) ci_obj[3] else NA_real_
    message(sprintf("[%s] [%s] AUC=%.4f (%.4f-%.4f DeLong)", SCRIPT_NAME, cohort,
                    auc_val, auc_lo, auc_hi))
  }

  # ------- Calibration --------------------------------------------------------
  old_warn <- getOption("warn"); options(warn = 0)
  probs    <- tryCatch(predict(fit_uni, type = "response"), error = function(e) NULL)
  options(warn = old_warn)
  cal_res  <- if (!is.null(probs)) logistic_calibration(df$pcr, probs) else list(intercept=NA, slope=NA)
  brier    <- if (!is.null(probs)) brier_score(df$pcr, probs) else NA_real_

  # ------- Adjusted logistic (age + ER if available) --------------------------
  adj_vars  <- c("age", "er_status")
  adj_avail <- adj_vars[sapply(adj_vars, function(v) {
    v %in% names(df) && mean(!is.na(df[[v]])) >= 0.8
  })]

  or_adj <- NA_real_; or_adj_lo <- NA_real_; or_adj_hi <- NA_real_
  p_adj  <- NA_real_; adj_used  <- "none"

  if (length(adj_avail) > 0) {
    df_adj  <- df[complete.cases(df[, c("pcr", "score_z", adj_avail)]), ]
    fml_adj <- as.formula(paste("pcr ~ score_z +", paste(adj_avail, collapse = " + ")))
    old_warn <- getOption("warn"); options(warn = 0)
    fit_adj <- tryCatch(glm(fml_adj, data = df_adj, family = binomial),
                        error = function(e) NULL,
                        warning = function(w) glm(fml_adj, data = df_adj, family = binomial))
    options(warn = old_warn)
    if (!is.null(fit_adj)) {
      sm_adj    <- summary(fit_adj)
      or_adj    <- exp(coef(fit_adj)[["score_z"]])
      se_adj    <- sm_adj$coefficients["score_z", "Std. Error"]
      or_adj_lo <- exp(coef(fit_adj)[["score_z"]] - 1.96 * se_adj)
      or_adj_hi <- exp(coef(fit_adj)[["score_z"]] + 1.96 * se_adj)
      p_adj     <- sm_adj$coefficients["score_z", "Pr(>|z|)"]
      adj_used  <- paste(adj_avail, collapse = "+")
      message(sprintf("[%s] [%s] Adj OR=%.3f (%.3f-%.3f) p=%.4g [%s]",
                      SCRIPT_NAME, cohort, or_adj, or_adj_lo, or_adj_hi, p_adj, adj_used))
    }
  }

  # ------- CORE-NACT incremental value: ΔAUC (ER+HER2 baseline vs +score) ----
  nact_base_vars <- c("er_status", "her2_status")
  nact_avail     <- nact_base_vars[sapply(nact_base_vars, function(v) {
    v %in% names(df) && mean(!is.na(df[[v]])) >= 0.8
  })]

  delta_auc    <- NA_real_; delta_auc_lo <- NA_real_; delta_auc_hi <- NA_real_
  auc_baseline <- NA_real_; auc_full     <- NA_real_
  lrt_p        <- NA_real_; nact_used    <- "none"

  if (length(nact_avail) > 0) {
    df_nact <- df[complete.cases(df[, c("pcr", "score_z", nact_avail)]), ]
    nact_used <- paste(nact_avail, collapse = "+")

    fml_base <- as.formula(paste("pcr ~", paste(nact_avail, collapse = " + ")))
    fml_full <- as.formula(paste("pcr ~", paste(c(nact_avail, "score_z"), collapse = " + ")))

    old_warn <- getOption("warn"); options(warn = 0)
    fit_base <- tryCatch(glm(fml_base, data = df_nact, family = binomial),
                         error = function(e) NULL, warning = function(w) glm(fml_base, data = df_nact, family = binomial))
    fit_full <- tryCatch(glm(fml_full, data = df_nact, family = binomial),
                         error = function(e) NULL, warning = function(w) glm(fml_full, data = df_nact, family = binomial))
    options(warn = old_warn)

    if (!is.null(fit_base) && !is.null(fit_full)) {
      # AUC baseline and full
      p_base <- predict(fit_base, type = "response")
      p_full <- predict(fit_full, type = "response")
      old_warn <- getOption("warn"); options(warn = 0)
      roc_base <- tryCatch(pROC::roc(df_nact$pcr, p_base, quiet = TRUE), error = function(e) NULL)
      roc_full <- tryCatch(pROC::roc(df_nact$pcr, p_full, quiet = TRUE), error = function(e) NULL)
      options(warn = old_warn)
      auc_baseline <- if (!is.null(roc_base)) as.numeric(pROC::auc(roc_base)) else NA_real_
      auc_full     <- if (!is.null(roc_full)) as.numeric(pROC::auc(roc_full)) else NA_real_
      delta_auc    <- auc_full - auc_baseline

      # Bootstrap ΔAUC
      set.seed(as.integer(FREEZE$seed_folds))
      delta_boot <- numeric(as.integer(FREEZE$bootstrap_n))
      old_warn <- getOption("warn"); options(warn = 0)
      n_nact <- nrow(df_nact)
      for (bi in seq_len(as.integer(FREEZE$bootstrap_n))) {
        idx2 <- sample(n_nact, n_nact, replace = TRUE)
        d2   <- df_nact[idx2, ]
        rb   <- tryCatch(as.numeric(pROC::auc(pROC::roc(d2$pcr, predict(
                  glm(fml_base, data = d2, family = binomial), type = "response"), quiet = TRUE))),
                  error = function(e) NA_real_)
        rf   <- tryCatch(as.numeric(pROC::auc(pROC::roc(d2$pcr, predict(
                  glm(fml_full, data = d2, family = binomial), type = "response"), quiet = TRUE))),
                  error = function(e) NA_real_)
        delta_boot[bi] <- rf - rb
      }
      options(warn = old_warn)
      delta_auc_lo <- quantile(delta_boot, 0.025, na.rm = TRUE)
      delta_auc_hi <- quantile(delta_boot, 0.975, na.rm = TRUE)

      # LRT
      old_warn <- getOption("warn"); options(warn = 0)
      lrt_test <- tryCatch(anova(fit_base, fit_full, test = "LRT"), error = function(e) NULL)
      options(warn = old_warn)
      lrt_p <- if (!is.null(lrt_test)) lrt_test[2, "Pr(>Chi)"] else NA_real_

      message(sprintf("[%s] [%s] ΔAUC=%.4f (%.4f-%.4f) LRT p=%.4g [%s]",
                      SCRIPT_NAME, cohort, delta_auc, delta_auc_lo, delta_auc_hi,
                      lrt_p, nact_used))
    }
  }

  # ------- Collect results ---------------------------------------------------
  res_row <- tibble(
    cohort          = cohort,
    n_total         = n_total,
    n_pcr1          = n_pcr1,
    pcr_rate        = round(pcr_rate, 4),
    or_uni          = round(or_uni, 4),
    or_uni_lo_wald  = round(or_lo_wald, 4),
    or_uni_hi_wald  = round(or_hi_wald, 4),
    or_uni_lo_boot  = round(boot_ci$ci_low, 4),
    or_uni_hi_boot  = round(boot_ci$ci_high, 4),
    p_uni           = signif(p_uni, 4),
    logOR_uni       = round(log(or_uni), 6),
    se_logOR_uni    = round(se_logOR, 6),
    auc             = round(auc_val, 4),
    auc_lo95        = round(auc_lo, 4),
    auc_hi95        = round(auc_hi, 4),
    cal_intercept   = cal_res$intercept,
    cal_slope       = cal_res$slope,
    brier_score     = brier,
    or_adj          = round(or_adj, 4),
    or_adj_lo95     = round(or_adj_lo, 4),
    or_adj_hi95     = round(or_adj_hi, 4),
    p_adj           = signif(p_adj, 4),
    adj_vars_used   = adj_used
  )
  results_list[[cohort]] <- res_row

  # ------- ER-stratified analysis (if ER available) --------------------------
  er_strat_rows <- list()
  if ("er_status" %in% names(df) && mean(!is.na(df$er_status)) >= 0.5) {
    for (er_val in c("positive", "negative")) {
      df_er <- df[!is.na(df$er_status) & df$er_status == er_val, ]
      if (nrow(df_er) < 20 || sum(df_er$pcr == 1L) < 5 || sum(df_er$pcr == 0L) < 5) {
        message(sprintf("[%s] [%s] ER=%s: N=%d too small for stratified analysis — skipping",
                        SCRIPT_NAME, cohort, er_val, nrow(df_er)))
        next
      }
      old_warn <- getOption("warn"); options(warn = 0)
      fit_er <- tryCatch(
        glm(pcr ~ score_z, data = df_er, family = binomial),
        error = function(e) NULL,
        warning = function(w) glm(pcr ~ score_z, data = df_er, family = binomial)
      )
      options(warn = old_warn)
      if (!is.null(fit_er)) {
        sm_er     <- summary(fit_er)
        or_er     <- exp(coef(fit_er)[["score_z"]])
        se_er     <- sm_er$coefficients["score_z", "Std. Error"]
        or_er_lo  <- exp(coef(fit_er)[["score_z"]] - 1.96 * se_er)
        or_er_hi  <- exp(coef(fit_er)[["score_z"]] + 1.96 * se_er)
        p_er      <- sm_er$coefficients["score_z", "Pr(>|z|)"]
        message(sprintf("[%s] [%s] ER=%s: N=%d pCR=%d OR=%.3f (%.3f-%.3f) p=%.4g",
                        SCRIPT_NAME, cohort, er_val, nrow(df_er),
                        sum(df_er$pcr == 1L), or_er, or_er_lo, or_er_hi, p_er))
        er_strat_rows[[er_val]] <- tibble(
          cohort    = cohort,
          er_strat  = er_val,
          n_total   = nrow(df_er),
          n_pcr1    = sum(df_er$pcr == 1L),
          or_per_1SD = round(or_er, 4),
          or_lo95   = round(or_er_lo, 4),
          or_hi95   = round(or_er_hi, 4),
          p_wald    = signif(p_er, 4)
        )
      }
    }
  }
  if (length(er_strat_rows) > 0) {
    er_strat_list_all[[cohort]] <- bind_rows(er_strat_rows)
  }

  # pcr_validation_table row (minimum fields per Memorial §11.12)
  genes_total <- if (file.exists(file.path(PATHS$results$corepam, "CorePAM_weights.csv"))) {
    nrow(readr::read_csv(file.path(PATHS$results$corepam, "CorePAM_weights.csv"),
                         show_col_types = FALSE))
  } else NA_integer_
  val_table_list[[cohort]] <- tibble(
    cohort              = cohort,
    endpoint            = "pCR",
    n_total             = n_total,
    n_pcr1              = n_pcr1,
    n_rd                = n_rd,
    genes_total         = genes_total,
    genes_present       = n_genes_present,
    genes_missing       = if (!is.na(n_genes_present) && !is.na(genes_total)) genes_total - n_genes_present else NA,
    coverage_pct        = if (!is.na(n_genes_present) && !is.na(genes_total)) round(100 * n_genes_present / genes_total, 1) else NA,
    or_per_1SD          = round(or_uni, 4),
    or_ci_lo95          = round(or_lo_wald, 4),
    or_ci_hi95          = round(or_hi_wald, 4),
    ci_method           = "Wald (95%)",
    p_value             = signif(p_uni, 4),
    auc                 = round(auc_val, 4),
    auc_lo95            = round(auc_lo, 4),
    auc_hi95            = round(auc_hi, 4),
    baseline_covariates = adj_used
  )

  # incremental value row
  incr_list[[cohort]] <- tibble(
    cohort         = cohort,
    nact_baseline  = nact_used,
    n_complete     = if (length(nact_avail) > 0) nrow(df[complete.cases(df[, c("pcr","score_z",nact_avail)]),]) else NA_integer_,
    auc_baseline   = round(auc_baseline, 4),
    auc_full       = round(auc_full, 4),
    delta_auc      = round(delta_auc, 4),
    delta_auc_lo95 = round(delta_auc_lo, 4),
    delta_auc_hi95 = round(delta_auc_hi, 4),
    lrt_p          = signif(lrt_p, 4)
  )

  # pCR definition audit (row — canonical values from literature/manifest)
  manifest_row <- if (file.exists(manifest_path) && exists("manifest")) {
    manifest[manifest$cohort == cohort, ]
  } else NULL
  audit_list[[cohort]] <- tibble(
    cohort               = cohort,
    geo_accession        = if (!is.null(manifest_row) && nrow(manifest_row) > 0) manifest_row$geo_accession else NA_character_,
    pcr_definition_text  = if (!is.null(manifest_row) && nrow(manifest_row) > 0) manifest_row$pcr_definition_text else NA_character_,
    insitu_included      = if (!is.null(manifest_row) && nrow(manifest_row) > 0) manifest_row$insitu_included else NA_character_,
    n_pcr1               = n_pcr1,
    n_total              = n_total,
    pcr_rate_pct         = round(100 * pcr_rate, 1),
    mapping_note         = "auto_detect: verified against GEO pData on first run"
  )
}

# --------------------------------------------------------------------------
# Save all outputs
# --------------------------------------------------------------------------
if (length(results_list) == 0) stop(sprintf("[%s] No cohort results generated.", SCRIPT_NAME))

dir.create(PATHS$results$pcr, showWarnings = FALSE, recursive = TRUE)

# ER-stratified results
if (length(er_strat_list_all) > 0) {
  er_strat_df <- bind_rows(er_strat_list_all)
  out_er_strat <- file.path(PATHS$results$pcr, "pcr_er_stratified.csv")
  readr::write_csv(er_strat_df, out_er_strat)
  h_er <- sha256_file(out_er_strat)
  registry_append("META_PCR", "pcr_er_stratified", out_er_strat, h_er, "ok", SCRIPT_NAME,
                  file.info(out_er_strat)$size / 1e6)
  message(sprintf("[%s] ER-stratified table saved: %s rows", SCRIPT_NAME, nrow(er_strat_df)))
}

# Full results
res_df <- bind_rows(results_list)
readr::write_csv(res_df, out_path)
h1 <- sha256_file(out_path)
registry_append("META_PCR", "pcr_results_by_cohort", out_path, h1, "ok", SCRIPT_NAME,
                file.info(out_path)$size / 1e6)

# pcr_validation_table (minimum fields — T_pCR_2)
val_df <- bind_rows(val_table_list)
readr::write_csv(val_df, out_val_table)
h2 <- sha256_file(out_val_table); sz2 <- file.info(out_val_table)$size / 1e6
registry_append("META_PCR", "pcr_validation_table", out_val_table, h2, "ok", SCRIPT_NAME, sz2)

# Incremental value (T_pCR_3)
incr_df <- bind_rows(incr_list)
readr::write_csv(incr_df, out_incr)
h3 <- sha256_file(out_incr); sz3 <- file.info(out_incr)$size / 1e6
registry_append("META_PCR", "pcr_incremental_table", out_incr, h3, "ok", SCRIPT_NAME, sz3)

# pCR definition audit
audit_df <- bind_rows(audit_list)
readr::write_csv(audit_df, out_audit)
h4 <- sha256_file(out_audit); sz4 <- file.info(out_audit)$size / 1e6
registry_append("META_PCR", "audit_pcr_definition", out_audit, h4, "ok", SCRIPT_NAME, sz4)

# Cohort characteristics (T_pCR_1)
char_df <- bind_rows(char_list)
readr::write_csv(char_df, out_char)
h5 <- sha256_file(out_char); sz5 <- file.info(out_char)$size / 1e6
registry_append("META_PCR", "T_pCR_1_cohort_char", out_char, h5, "ok", SCRIPT_NAME, sz5)

message(sprintf("\n[%s] === RESULTS SUMMARY ===", SCRIPT_NAME))
print(res_df[, c("cohort", "n_total", "n_pcr1", "or_uni", "or_uni_lo_wald",
                  "or_uni_hi_wald", "p_uni", "auc")])
message(sprintf("[%s] COMPLETED — saved 5 output files to %s",
                SCRIPT_NAME, PATHS$results$pcr))
