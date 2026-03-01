# =============================================================================
# SCRIPT: 21_pCR_logistic_analysis.R
# PURPOSE: CorePAM score association with pCR (pathologic complete response
#          to NACT) across all pCR-block cohorts:
#            - Logistic regression: glm(pcr ~ score_z, family=binomial)
#            - OR (95% CI Wald + bootstrap B=1000)
#            - AUC (pROC DeLong 95% CI)
#            - Logistic with age/ER adjustment where available
#          Output: results/pcr/pCR_results_by_cohort.csv
#
# PROJETO: Core-PAM (Memorial v6.1 / Freeze Core-PAM — pCR block)
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "21_pCR_logistic_analysis.R"

suppressPackageStartupMessages({
  library(pROC)
})

set.seed(FREEZE$seed_folds)
FORCE <- as.logical(Sys.getenv("FORCE_RERUN", "FALSE"))

out_path <- file.path(PATHS$results$pcr, "pCR_results_by_cohort.csv")
if (!FORCE && file.exists(out_path)) {
  message(sprintf("[%s] Output exists — skipping. Set FORCE_RERUN=TRUE to rerun.",
                  SCRIPT_NAME))
  quit(save = "no", status = 0)
}

# --------------------------------------------------------------------------
# pCR cohort list
# --------------------------------------------------------------------------
PCR_COHORTS <- c("GSE25066", "GSE20194", "GSE32646", "ISPY1")

# --------------------------------------------------------------------------
# Helper: bootstrap OR CI
# --------------------------------------------------------------------------
bootstrap_or_ci <- function(pcr, score_z, n_boot = 1000, seed = 42) {
  set.seed(seed)
  n <- length(pcr)
  or_boot <- numeric(n_boot)
  old_warn <- getOption("warn"); options(warn = 0)
  for (i in seq_len(n_boot)) {
    idx <- sample(n, n, replace = TRUE)
    fit_b <- tryCatch(
      glm(pcr[idx] ~ score_z[idx], family = binomial),
      error   = function(e) NULL,
      warning = function(w) NULL
    )
    or_boot[i] <- if (!is.null(fit_b)) exp(coef(fit_b)[2]) else NA_real_
  }
  options(warn = old_warn)
  list(
    ci_low  = quantile(or_boot, 0.025, na.rm = TRUE),
    ci_high = quantile(or_boot, 0.975, na.rm = TRUE)
  )
}

# --------------------------------------------------------------------------
# Main analysis loop
# --------------------------------------------------------------------------
results_list <- list()

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
  pcr_rate <- n_pcr1 / n_total

  message(sprintf("[%s] [%s] N=%d | pCR=%d (%.1f%%)",
                  SCRIPT_NAME, cohort, n_total, n_pcr1, 100 * pcr_rate))

  if (n_total < 30 || n_pcr1 < 5 || (n_total - n_pcr1) < 5) {
    message(sprintf("[%s] [%s] Insufficient events for logistic regression — skipping.",
                    SCRIPT_NAME, cohort))
    next
  }

  # ------- Univariate logistic regression -----------------------------------
  old_warn <- getOption("warn"); options(warn = 0)
  fit_uni <- tryCatch(
    glm(pcr ~ score_z, data = df, family = binomial),
    error   = function(e) { message(sprintf("[%s] [%s] glm error: %s", SCRIPT_NAME, cohort, conditionMessage(e))); NULL },
    warning = function(w) glm(pcr ~ score_z, data = df, family = binomial)
  )
  options(warn = old_warn)

  if (is.null(fit_uni)) next

  sm_uni   <- summary(fit_uni)
  or_uni   <- exp(coef(fit_uni)[["score_z"]])
  se_logOR <- sm_uni$coefficients["score_z", "Std. Error"]
  z_stat   <- sm_uni$coefficients["score_z", "z value"]
  p_uni    <- sm_uni$coefficients["score_z", "Pr(>|z|)"]
  or_lo_wald <- exp(coef(fit_uni)[["score_z"]] - 1.96 * se_logOR)
  or_hi_wald <- exp(coef(fit_uni)[["score_z"]] + 1.96 * se_logOR)

  message(sprintf("[%s] [%s] OR=%.3f (Wald: %.3f-%.3f) p=%.4g",
                  SCRIPT_NAME, cohort, or_uni, or_lo_wald, or_hi_wald, p_uni))

  # Bootstrap CI
  message(sprintf("[%s] [%s] Bootstrap OR CI (B=%d)...", SCRIPT_NAME, cohort,
                  FREEZE$bootstrap_n))
  boot_ci <- bootstrap_or_ci(df$pcr, df$score_z,
                              n_boot = as.integer(FREEZE$bootstrap_n),
                              seed   = as.integer(FREEZE$seed_folds))

  # ------- AUC (DeLong) -------------------------------------------------------
  old_warn <- getOption("warn"); options(warn = 0)
  roc_obj <- tryCatch(
    pROC::roc(df$pcr, df$score_z, direction = "<", quiet = TRUE),
    error = function(e) NULL
  )
  options(warn = old_warn)

  if (is.null(roc_obj)) {
    auc_val <- NA_real_; auc_lo <- NA_real_; auc_hi <- NA_real_
  } else {
    auc_val <- as.numeric(pROC::auc(roc_obj))
    old_warn <- getOption("warn"); options(warn = 0)
    ci_obj  <- tryCatch(
      pROC::ci.auc(roc_obj, method = "delong"),
      error = function(e) NULL
    )
    options(warn = old_warn)
    if (is.null(ci_obj)) {
      auc_lo <- NA_real_; auc_hi <- NA_real_
    } else {
      auc_lo <- ci_obj[1]; auc_hi <- ci_obj[3]
    }
    message(sprintf("[%s] [%s] AUC=%.4f (%.4f-%.4f DeLong)",
                    SCRIPT_NAME, cohort, auc_val, auc_lo, auc_hi))
  }

  # ------- Adjusted logistic (age + ER if available) -------------------------
  adj_vars  <- c("age", "er_status")
  adj_avail <- adj_vars[sapply(adj_vars, function(v) {
    v %in% names(df) && mean(!is.na(df[[v]])) >= 0.8
  })]

  or_adj   <- NA_real_; or_adj_lo <- NA_real_; or_adj_hi <- NA_real_
  p_adj    <- NA_real_; adj_used  <- "none"

  if (length(adj_avail) > 0) {
    df_adj <- df[complete.cases(df[, c("pcr", "score_z", adj_avail)]), ]
    fml_adj <- as.formula(paste("pcr ~ score_z +", paste(adj_avail, collapse = " + ")))
    old_warn <- getOption("warn"); options(warn = 0)
    fit_adj <- tryCatch(
      glm(fml_adj, data = df_adj, family = binomial),
      error   = function(e) NULL,
      warning = function(w) glm(fml_adj, data = df_adj, family = binomial)
    )
    options(warn = old_warn)
    if (!is.null(fit_adj)) {
      sm_adj  <- summary(fit_adj)
      or_adj  <- exp(coef(fit_adj)[["score_z"]])
      se_adj  <- sm_adj$coefficients["score_z", "Std. Error"]
      or_adj_lo <- exp(coef(fit_adj)[["score_z"]] - 1.96 * se_adj)
      or_adj_hi <- exp(coef(fit_adj)[["score_z"]] + 1.96 * se_adj)
      p_adj   <- sm_adj$coefficients["score_z", "Pr(>|z|)"]
      adj_used <- paste(adj_avail, collapse = "+")
      message(sprintf("[%s] [%s] Adj OR=%.3f (%.3f-%.3f) p=%.4g [vars: %s]",
                      SCRIPT_NAME, cohort, or_adj, or_adj_lo, or_adj_hi, p_adj, adj_used))
    }
  }

  # ------- Collect result row ------------------------------------------------
  results_list[[cohort]] <- tibble(
    cohort         = cohort,
    n_total        = n_total,
    n_pcr1         = n_pcr1,
    pcr_rate       = round(pcr_rate, 4),
    or_uni         = round(or_uni, 4),
    or_uni_lo_wald = round(or_lo_wald, 4),
    or_uni_hi_wald = round(or_hi_wald, 4),
    or_uni_lo_boot = round(boot_ci$ci_low, 4),
    or_uni_hi_boot = round(boot_ci$ci_high, 4),
    p_uni          = signif(p_uni, 4),
    logOR_uni      = round(log(or_uni), 6),
    se_logOR_uni   = round(se_logOR, 6),
    auc            = round(auc_val, 4),
    auc_lo95       = round(auc_lo, 4),
    auc_hi95       = round(auc_hi, 4),
    or_adj         = round(or_adj, 4),
    or_adj_lo95    = round(or_adj_lo, 4),
    or_adj_hi95    = round(or_adj_hi, 4),
    p_adj          = signif(p_adj, 4),
    adj_vars_used  = adj_used
  )
}

# --------------------------------------------------------------------------
# Save combined results
# --------------------------------------------------------------------------
if (length(results_list) == 0) stop(sprintf("[%s] No cohort results generated.", SCRIPT_NAME))

res_df <- bind_rows(results_list)

dir.create(PATHS$results$pcr, showWarnings = FALSE, recursive = TRUE)
readr::write_csv(res_df, out_path)
h  <- sha256_file(out_path)
sz <- file.info(out_path)$size / 1e6

# Register per-cohort
for (cohort in res_df$cohort) {
  registry_append(cohort, "pcr_logistic_results", out_path, h, "ok", SCRIPT_NAME, sz)
}

message(sprintf("\n[%s] === RESULTS SUMMARY ===", SCRIPT_NAME))
print(res_df[, c("cohort", "n_total", "n_pcr1", "or_uni", "or_uni_lo_boot",
                  "or_uni_hi_boot", "p_uni", "auc")])
message(sprintf("[%s] COMPLETED — saved: %s", SCRIPT_NAME, out_path))
