# =============================================================================
# SCRIPT: 22b_meta_pCR_with_ispy2.R
# PURPOSE: Two-version meta-analysis of CorePAM log(OR) for pCR:
#   Version A: GEO cohorts only (GSE25066, GSE20194, GSE32646, ISPY1)
#              → same as 22_meta_pCR.R but re-labelled as "without_ispy2"
#   Version B: GEO cohorts + I-SPY2 (univariate OR per 1 SD)
#              → "with_ispy2" (exploratory sensitivity)
#
# Outputs (results/pcr/):
#   meta_pCR_without_ispy2.csv    — RE+FE meta, 4 GEO cohorts
#   meta_pCR_with_ispy2.csv       — RE+FE meta, 4 GEO + I-SPY2
#   meta_pCR_cohort_weights_with_ispy2.csv — per-cohort RE weights (5 cohorts)
#
# PROJETO: Core-PAM (Memorial v6.1 / Freeze Core-PAM — pCR block)
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "22b_meta_pCR_with_ispy2.R"

suppressPackageStartupMessages(library(metafor))

FORCE <- as.logical(Sys.getenv("FORCE_RERUN", "FALSE"))

out_no_ispy2   <- file.path(PATHS$results$pcr, "meta_pCR_without_ispy2.csv")
out_with_ispy2 <- file.path(PATHS$results$pcr, "meta_pCR_with_ispy2.csv")
out_weights    <- file.path(PATHS$results$pcr, "meta_pCR_cohort_weights_with_ispy2.csv")

if (!FORCE && all(file.exists(c(out_no_ispy2, out_with_ispy2, out_weights)))) {
  message(sprintf("[%s] Outputs exist — skipping. Set FORCE_RERUN=TRUE to rerun.",
                  SCRIPT_NAME))
  quit(save = "no", status = 0)
}

# --------------------------------------------------------------------------
# 1) Load GEO cohort results (from 22_meta_pCR.R input)
# --------------------------------------------------------------------------
pcr_res_path <- file.path(PATHS$results$pcr, "pCR_results_by_cohort.csv")
if (!file.exists(pcr_res_path))
  stop(sprintf("[%s] pCR_results_by_cohort.csv not found. Run 21_pCR_logistic_analysis.R first.",
               SCRIPT_NAME))

pcr_geo <- strict_csv(pcr_res_path)
pcr_geo <- pcr_geo[!is.na(pcr_geo$logOR_uni) & !is.na(pcr_geo$se_logOR_uni) &
                   pcr_geo$se_logOR_uni > 0, ]
message(sprintf("[%s] GEO cohorts: %d | cohorts: %s",
                SCRIPT_NAME, nrow(pcr_geo), paste(pcr_geo$cohort, collapse = ", ")))

# --------------------------------------------------------------------------
# 2) Load I-SPY2 univariate result
# --------------------------------------------------------------------------
ispy2_path <- file.path(PATHS$results$pcr, "ispy2_pcr_meta_input.csv")
if (!file.exists(ispy2_path))
  stop(sprintf("[%s] ispy2_pcr_meta_input.csv not found. Run 21_ispy2_analysis.R first.",
               SCRIPT_NAME))

ispy2_row <- strict_csv(ispy2_path)
message(sprintf("[%s] I-SPY2: logOR=%.4f SE=%.4f OR=%.3f n=%d",
                SCRIPT_NAME, ispy2_row$logOR_uni, ispy2_row$se_logOR_uni,
                ispy2_row$or_uni, ispy2_row$n_total))

# --------------------------------------------------------------------------
# 3) Helper: run RE + FE meta-analysis and return summary tibble
# --------------------------------------------------------------------------
run_meta <- function(cohort_names, logOR_vals, se_vals, n_total_vals, n_pcr1_vals,
                     version_label) {
  df_in <- data.frame(
    cohort   = cohort_names,
    logOR    = logOR_vals,
    se_logOR = se_vals,
    n_total  = n_total_vals,
    n_pcr1   = n_pcr1_vals,
    stringsAsFactors = FALSE
  )

  message(sprintf("[%s] Meta v=%s | k=%d cohorts", SCRIPT_NAME, version_label, nrow(df_in)))

  old_warn <- getOption("warn"); options(warn = 0)
  meta_re <- metafor::rma(yi = df_in$logOR, sei = df_in$se_logOR,
                          method = "DL", slab = df_in$cohort)
  meta_fe <- metafor::rma(yi = df_in$logOR, sei = df_in$se_logOR,
                          method = "FE", slab = df_in$cohort)
  options(warn = old_warn)

  message(sprintf("[%s]   RE: OR=%.3f (%.3f–%.3f) p=%.4g I²=%.1f%% τ²=%.4f",
                  SCRIPT_NAME,
                  exp(meta_re$beta[1]), exp(meta_re$ci.lb), exp(meta_re$ci.ub),
                  meta_re$pval, meta_re$I2, meta_re$tau2))
  message(sprintf("[%s]   FE: OR=%.3f (%.3f–%.3f)",
                  SCRIPT_NAME,
                  exp(meta_fe$beta[1]), exp(meta_fe$ci.lb), exp(meta_fe$ci.ub)))

  # Per-cohort weights (RE)
  wi     <- 1 / (df_in$se_logOR^2 + meta_re$tau2)
  wi_pct <- 100 * wi / sum(wi)

  weights_df <- tibble(
    version       = version_label,
    cohort        = df_in$cohort,
    n_total       = df_in$n_total,
    n_pcr1        = df_in$n_pcr1,
    or_uni        = round(exp(df_in$logOR), 4),
    or_lo95       = round(exp(df_in$logOR - 1.96 * df_in$se_logOR), 4),
    or_hi95       = round(exp(df_in$logOR + 1.96 * df_in$se_logOR), 4),
    logOR         = round(df_in$logOR, 6),
    se_logOR      = round(df_in$se_logOR, 6),
    weight_re_pct = round(wi_pct, 2)
  )

  # Summary table
  meta_summary <- tibble(
    version            = version_label,
    model              = c("RE", "FE"),
    estimator          = c("DerSimonian-Laird", "IVW"),
    k_cohorts          = nrow(df_in),
    pooled_logOR       = round(c(meta_re$beta[1],     meta_fe$beta[1]),   6),
    pooled_se_logOR    = round(c(meta_re$se[1],        meta_fe$se[1]),     6),
    pooled_OR          = round(c(exp(meta_re$beta[1]), exp(meta_fe$beta[1])), 4),
    pooled_OR_lo95     = round(c(exp(meta_re$ci.lb),   exp(meta_fe$ci.lb)), 4),
    pooled_OR_hi95     = round(c(exp(meta_re$ci.ub),   exp(meta_fe$ci.ub)), 4),
    p_pooled           = signif(c(meta_re$pval,        meta_fe$pval), 4),
    I2_pct             = round(c(meta_re$I2,           NA_real_), 2),
    tau2               = round(c(meta_re$tau2,         NA_real_), 6),
    Q_stat             = round(c(meta_re$QE,           meta_re$QE), 3),
    Q_df               = c(meta_re$k - 1,              meta_re$k - 1),
    p_heterogeneity    = signif(c(meta_re$QEp,         meta_re$QEp), 4)
  )

  list(summary = meta_summary, weights = weights_df,
       re = meta_re, fe = meta_fe)
}

# --------------------------------------------------------------------------
# 4) Version A: GEO cohorts only (without I-SPY2)
# --------------------------------------------------------------------------
res_no_ispy2 <- run_meta(
  cohort_names = pcr_geo$cohort,
  logOR_vals   = pcr_geo$logOR_uni,
  se_vals      = pcr_geo$se_logOR_uni,
  n_total_vals = pcr_geo$n_total,
  n_pcr1_vals  = pcr_geo$n_pcr1,
  version_label = "without_ISPY2"
)

# --------------------------------------------------------------------------
# 5) Version B: GEO cohorts + I-SPY2
# --------------------------------------------------------------------------
res_with_ispy2 <- run_meta(
  cohort_names = c(pcr_geo$cohort, ispy2_row$cohort),
  logOR_vals   = c(pcr_geo$logOR_uni,   ispy2_row$logOR_uni),
  se_vals      = c(pcr_geo$se_logOR_uni, ispy2_row$se_logOR_uni),
  n_total_vals = c(pcr_geo$n_total,     ispy2_row$n_total),
  n_pcr1_vals  = c(pcr_geo$n_pcr1,     ispy2_row$n_pcr1),
  version_label = "with_ISPY2"
)

# --------------------------------------------------------------------------
# 6) Save outputs
# --------------------------------------------------------------------------
dir.create(PATHS$results$pcr, showWarnings = FALSE, recursive = TRUE)

readr::write_csv(res_no_ispy2$summary, out_no_ispy2)
readr::write_csv(res_with_ispy2$summary, out_with_ispy2)

# Combined weights (both versions)
all_weights <- dplyr::bind_rows(res_no_ispy2$weights, res_with_ispy2$weights)
readr::write_csv(all_weights, out_weights)

for (f in c(out_no_ispy2, out_with_ispy2, out_weights)) {
  h  <- sha256_file(f)
  sz <- file.info(f)$size / 1e6
  registry_append("META_PCR_ISPY2", basename(f), f, h, "ok", SCRIPT_NAME, sz)
  message(sprintf("[%s] Saved: %s", SCRIPT_NAME, f))
}

# --------------------------------------------------------------------------
# 7) Print comparison
# --------------------------------------------------------------------------
cat("\n=== META-ANALYSIS COMPARISON ===\n")
cat(sprintf("Without I-SPY2 (k=%d): OR_RE = %.3f (%.3f–%.3f) I²=%.1f%%\n",
            res_no_ispy2$re$k,
            exp(res_no_ispy2$re$beta[1]),
            exp(res_no_ispy2$re$ci.lb), exp(res_no_ispy2$re$ci.ub),
            res_no_ispy2$re$I2))
cat(sprintf("With    I-SPY2 (k=%d): OR_RE = %.3f (%.3f–%.3f) I²=%.1f%%\n",
            res_with_ispy2$re$k,
            exp(res_with_ispy2$re$beta[1]),
            exp(res_with_ispy2$re$ci.lb), exp(res_with_ispy2$re$ci.ub),
            res_with_ispy2$re$I2))

message(sprintf("[%s] COMPLETED.", SCRIPT_NAME))
