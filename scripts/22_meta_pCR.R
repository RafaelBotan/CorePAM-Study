# =============================================================================
# SCRIPT: 22_meta_pCR.R
# PURPOSE: Random-effects meta-analysis of CorePAM log(OR) across pCR cohorts.
#          Uses DerSimonian-Laird estimator (metafor::rma) on log(OR) ± SE.
#          Outputs: results/pcr/meta_pCR_results.csv
#                   results/pcr/meta_pCR_cohort_weights.csv
#
# PROJECT: Core-PAM (Memorial v6.1 / Freeze Core-PAM — pCR block)
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "22_meta_pCR.R"

suppressPackageStartupMessages({
  library(metafor)
})

FORCE <- as.logical(Sys.getenv("FORCE_RERUN", "FALSE"))

out_meta    <- file.path(PATHS$results$pcr, "meta_pCR_results.csv")
out_cohort  <- file.path(PATHS$results$pcr, "meta_pCR_cohort_weights.csv")
out_meta_or <- file.path(PATHS$results$pcr, "pcr_meta_OR.csv")   # canonical name per §11.9

if (!FORCE && file.exists(out_meta) && file.exists(out_cohort) && file.exists(out_meta_or)) {
  message(sprintf("[%s] Outputs exist — skipping. Set FORCE_RERUN=TRUE to rerun.",
                  SCRIPT_NAME))
  quit(save = "no", status = 0)
}

# --------------------------------------------------------------------------
# 1) Load per-cohort pCR logistic results
# --------------------------------------------------------------------------
pcr_res_path <- file.path(PATHS$results$pcr, "pCR_results_by_cohort.csv")
if (!file.exists(pcr_res_path)) {
  stop(sprintf("[%s] pCR_results_by_cohort.csv not found. Run 21_pCR_logistic_analysis.R first.",
               SCRIPT_NAME))
}

pcr_df <- strict_csv(pcr_res_path)
message(sprintf("[%s] Loaded pCR results: %d cohorts", SCRIPT_NAME, nrow(pcr_df)))
print(pcr_df[, c("cohort", "n_total", "n_pcr1", "logOR_uni", "se_logOR_uni", "or_uni", "p_uni")])

# Filter to cohorts with valid estimates
pcr_df <- pcr_df[!is.na(pcr_df$logOR_uni) & !is.na(pcr_df$se_logOR_uni) &
                 pcr_df$se_logOR_uni > 0, ]
message(sprintf("[%s] Cohorts for meta-analysis: %d", SCRIPT_NAME, nrow(pcr_df)))

if (nrow(pcr_df) < 2) {
  stop(sprintf("[%s] At least 2 cohorts required for meta-analysis. Found: %d",
               SCRIPT_NAME, nrow(pcr_df)))
}

# --------------------------------------------------------------------------
# 2) Random-effects meta-analysis (DerSimonian-Laird)
# --------------------------------------------------------------------------
message(sprintf("[%s] Fitting RE meta-analysis (DerSimonian-Laird)...", SCRIPT_NAME))

old_warn <- getOption("warn"); options(warn = 0)
meta_re <- metafor::rma(
  yi   = pcr_df$logOR_uni,
  sei  = pcr_df$se_logOR_uni,
  method = "DL",
  slab = pcr_df$cohort
)
options(warn = old_warn)

message(sprintf("[%s] RE pooled log(OR) = %.4f (SE=%.4f)",
                SCRIPT_NAME, meta_re$beta[1], meta_re$se[1]))
message(sprintf("[%s] RE pooled OR = %.3f (95%% CI: %.3f–%.3f) p=%.4g",
                SCRIPT_NAME, exp(meta_re$beta[1]),
                exp(meta_re$ci.lb), exp(meta_re$ci.ub),
                meta_re$pval))
message(sprintf("[%s] I²=%.1f%% | τ²=%.4f | Q=%.2f (df=%d) p_het=%.4g",
                SCRIPT_NAME, meta_re$I2, meta_re$tau2,
                meta_re$QE, meta_re$k - 1, meta_re$QEp))

# --------------------------------------------------------------------------
# 3) Fixed-effects meta-analysis (IVW) — for sensitivity comparison
# --------------------------------------------------------------------------
message(sprintf("[%s] Fitting FE meta-analysis (IVW)...", SCRIPT_NAME))

old_warn <- getOption("warn"); options(warn = 0)
meta_fe <- metafor::rma(
  yi     = pcr_df$logOR_uni,
  sei    = pcr_df$se_logOR_uni,
  method = "FE",
  slab   = pcr_df$cohort
)
options(warn = old_warn)

message(sprintf("[%s] FE pooled OR = %.3f (95%% CI: %.3f–%.3f)",
                SCRIPT_NAME, exp(meta_fe$beta[1]),
                exp(meta_fe$ci.lb), exp(meta_fe$ci.ub)))

# --------------------------------------------------------------------------
# 4) Cohort-level weights (RE model)
# --------------------------------------------------------------------------
wi        <- 1 / (pcr_df$se_logOR_uni^2 + meta_re$tau2)
wi_pct    <- 100 * wi / sum(wi)

cohort_df <- tibble(
  cohort        = pcr_df$cohort,
  n_total       = pcr_df$n_total,
  n_pcr1        = pcr_df$n_pcr1,
  or_uni        = round(pcr_df$or_uni, 4),
  or_lo95       = round(exp(pcr_df$logOR_uni - 1.96 * pcr_df$se_logOR_uni), 4),
  or_hi95       = round(exp(pcr_df$logOR_uni + 1.96 * pcr_df$se_logOR_uni), 4),
  p_uni         = signif(pcr_df$p_uni, 4),
  logOR         = round(pcr_df$logOR_uni, 6),
  se_logOR      = round(pcr_df$se_logOR_uni, 6),
  weight_re_pct = round(wi_pct, 2)
)

# --------------------------------------------------------------------------
# 5) Summary meta results table
# --------------------------------------------------------------------------
meta_df <- tibble(
  model              = c("RE", "FE"),
  estimator          = c("DerSimonian-Laird", "IVW"),
  k_cohorts          = nrow(pcr_df),
  pooled_logOR       = round(c(meta_re$beta[1], meta_fe$beta[1]), 6),
  pooled_se_logOR    = round(c(meta_re$se[1],   meta_fe$se[1]),   6),
  pooled_OR          = round(c(exp(meta_re$beta[1]), exp(meta_fe$beta[1])), 4),
  pooled_OR_lo95     = round(c(exp(meta_re$ci.lb),   exp(meta_fe$ci.lb)), 4),
  pooled_OR_hi95     = round(c(exp(meta_re$ci.ub),   exp(meta_fe$ci.ub)), 4),
  p_pooled           = signif(c(meta_re$pval, meta_fe$pval), 4),
  I2_pct             = round(c(meta_re$I2, NA_real_), 2),
  tau2               = round(c(meta_re$tau2, NA_real_), 6),
  Q_stat             = round(c(meta_re$QE, meta_re$QE), 3),
  Q_df               = c(meta_re$k - 1, meta_re$k - 1),
  p_heterogeneity    = signif(c(meta_re$QEp, meta_re$QEp), 4)
)

# --------------------------------------------------------------------------
# 6) Save outputs
# --------------------------------------------------------------------------
dir.create(PATHS$results$pcr, showWarnings = FALSE, recursive = TRUE)

readr::write_csv(meta_df, out_meta)
h1 <- sha256_file(out_meta)
registry_append("META_PCR", "meta_pcr_results", out_meta, h1, "ok", SCRIPT_NAME,
                file.info(out_meta)$size / 1e6)

readr::write_csv(cohort_df, out_cohort)
h2 <- sha256_file(out_cohort)
registry_append("META_PCR", "meta_pcr_cohort_weights", out_cohort, h2, "ok", SCRIPT_NAME,
                file.info(out_cohort)$size / 1e6)

# Save pcr_meta_OR.csv (canonical filename per Memorial §11.9)
readr::write_csv(meta_df, out_meta_or)
h3 <- sha256_file(out_meta_or)
registry_append("META_PCR", "pcr_meta_OR", out_meta_or, h3, "ok", SCRIPT_NAME,
                file.info(out_meta_or)$size / 1e6)

message(sprintf("[%s] Saved meta results: %s", SCRIPT_NAME, out_meta))
message(sprintf("[%s] Saved pcr_meta_OR: %s", SCRIPT_NAME, out_meta_or))
message(sprintf("[%s] Saved cohort weights: %s", SCRIPT_NAME, out_cohort))
message(sprintf("[%s] COMPLETED | RE OR=%.3f (%.3f-%.3f) I²=%.1f%% τ²=%.4f",
                SCRIPT_NAME,
                exp(meta_re$beta[1]), exp(meta_re$ci.lb), exp(meta_re$ci.ub),
                meta_re$I2, meta_re$tau2))
