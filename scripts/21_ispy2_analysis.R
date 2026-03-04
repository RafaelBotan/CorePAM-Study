# =============================================================================
# SCRIPT: 21_ispy2_analysis.R
# PURPOSE: CorePAM score association with pCR in I-SPY2 (GSE194040).
#          Exploratory External Cohort — supplementary analysis only.
#
# Three pre-specified analyses:
#   (i)   Univariate:       pCR ~ score_z
#   (ii)  Clinically adj.:  pCR ~ score_z + hr + her2 + arm  (arm = categorical)
#   (iii) Control-only:     pCR ~ score_z  (subset: arm == "Paclitaxel")
#
# OUTPUTS (all in results/pcr/):
#   ispy2_results.csv            — full results (3 analyses + AUC + bootstrap CI)
#   ispy2_validation_table.csv   — pcr_validation_table.csv-compatible row
#   ispy2_pcr_meta_input.csv     — logOR + SE for meta-analysis integration
#
# NOTE: No CORE-A / CORE-NACT baseline — age NOT available in GSE194040 metadata.
#       Adjusted model uses HR+HER2+arm as the clinical baseline.
#
# PROJECT: Core-PAM (Memorial v6.1 / Freeze Core-PAM — pCR block)
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "21_ispy2_analysis.R"
COHORT      <- "ISPY2"

suppressPackageStartupMessages(library(pROC))
set.seed(FREEZE$seed_folds)

FORCE <- as.logical(Sys.getenv("FORCE_RERUN", "FALSE"))

out_full  <- file.path(PATHS$results$pcr, "ispy2_results.csv")
out_val   <- file.path(PATHS$results$pcr, "ispy2_validation_table.csv")
out_meta  <- file.path(PATHS$results$pcr, "ispy2_pcr_meta_input.csv")

if (!FORCE && all(file.exists(c(out_full, out_val, out_meta)))) {
  message(sprintf("[%s] Outputs exist — skipping. Set FORCE_RERUN=TRUE to rerun.",
                  SCRIPT_NAME))
  quit(save = "no", status = 0)
}

# --------------------------------------------------------------------------
# 1) Load analysis-ready data
# --------------------------------------------------------------------------
ready_path <- file.path(proc_pcr_cohort(COHORT), "analysis_ready.parquet")
if (!file.exists(ready_path))
  stop(sprintf("[%s] analysis_ready.parquet not found. Run 20_prepare_pCR_ISPY2.R first.",
               SCRIPT_NAME))

old_warn <- getOption("warn"); options(warn = 0)
df <- arrow::read_parquet(ready_path)
options(warn = old_warn)

message(sprintf("[%s] Loaded: %d samples | pCR+=%d (%.1f%%)",
                SCRIPT_NAME, nrow(df),
                sum(df$pcr == 1L, na.rm = TRUE),
                100 * mean(df$pcr == 1L, na.rm = TRUE)))

# --------------------------------------------------------------------------
# 2) Clean arm names
#    Normalize hyphen vs space variants (e.g., AMG-386 → AMG 386)
#    and collapse exact duplicates after trimming.
# --------------------------------------------------------------------------
df$arm_clean <- trimws(gsub("\\s+", " ", gsub("-", " ", df$arm)))
# Restore "T-DM1" branding (Trastuzumab emtansine — should keep hyphen)
df$arm_clean <- gsub("T DM1", "T-DM1", df$arm_clean)

arm_tab <- sort(table(df$arm_clean), decreasing = TRUE)
message(sprintf("[%s] Arm distribution (%d levels):", SCRIPT_NAME, length(arm_tab)))
for (a in names(arm_tab)) message(sprintf("  %-55s n=%d", a, arm_tab[a]))

# Reference arm = Paclitaxel (control)
df$arm_f <- relevel(factor(df$arm_clean), ref = "Paclitaxel")
df$hr_f   <- factor(df$hr,   levels = c(0L, 1L), labels = c("HR-", "HR+"))
df$her2_f <- factor(df$her2, levels = c(0L, 1L), labels = c("HER2-", "HER2+"))

message(sprintf("[%s] HR+: %d (%.1f%%) | HER2+: %d (%.1f%%) | Paclitaxel control: %d",
                SCRIPT_NAME,
                sum(df$hr == 1L, na.rm = TRUE),
                100 * mean(df$hr == 1L, na.rm = TRUE),
                sum(df$her2 == 1L, na.rm = TRUE),
                100 * mean(df$her2 == 1L, na.rm = TRUE),
                sum(df$arm_clean == "Paclitaxel", na.rm = TRUE)))

# --------------------------------------------------------------------------
# 3) Helpers
# --------------------------------------------------------------------------
bootstrap_or_ci <- function(pcr, score_z, n_boot = 1000, seed = 42) {
  set.seed(seed)
  n <- length(pcr); or_boot <- numeric(n_boot)
  old_warn <- getOption("warn"); options(warn = 0)
  for (i in seq_len(n_boot)) {
    idx   <- sample(n, n, replace = TRUE)
    fit_b <- tryCatch(
      glm(pcr[idx] ~ score_z[idx], family = binomial),
      error   = function(e) NULL,
      warning = function(w) suppressWarnings(glm(pcr[idx] ~ score_z[idx], family = binomial))
    )
    or_boot[i] <- if (!is.null(fit_b)) exp(coef(fit_b)[2]) else NA_real_
  }
  options(warn = old_warn)
  list(ci_lo = quantile(or_boot, 0.025, na.rm = TRUE),
       ci_hi = quantile(or_boot, 0.975, na.rm = TRUE))
}

run_auc <- function(pcr, prob) {
  tryCatch({
    r <- suppressMessages(suppressWarnings(
      pROC::roc(response = as.integer(pcr), predictor = as.numeric(prob))
    ))
    auc_val <- as.numeric(suppressMessages(pROC::auc(r)))
    ci <- tryCatch(
      as.numeric(suppressMessages(pROC::ci.auc(r, conf.level = 0.95))),
      error = function(e) c(NA_real_, auc_val, NA_real_)
    )
    list(auc = auc_val, ci_lo = ci[1], ci_hi = ci[3])
  }, error = function(e) {
    message(sprintf("run_auc error: %s", e$message))
    list(auc = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_)
  })
}

# --------------------------------------------------------------------------
# 4) Analysis (i): Univariate — pCR ~ score_z
# --------------------------------------------------------------------------
message(sprintf("[%s] Analysis (i): Univariate pCR ~ score_z", SCRIPT_NAME))

old_warn <- getOption("warn"); options(warn = 0)
fit_uni  <- glm(pcr ~ score_z, data = df, family = binomial())
options(warn = old_warn)

logOR_uni <- coef(fit_uni)["score_z"]
se_uni    <- sqrt(diag(vcov(fit_uni)))["score_z"]
or_uni    <- exp(logOR_uni)
or_uni_lo <- exp(logOR_uni - 1.96 * se_uni)
or_uni_hi <- exp(logOR_uni + 1.96 * se_uni)
p_uni     <- summary(fit_uni)$coefficients["score_z", "Pr(>|z|)"]

boot_uni  <- bootstrap_or_ci(df$pcr, df$score_z)
# Use model frame to ensure pcr/prob vectors are aligned (excludes NA score_z rows)
pcr_uni   <- as.integer(model.frame(fit_uni)[[1]])
prob_uni  <- predict(fit_uni, type = "response")
auc_uni   <- run_auc(pcr_uni, prob_uni)

message(sprintf("[%s]   OR=%.3f (%.3f–%.3f) p=%.4g AUC=%.3f",
                SCRIPT_NAME, or_uni, or_uni_lo, or_uni_hi, p_uni, auc_uni$auc))

# --------------------------------------------------------------------------
# 5) Analysis (ii): Clinically adjusted — pCR ~ score_z + hr + her2 + arm
# --------------------------------------------------------------------------
message(sprintf("[%s] Analysis (ii): Adjusted pCR ~ score_z + hr + her2 + arm", SCRIPT_NAME))

# Subset to samples with complete hr, her2, arm
df_adj <- df[!is.na(df$hr) & !is.na(df$her2) & !is.na(df$arm_clean), ]
message(sprintf("[%s]   N for adjusted model: %d", SCRIPT_NAME, nrow(df_adj)))

old_warn <- getOption("warn"); options(warn = 0)
fit_adj  <- tryCatch(
  glm(pcr ~ score_z + hr_f + her2_f + arm_f, data = df_adj, family = binomial()),
  warning = function(w) {
    message(sprintf("[%s]   Warning in adjusted fit: %s", SCRIPT_NAME, conditionMessage(w)))
    suppressWarnings(glm(pcr ~ score_z + hr_f + her2_f + arm_f, data = df_adj, family = binomial()))
  }
)
options(warn = old_warn)

logOR_adj <- coef(fit_adj)["score_z"]
se_adj    <- sqrt(diag(vcov(fit_adj)))["score_z"]
or_adj    <- exp(logOR_adj)
or_adj_lo <- exp(logOR_adj - 1.96 * se_adj)
or_adj_hi <- exp(logOR_adj + 1.96 * se_adj)
p_adj     <- summary(fit_adj)$coefficients["score_z", "Pr(>|z|)"]

pcr_adj   <- as.integer(model.frame(fit_adj)[[1]])
prob_adj  <- predict(fit_adj, type = "response")
auc_adj   <- run_auc(pcr_adj, prob_adj)

message(sprintf("[%s]   OR_adj=%.3f (%.3f–%.3f) p=%.4g AUC=%.3f",
                SCRIPT_NAME, or_adj, or_adj_lo, or_adj_hi, p_adj, auc_adj$auc))

# --------------------------------------------------------------------------
# 6) Analysis (iii): Paclitaxel control-arm only — pCR ~ score_z
# --------------------------------------------------------------------------
message(sprintf("[%s] Analysis (iii): Paclitaxel control-only", SCRIPT_NAME))

df_ctrl <- df[df$arm_clean == "Paclitaxel" & !is.na(df$arm_clean), ]
message(sprintf("[%s]   N control arm: %d | pCR+=%d (%.1f%%)",
                SCRIPT_NAME, nrow(df_ctrl),
                sum(df_ctrl$pcr == 1L, na.rm = TRUE),
                100 * mean(df_ctrl$pcr == 1L, na.rm = TRUE)))

or_ctrl <- NA_real_; or_ctrl_lo <- NA_real_; or_ctrl_hi <- NA_real_
p_ctrl  <- NA_real_; auc_ctrl   <- list(auc = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_)
logOR_ctrl <- NA_real_; se_ctrl <- NA_real_

if (nrow(df_ctrl) >= 20) {
  old_warn <- getOption("warn"); options(warn = 0)
  fit_ctrl  <- tryCatch(
    glm(pcr ~ score_z, data = df_ctrl, family = binomial()),
    error   = function(e) NULL,
    warning = function(w) suppressWarnings(glm(pcr ~ score_z, data = df_ctrl, family = binomial()))
  )
  options(warn = old_warn)

  if (!is.null(fit_ctrl)) {
    logOR_ctrl <- coef(fit_ctrl)["score_z"]
    se_ctrl    <- sqrt(diag(vcov(fit_ctrl)))["score_z"]
    or_ctrl    <- exp(logOR_ctrl)
    or_ctrl_lo <- exp(logOR_ctrl - 1.96 * se_ctrl)
    or_ctrl_hi <- exp(logOR_ctrl + 1.96 * se_ctrl)
    p_ctrl     <- summary(fit_ctrl)$coefficients["score_z", "Pr(>|z|)"]
    prob_ctrl  <- predict(fit_ctrl, type = "response")
    auc_ctrl   <- run_auc(df_ctrl$pcr, prob_ctrl)
    message(sprintf("[%s]   OR_ctrl=%.3f (%.3f–%.3f) p=%.4g AUC=%.3f",
                    SCRIPT_NAME, or_ctrl, or_ctrl_lo, or_ctrl_hi, p_ctrl, auc_ctrl$auc))
  }
}

# --------------------------------------------------------------------------
# 7) Save full results
# --------------------------------------------------------------------------
dir.create(PATHS$results$pcr, showWarnings = FALSE, recursive = TRUE)

results_df <- data.frame(
  cohort             = COHORT,
  analysis           = c("univariate", "adjusted_hr_her2_arm", "control_paclitaxel_only"),
  n_total            = c(nrow(df), nrow(df_adj), nrow(df_ctrl)),
  n_pcr1             = c(sum(df$pcr == 1L, na.rm = TRUE),
                         sum(df_adj$pcr == 1L, na.rm = TRUE),
                         sum(df_ctrl$pcr == 1L, na.rm = TRUE)),
  pcr_rate           = c(round(mean(df$pcr == 1L, na.rm = TRUE), 4),
                         round(mean(df_adj$pcr == 1L, na.rm = TRUE), 4),
                         round(mean(df_ctrl$pcr == 1L, na.rm = TRUE), 4)),
  genes_total        = nrow(strict_csv(file.path(PATHS$results$corepam, "CorePAM_weights.csv")) |>
                              dplyr::filter(weight != 0)),
  genes_present      = df$genes_present[1],
  coverage_pct       = round(df$genes_present[1] / 24 * 100, 1),
  logOR              = c(logOR_uni, logOR_adj, logOR_ctrl),
  se_logOR           = c(se_uni, se_adj, se_ctrl),
  or_per_1SD         = c(or_uni, or_adj, or_ctrl),
  or_ci_lo95_wald    = c(or_uni_lo, or_adj_lo, or_ctrl_lo),
  or_ci_hi95_wald    = c(or_uni_hi, or_adj_hi, or_ctrl_hi),
  or_ci_lo95_boot    = c(boot_uni$ci_lo, NA_real_, NA_real_),
  or_ci_hi95_boot    = c(boot_uni$ci_hi, NA_real_, NA_real_),
  p_value            = c(p_uni, p_adj, p_ctrl),
  auc                = c(auc_uni$auc, auc_adj$auc, auc_ctrl$auc),
  auc_ci_lo95        = c(auc_uni$ci_lo, auc_adj$ci_lo, auc_ctrl$ci_lo),
  auc_ci_hi95        = c(auc_uni$ci_hi, auc_adj$ci_hi, auc_ctrl$ci_hi),
  adj_vars           = c("none", "hr+her2+arm", "none_control_arm_only"),
  stringsAsFactors   = FALSE
)

readr::write_csv(results_df, out_full)
message(sprintf("[%s] Saved: %s", SCRIPT_NAME, out_full))

# --------------------------------------------------------------------------
# 8) Validation table row (pcr_validation_table.csv-compatible)
# --------------------------------------------------------------------------
val_row <- data.frame(
  cohort            = COHORT,
  endpoint          = "pCR",
  n_total           = nrow(df),
  n_pcr1            = sum(df$pcr == 1L, na.rm = TRUE),
  n_rd              = sum(df$pcr == 0L, na.rm = TRUE),
  genes_total       = 24L,
  genes_present     = df$genes_present[1],
  genes_missing     = 24L - df$genes_present[1],
  coverage_pct      = round(df$genes_present[1] / 24 * 100, 1),
  or_per_1SD        = or_uni,
  or_ci_lo95        = or_uni_lo,
  or_ci_hi95        = or_uni_hi,
  ci_method         = "Wald (95%)",
  p_value           = p_uni,
  auc               = auc_uni$auc,
  auc_lo95          = auc_uni$ci_lo,
  auc_hi95          = auc_uni$ci_hi,
  baseline_covariates = "none (age not available)",
  role              = "Exploratory External Cohort",
  stringsAsFactors  = FALSE
)
readr::write_csv(val_row, out_val)
message(sprintf("[%s] Saved: %s", SCRIPT_NAME, out_val))

# --------------------------------------------------------------------------
# 9) Meta-analysis input (univariate logOR + SE for pooling)
# --------------------------------------------------------------------------
meta_input <- data.frame(
  cohort      = COHORT,
  logOR_uni   = logOR_uni,
  se_logOR_uni = se_uni,
  or_uni      = or_uni,
  n_total     = nrow(df),
  n_pcr1      = sum(df$pcr == 1L, na.rm = TRUE),
  stringsAsFactors = FALSE
)
readr::write_csv(meta_input, out_meta)
message(sprintf("[%s] Saved: %s", SCRIPT_NAME, out_meta))

# --------------------------------------------------------------------------
# 10) SHA-256 registry
# --------------------------------------------------------------------------
for (f in c(out_full, out_val, out_meta)) {
  h  <- sha256_file(f)
  sz <- file.info(f)$size / 1e6
  registry_append(COHORT, paste0("ispy2_", basename(f)), f, h, "ok", SCRIPT_NAME, sz)
}

message(sprintf("[%s] COMPLETED — I-SPY2 | N=%d | OR_uni=%.3f (p=%.4g) | AUC=%.3f",
                SCRIPT_NAME, nrow(df), or_uni, p_uni, auc_uni$auc))
message(sprintf("[%s] OR_adj=%.3f (p=%.4g) | OR_ctrl=%.3f (p=%.4g)",
                SCRIPT_NAME, or_adj, p_adj, or_ctrl, p_ctrl))
