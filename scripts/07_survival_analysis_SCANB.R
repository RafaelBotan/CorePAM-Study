# =============================================================================
# SCRIPT: 07_survival_analysis_SCANB.R
# PURPOSE: CorePAM survival analysis in the SCAN-B cohort:
#          Univariate Cox + CORE-A, C-index bootstrap, KM (median + quartiles),
#          60m calibration. Follows Memorial v6.1 sec.8.
# COHORT:  SCANB | Primary endpoint: OS
# PROJECT: Core-PAM (Memorial v6.1 / Freeze Core-PAM)
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "07_survival_analysis_SCANB.R"
COHORT      <- "SCANB"
ENDPOINT    <- "OS"

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
})

set.seed(FREEZE$seed_folds)
message(sprintf("[%s] Starting survival analysis %s | %s", SCRIPT_NAME, COHORT, ENDPOINT))

# --------------------------------------------------------------------------
# Internal helpers
# --------------------------------------------------------------------------
bootstrap_cindex <- function(time, event, score_z, n_boot = 1000, seed = 42) {
  set.seed(seed)
  n <- length(time)
  cvals <- numeric(n_boot)
  old_warn <- getOption("warn"); options(warn = 0)
  for (i in seq_len(n_boot)) {
    idx    <- sample(n, n, replace = TRUE)
    cx_raw <- tryCatch(
      concordance(Surv(time[idx], event[idx]) ~ score_z[idx])$concordance,
      error = function(e) NA_real_
    )
    cvals[i] <- max(cx_raw, 1 - cx_raw, na.rm = TRUE)  # Cadj: sign-invariant
  }
  options(warn = old_warn)
  c_raw <- concordance(Surv(time, event) ~ score_z)$concordance
  list(
    c_index = max(c_raw, 1 - c_raw),
    ci_low  = quantile(cvals, 0.025, na.rm = TRUE),
    ci_high = quantile(cvals, 0.975, na.rm = TRUE)
  )
}

reverse_km_median <- function(time, event) {
  old_warn <- getOption("warn"); options(warn = 0)
  fit <- survfit(Surv(time, 1 - event) ~ 1)
  med <- summary(fit)$table["median"]
  options(warn = old_warn)
  as.numeric(med)
}

# --------------------------------------------------------------------------
# 1) Load data
# --------------------------------------------------------------------------
ready_path <- file.path(proc_cohort(COHORT), "analysis_ready.parquet")
df         <- strict_parquet(ready_path)

time_col  <- "os_time_months"
event_col <- "os_event"
stopifnot(all(c(time_col, event_col, "score_z") %in% names(df)))

df <- df[!is.na(df[[time_col]]) & df[[time_col]] > 0 &
           !is.na(df[[event_col]]) & !is.na(df$score_z), ]
message(sprintf("[%s] N analysis: %d | Events: %d", SCRIPT_NAME, nrow(df), sum(df[[event_col]])))

# --------------------------------------------------------------------------
# 2) Univariate Cox: Surv ~ score_z
# --------------------------------------------------------------------------
old_warn <- getOption("warn"); options(warn = 0)
cox_uni <- coxph(Surv(df[[time_col]], df[[event_col]]) ~ score_z, data = df)
options(warn = old_warn)

sm_uni <- summary(cox_uni)
hr_uni  <- sm_uni$conf.int[1, "exp(coef)"]
lo_uni  <- sm_uni$conf.int[1, "lower .95"]
hi_uni  <- sm_uni$conf.int[1, "upper .95"]
p_uni   <- sm_uni$coefficients[1, "Pr(>|z|)"]
message(sprintf("[%s] Cox uni HR=%.3f (%.3f-%.3f) p=%.4g", SCRIPT_NAME, hr_uni, lo_uni, hi_uni, p_uni))

# --------------------------------------------------------------------------
# 3) Multivariate Cox CORE-A + score_z (Age + ER + score_z)
# --------------------------------------------------------------------------
corea_vars <- c("age", "er_status")
corea_avail <- intersect(corea_vars, names(df))
corea_complete <- length(corea_avail) > 0

hr_multi <- lo_multi <- hi_multi <- p_multi <- NA_real_
corea_used <- paste(corea_avail, collapse = "+")

if (corea_complete) {
  df_m <- df[complete.cases(df[, c(corea_avail, "score_z", time_col, event_col)]), ]
  fml  <- as.formula(paste0("Surv(", time_col, ",", event_col, ") ~ score_z + ",
                             paste(corea_avail, collapse = " + ")))
  old_warn <- getOption("warn"); options(warn = 0)
  cox_multi <- tryCatch(coxph(fml, data = df_m), error = function(e) NULL)
  options(warn = old_warn)
  if (!is.null(cox_multi)) {
    sm_m     <- summary(cox_multi)
    hr_multi <- sm_m$conf.int["score_z", "exp(coef)"]
    lo_multi <- sm_m$conf.int["score_z", "lower .95"]
    hi_multi <- sm_m$conf.int["score_z", "upper .95"]
    p_multi  <- sm_m$coefficients["score_z", "Pr(>|z|)"]
    message(sprintf("[%s] Cox multi (CORE-A) HR=%.3f (%.3f-%.3f) p=%.4g",
                    SCRIPT_NAME, hr_multi, lo_multi, hi_multi, p_multi))
  }
}

# --------------------------------------------------------------------------
# 4) C-index with 95% CI via bootstrap (n=1000)
# --------------------------------------------------------------------------
message(sprintf("[%s] Bootstrap C-index (n=%d)...", SCRIPT_NAME, FREEZE$bootstrap_n))
boot_res <- bootstrap_cindex(
  df[[time_col]], df[[event_col]], df$score_z,
  n_boot = as.integer(FREEZE$bootstrap_n),
  seed   = as.integer(FREEZE$seed_folds)
)
message(sprintf("[%s] C-index=%.4f (%.4f-%.4f)", SCRIPT_NAME,
                boot_res$c_index, boot_res$ci_low, boot_res$ci_high))

# Median follow-up (Reverse KM)
fu_median <- reverse_km_median(df[[time_col]], df[[event_col]])
message(sprintf("[%s] Median follow-up (Reverse KM): %.1f months", SCRIPT_NAME, fu_median))

# --------------------------------------------------------------------------
# 5) KM: intra-cohort median (primary) + quartiles (sensitivity)
# --------------------------------------------------------------------------
# figures/{section}/{lang}/{ext}/ — dirs already created by 00_setup.R

km_cutpoint <- median(df$score_z, na.rm = TRUE)
df$risk_group_median <- ifelse(df$score_z >= km_cutpoint, "High", "Low")

# Cox HR for dichotomized groups (High vs Low) — used for KM annotation
cox_km <- coxph(Surv(df[[time_col]], df[[event_col]]) ~ I(risk_group_median == "High"), data = df)
sm_km  <- summary(cox_km)
hr_km  <- sm_km$conf.int[1, "exp(coef)"]
lo_km  <- sm_km$conf.int[1, "lower .95"]
hi_km  <- sm_km$conf.int[1, "upper .95"]
p_km   <- sm_km$coefficients[1, "Pr(>|z|)"]

old_warn <- getOption("warn"); options(warn = 0)
km_fit_med <- survfit(
  Surv(df[[time_col]], df[[event_col]]) ~ risk_group_median,
  data = df
)
options(warn = old_warn)

km_plot_med <- ggsurvplot(
  km_fit_med,
  data          = df,
  risk.table    = TRUE,
  pval          = FALSE,
  conf.int      = TRUE,
  palette       = c("#E74C3C", "#2980B9"),
  title         = sprintf("KM CorePAM — %s | %s | Intra-cohort median cutpoint (pre-specified, not optimized)", COHORT, ENDPOINT),
  xlab          = "Time (months)",
  ylab          = "Overall survival",
  legend.labs   = c("High risk", "Low risk"),
  ggtheme       = theme_classic()
)

# Add High-vs-Low HR annotation (positioned at 5 months, bottom of plot)
km_hr_lbl <- sprintf("High vs Low: HR = %.2f (%.2f\u2013%.2f), p = %s",
                     hr_km, lo_km, hi_km,
                     formatC(p_km, format = "e", digits = 1))
km_plot_med$plot <- km_plot_med$plot +
  ggplot2::annotate("text",
                    x = 5, y = 0.10,
                    label = km_hr_lbl,
                    hjust = 0, size = 3.0, colour = "grey20")

km_pdf <- file.path(PATHS$figures$main_en_pdf, sprintf("Fig3_KM_%s_%s_CorePAM.pdf", COHORT, ENDPOINT))
km_png <- file.path(PATHS$figures$main_en_png, sprintf("Fig3_KM_%s_%s_CorePAM.png", COHORT, ENDPOINT))

old_warn <- getOption("warn"); options(warn = 0)
pdf(km_pdf, width = 8, height = 6)
print(km_plot_med)
dev.off()
png(km_png, width = 800, height = 600, res = 100)
print(km_plot_med)
dev.off()
options(warn = old_warn)

h_km_pdf <- sha256_file(km_pdf)
h_km_png <- sha256_file(km_png)
registry_append(COHORT, "figure_km_main", km_pdf, h_km_pdf, "ok", SCRIPT_NAME,
                file.info(km_pdf)$size / 1e6)
registry_append(COHORT, "figure_km_main_png", km_png, h_km_png, "ok", SCRIPT_NAME,
                file.info(km_png)$size / 1e6)

# KM quartiles (sensitivity)
quart_cuts <- quantile(df$score_z, probs = c(0.25, 0.75), na.rm = TRUE)
df$risk_group_quartile <- cut(df$score_z,
                               breaks = c(-Inf, quart_cuts[1], quart_cuts[2], Inf),
                               labels = c("Q1_Low", "Q2_Mid", "Q3_High"))

old_warn <- getOption("warn"); options(warn = 0)
km_fit_q <- survfit(
  Surv(df[[time_col]], df[[event_col]]) ~ risk_group_quartile,
  data = df
)
km_plot_q <- ggsurvplot(
  km_fit_q,
  data        = df,
  risk.table  = TRUE,
  pval        = TRUE,
  conf.int    = FALSE,
  title       = sprintf("KM CorePAM — %s | %s | Quartiles (sensitivity)", COHORT, ENDPOINT),
  xlab        = "Time (months)",
  ylab        = "Overall survival",
  ggtheme     = theme_classic()
)
options(warn = old_warn)

km_q_pdf <- file.path(PATHS$figures$supp_en_pdf, sprintf("FigS_KM_%s_%s_Quartis_CorePAM.pdf", COHORT, ENDPOINT))

old_warn <- getOption("warn"); options(warn = 0)
pdf(km_q_pdf, width = 8, height = 6)
print(km_plot_q)
dev.off()
options(warn = old_warn)

h_km_q <- sha256_file(km_q_pdf)
registry_append(COHORT, "figure_km_quartile", km_q_pdf, h_km_q, "ok", SCRIPT_NAME,
                file.info(km_q_pdf)$size / 1e6)

# --------------------------------------------------------------------------
# 6) Output results per cohort
# --------------------------------------------------------------------------
supp_dir <- PATHS$results$supp
dir.create(supp_dir, showWarnings = FALSE, recursive = TRUE)

# Read OOF C-index from model selection summary (derived in script 05)
oof_json_path <- file.path(PATHS$results$corepam, "CorePAM_training_card.json")
c_oof <- c_oof_pam50full <- c_oof_gap <- NA_real_
if (file.exists(oof_json_path)) {
  old_warn <- getOption("warn"); options(warn = 0)
  oof_summary <- tryCatch(jsonlite::fromJSON(oof_json_path), error = function(e) NULL)
  options(warn = old_warn)
  if (!is.null(oof_summary)) {
    c_oof          <- round(oof_summary$selected$c_adj,  4)
    c_oof_pam50full <- round(oof_summary$selected$c_max, 4)
    c_oof_gap       <- round(oof_summary$selected$gap,   4)
    message(sprintf("[%s] OOF C-index: CorePAM=%.4f | PAM50-full=%.4f | gap=%.4f",
                    SCRIPT_NAME, c_oof, c_oof_pam50full, c_oof_gap))
  }
}

res_df <- tibble(
  cohort             = COHORT,
  endpoint           = ENDPOINT,
  n_samples          = nrow(df),
  n_events           = sum(df[[event_col]]),
  fu_median_months   = round(fu_median, 1),
  hr_uni             = round(hr_uni, 4),
  hr_uni_lo95        = round(lo_uni, 4),
  hr_uni_hi95        = round(hi_uni, 4),
  p_uni              = signif(p_uni, 4),
  hr_multi           = round(hr_multi, 4),
  hr_multi_lo95      = round(lo_multi, 4),
  hr_multi_hi95      = round(hi_multi, 4),
  p_multi            = signif(p_multi, 4),
  corea_vars_used    = corea_used,
  c_index            = round(boot_res$c_index, 4),  # in-sample (training, no optimism correction)
  c_index_lo95       = round(boot_res$ci_low, 4),
  c_index_hi95       = round(boot_res$ci_high, 4),
  c_oof              = c_oof,           # OOF C-index used for model selection (script 05)
  c_oof_pam50full    = c_oof_pam50full, # PAM50-full OOF C-index
  c_oof_gap          = c_oof_gap,       # gap = c_oof_pam50full - c_oof (must be < delta_c)
  loghr_uni          = round(log(hr_uni), 6),
  se_loghr_uni       = round((log(hi_uni) - log(lo_uni)) / (2 * 1.96), 6)
)

supp_path <- file.path(supp_dir, sprintf("survival_results_%s.csv", COHORT))
readr::write_csv(res_df, supp_path)
h_res <- sha256_file(supp_path)
registry_append(COHORT, "survival_results", supp_path, h_res, "ok", SCRIPT_NAME,
                file.info(supp_path)$size / 1e6)

message(sprintf("[%s] COMPLETED for %s | %s", SCRIPT_NAME, COHORT, ENDPOINT))
