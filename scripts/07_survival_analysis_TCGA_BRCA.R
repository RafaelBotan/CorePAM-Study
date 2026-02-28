# =============================================================================
# SCRIPT: 07_survival_analysis_TCGA_BRCA.R
# PURPOSE: Analise de sobrevida CorePAM na coorte TCGA-BRCA:
#          Cox univariado + CORE-A, C-index bootstrap, KM (mediana + quartis),
#          sensibilidade horizonte 24m. Segue Memorial v6.1 §8 + §9.4.
# COORTE:  TCGA_BRCA | Endpoint primario: OS
# PROJETO: Core-PAM (Memorial v6.1 / Freeze Core-PAM)
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "07_survival_analysis_TCGA_BRCA.R"
COHORT      <- "TCGA_BRCA"
ENDPOINT    <- "OS"

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
})

set.seed(FREEZE$seed_folds)
message(sprintf("[%s] Iniciando analise de sobrevida %s | %s", SCRIPT_NAME, COHORT, ENDPOINT))

HORIZON_SENS <- as.integer(FREEZE$horizon_sensitivity_months)  # 24 meses

# --------------------------------------------------------------------------
# Helpers internos
# --------------------------------------------------------------------------
bootstrap_cindex <- function(time, event, score_z, n_boot = 1000, seed = 42) {
  set.seed(seed)
  n <- length(time)
  cvals <- numeric(n_boot)
  old_warn <- getOption("warn"); options(warn = 0)
  for (i in seq_len(n_boot)) {
    idx <- sample(n, n, replace = TRUE)
    cx  <- tryCatch(
      concordance(Surv(time[idx], event[idx]) ~ score_z[idx])$concordance,
      error = function(e) NA_real_
    )
    cvals[i] <- cx
  }
  options(warn = old_warn)
  list(
    c_index = concordance(Surv(time, event) ~ score_z)$concordance,
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
# 1) Carregar dados
# --------------------------------------------------------------------------
ready_path <- file.path(proc_cohort(COHORT), "analysis_ready.parquet")
df         <- strict_parquet(ready_path)

time_col  <- "os_time"
event_col <- "os_event"
stopifnot(all(c(time_col, event_col, "score_z") %in% names(df)))

df <- df[!is.na(df[[time_col]]) & df[[time_col]] > 0 &
           !is.na(df[[event_col]]) & !is.na(df$score_z), ]
message(sprintf("[%s] N analise: %d | Eventos: %d", SCRIPT_NAME, nrow(df), sum(df[[event_col]])))

# --------------------------------------------------------------------------
# 2) Cox univariado
# --------------------------------------------------------------------------
old_warn <- getOption("warn"); options(warn = 0)
cox_uni <- coxph(Surv(df[[time_col]], df[[event_col]]) ~ score_z, data = df)
options(warn = old_warn)

sm_uni  <- summary(cox_uni)
hr_uni  <- sm_uni$conf.int[1, "exp(coef)"]
lo_uni  <- sm_uni$conf.int[1, "lower .95"]
hi_uni  <- sm_uni$conf.int[1, "upper .95"]
p_uni   <- sm_uni$coefficients[1, "Pr(>|z|)"]
message(sprintf("[%s] Cox uni HR=%.3f (%.3f-%.3f) p=%.4g", SCRIPT_NAME, hr_uni, lo_uni, hi_uni, p_uni))

# --------------------------------------------------------------------------
# 3) Cox multivariado CORE-A + score_z
# --------------------------------------------------------------------------
corea_vars  <- c("age", "er_status")
corea_avail <- intersect(corea_vars, names(df))
corea_used  <- paste(corea_avail, collapse = "+")

hr_multi <- lo_multi <- hi_multi <- p_multi <- NA_real_

if (length(corea_avail) > 0) {
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
# 4) C-index bootstrap
# --------------------------------------------------------------------------
message(sprintf("[%s] Bootstrap C-index (n=%d)...", SCRIPT_NAME, FREEZE$bootstrap_n))
boot_res <- bootstrap_cindex(
  df[[time_col]], df[[event_col]], df$score_z,
  n_boot = as.integer(FREEZE$bootstrap_n),
  seed   = as.integer(FREEZE$seed_folds)
)
message(sprintf("[%s] C-index=%.4f (%.4f-%.4f)", SCRIPT_NAME,
                boot_res$c_index, boot_res$ci_low, boot_res$ci_high))

fu_median <- reverse_km_median(df[[time_col]], df[[event_col]])
message(sprintf("[%s] Follow-up mediana (Reverse KM): %.1f meses", SCRIPT_NAME, fu_median))

# --------------------------------------------------------------------------
# 5) KM: mediana intra-coorte (principal)
# --------------------------------------------------------------------------
fig_dir <- PATHS$figures$main
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

km_cutpoint <- median(df$score_z, na.rm = TRUE)
df$risk_group_median <- ifelse(df$score_z >= km_cutpoint, "High", "Low")

old_warn <- getOption("warn"); options(warn = 0)
km_fit_med <- survfit(
  Surv(df[[time_col]], df[[event_col]]) ~ risk_group_median, data = df
)
km_plot_med <- ggsurvplot(
  km_fit_med,
  data        = df,
  risk.table  = TRUE,
  pval        = TRUE,
  conf.int    = TRUE,
  palette     = c("#E74C3C", "#2980B9"),
  title       = sprintf("KM CorePAM — %s | %s | Mediana cutpoint", COHORT, ENDPOINT),
  xlab        = "Tempo (meses)",
  ylab        = "Sobrevida geral",
  legend.labs = c("High", "Low"),
  ggtheme     = theme_classic()
)
options(warn = old_warn)

km_pdf <- file.path(fig_dir, sprintf("Fig3_KM_%s_%s_CorePAM.pdf", COHORT, ENDPOINT))
km_png <- file.path(fig_dir, sprintf("Fig3_KM_%s_%s_CorePAM.png", COHORT, ENDPOINT))

old_warn <- getOption("warn"); options(warn = 0)
pdf(km_pdf, width = 8, height = 6); print(km_plot_med); dev.off()
png(km_png, width = 800, height = 600, res = 100); print(km_plot_med); dev.off()
options(warn = old_warn)

registry_append(COHORT, "figure_km_main", km_pdf, sha256_file(km_pdf), "ok", SCRIPT_NAME,
                file.info(km_pdf)$size / 1e6)
registry_append(COHORT, "figure_km_main_png", km_png, sha256_file(km_png), "ok", SCRIPT_NAME,
                file.info(km_png)$size / 1e6)

# Quartis (sensibilidade)
quart_cuts <- quantile(df$score_z, probs = c(0.25, 0.75), na.rm = TRUE)
df$risk_group_quartile <- cut(df$score_z,
                               breaks = c(-Inf, quart_cuts[1], quart_cuts[2], Inf),
                               labels = c("Q1_Low", "Q2_Mid", "Q3_High"))

supp_fig_dir <- PATHS$figures$supp
dir.create(supp_fig_dir, showWarnings = FALSE, recursive = TRUE)

old_warn <- getOption("warn"); options(warn = 0)
km_fit_q <- survfit(
  Surv(df[[time_col]], df[[event_col]]) ~ risk_group_quartile, data = df
)
km_plot_q <- ggsurvplot(
  km_fit_q, data = df, risk.table = TRUE, pval = TRUE, conf.int = FALSE,
  title = sprintf("KM CorePAM — %s | %s | Quartis (sensibilidade)", COHORT, ENDPOINT),
  xlab = "Tempo (meses)", ylab = "Sobrevida geral", ggtheme = theme_classic()
)
options(warn = old_warn)

km_q_pdf <- file.path(supp_fig_dir, sprintf("FigS_KM_%s_%s_Quartis_CorePAM.pdf", COHORT, ENDPOINT))
old_warn <- getOption("warn"); options(warn = 0)
pdf(km_q_pdf, width = 8, height = 6); print(km_plot_q); dev.off()
options(warn = old_warn)
registry_append(COHORT, "figure_km_quartile", km_q_pdf, sha256_file(km_q_pdf), "ok",
                SCRIPT_NAME, file.info(km_q_pdf)$size / 1e6)

# --------------------------------------------------------------------------
# 6) Sensibilidade horizonte 24m (TCGA — follow-up curto; §9.4)
# --------------------------------------------------------------------------
message(sprintf("[%s] Sensibilidade: horizonte %dm", SCRIPT_NAME, HORIZON_SENS))

df_h <- df
df_h$os_time_h  <- pmin(df_h[[time_col]], HORIZON_SENS)
df_h$os_event_h <- ifelse(df_h[[time_col]] <= HORIZON_SENS, df_h[[event_col]], 0)
df_h <- df_h[df_h$os_time_h > 0, ]

old_warn <- getOption("warn"); options(warn = 0)
cox_h <- tryCatch(
  coxph(Surv(os_time_h, os_event_h) ~ score_z, data = df_h),
  error = function(e) NULL
)
options(warn = old_warn)

hr_h <- lo_h <- hi_h <- p_h <- NA_real_
if (!is.null(cox_h)) {
  sm_h  <- summary(cox_h)
  hr_h  <- sm_h$conf.int[1, "exp(coef)"]
  lo_h  <- sm_h$conf.int[1, "lower .95"]
  hi_h  <- sm_h$conf.int[1, "upper .95"]
  p_h   <- sm_h$coefficients[1, "Pr(>|z|)"]
  message(sprintf("[%s] Horizonte %dm: HR=%.3f (%.3f-%.3f) p=%.4g",
                  SCRIPT_NAME, HORIZON_SENS, hr_h, lo_h, hi_h, p_h))
}

# KM sensibilidade 24m
old_warn <- getOption("warn"); options(warn = 0)
km_fit_h <- survfit(Surv(os_time_h, os_event_h) ~ risk_group_median, data = df_h)
km_plot_h <- ggsurvplot(
  km_fit_h, data = df_h, risk.table = TRUE, pval = TRUE, conf.int = TRUE,
  title      = sprintf("KM CorePAM — %s | OS %dm (sens)", COHORT, HORIZON_SENS),
  xlab       = "Tempo (meses)", ylab = "Sobrevida geral",
  xlim       = c(0, HORIZON_SENS),
  palette    = c("#E74C3C", "#2980B9"),
  legend.labs = c("High", "Low"),
  ggtheme    = theme_classic()
)
options(warn = old_warn)

km_h_pdf <- file.path(supp_fig_dir, sprintf("FigS6_TCGA_%dm_Sensitivity.pdf", HORIZON_SENS))
km_h_png <- file.path(supp_fig_dir, sprintf("FigS6_TCGA_%dm_Sensitivity.png", HORIZON_SENS))

old_warn <- getOption("warn"); options(warn = 0)
pdf(km_h_pdf, width = 8, height = 6); print(km_plot_h); dev.off()
png(km_h_png, width = 800, height = 600, res = 100); print(km_plot_h); dev.off()
options(warn = old_warn)

registry_append(COHORT, "figure_km_sensitivity_24m", km_h_pdf, sha256_file(km_h_pdf), "ok",
                SCRIPT_NAME, file.info(km_h_pdf)$size / 1e6)
registry_append(COHORT, "figure_km_sensitivity_24m_png", km_h_png, sha256_file(km_h_png), "ok",
                SCRIPT_NAME, file.info(km_h_png)$size / 1e6)

# --------------------------------------------------------------------------
# 7) Output resultados
# --------------------------------------------------------------------------
supp_dir <- PATHS$results$supp
dir.create(supp_dir, showWarnings = FALSE, recursive = TRUE)

res_df <- tibble(
  cohort              = COHORT,
  endpoint            = ENDPOINT,
  n_samples           = nrow(df),
  n_events            = sum(df[[event_col]]),
  fu_median_months    = round(fu_median, 1),
  hr_uni              = round(hr_uni, 4),
  hr_uni_lo95         = round(lo_uni, 4),
  hr_uni_hi95         = round(hi_uni, 4),
  p_uni               = signif(p_uni, 4),
  hr_multi            = round(hr_multi, 4),
  hr_multi_lo95       = round(lo_multi, 4),
  hr_multi_hi95       = round(hi_multi, 4),
  p_multi             = signif(p_multi, 4),
  corea_vars_used     = corea_used,
  c_index             = round(boot_res$c_index, 4),
  c_index_lo95        = round(boot_res$ci_low, 4),
  c_index_hi95        = round(boot_res$ci_high, 4),
  loghr_uni           = round(log(hr_uni), 6),
  se_loghr_uni        = round((log(hi_uni) - log(lo_uni)) / (2 * 1.96), 6),
  hr_sens_24m         = round(hr_h, 4),
  hr_sens_24m_lo95    = round(lo_h, 4),
  hr_sens_24m_hi95    = round(hi_h, 4),
  p_sens_24m          = signif(p_h, 4)
)

supp_path <- file.path(supp_dir, sprintf("survival_results_%s.csv", COHORT))
readr::write_csv(res_df, supp_path)
h_res <- sha256_file(supp_path)
registry_append(COHORT, "survival_results", supp_path, h_res, "ok", SCRIPT_NAME,
                file.info(supp_path)$size / 1e6)

message(sprintf("[%s] CONCLUIDO para %s | %s", SCRIPT_NAME, COHORT, ENDPOINT))
