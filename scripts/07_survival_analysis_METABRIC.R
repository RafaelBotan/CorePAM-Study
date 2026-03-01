# =============================================================================
# SCRIPT: 07_survival_analysis_METABRIC.R
# PURPOSE: CorePAM survival analysis in the METABRIC cohort:
#          Univariate Cox + CORE-A, C-index bootstrap, KM (median + quartiles),
#          DSS (primary) + OS (sensitivity) + Fine-Gray (sensitivity).
#          Follows Memorial v6.1 sec.7.2 + sec.8.
# COHORT:  METABRIC | Primary endpoint: DSS
# PROJETO: Core-PAM (Memorial v6.1 / Freeze Core-PAM)
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "07_survival_analysis_METABRIC.R"
COHORT      <- "METABRIC"
ENDPOINT    <- "DSS"

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
})
source("Y:/Phd-Genomic-claude/scripts/00_colors.R")

# cmprsk for Fine-Gray (sensitivity) — load with warning suppression
old_warn <- getOption("warn"); options(warn = 0)
has_cmprsk <- requireNamespace("cmprsk", quietly = TRUE)
options(warn = old_warn)

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
    cvals[i] <- max(cx_raw, 1 - cx_raw, na.rm = TRUE)  # C_adj: invariant to score sign
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

x_cutoff <- function(km_fit, min_nrisk = 10) {
  s <- summary(km_fit)
  min_by_t <- tapply(s$n.risk, s$time, min)
  valid <- as.numeric(names(min_by_t[min_by_t >= min_nrisk]))
  if (length(valid) == 0) return(max(s$time))
  max(valid)
}

# --------------------------------------------------------------------------
# 1) Load data
# --------------------------------------------------------------------------
ready_path <- file.path(proc_cohort(COHORT), "analysis_ready.parquet")
df         <- strict_parquet(ready_path)

# Primary endpoint: DSS; sensitivity: OS
has_dss <- all(c("dss_time_months", "dss_event") %in% names(df))
if (!has_dss) {
  message(sprintf("[%s] WARNING: DSS not found; using OS as primary endpoint", SCRIPT_NAME))
  time_col  <- "os_time_months"
  event_col <- "os_event"
  ENDPOINT  <- "OS"
} else {
  time_col  <- "dss_time_months"
  event_col <- "dss_event"
}

stopifnot(all(c(time_col, event_col, "score_z") %in% names(df)))

df <- df[!is.na(df[[time_col]]) & df[[time_col]] > 0 &
           !is.na(df[[event_col]]) & !is.na(df$score_z), ]
message(sprintf("[%s] N analysis (%s): %d | Events: %d",
                SCRIPT_NAME, ENDPOINT, nrow(df), sum(df[[event_col]])))

# --------------------------------------------------------------------------
# 2) Univariate Cox (DSS)
# --------------------------------------------------------------------------
old_warn <- getOption("warn"); options(warn = 0)
cox_uni <- coxph(Surv(df[[time_col]], df[[event_col]]) ~ score_z, data = df)
options(warn = old_warn)

sm_uni  <- summary(cox_uni)
hr_uni  <- sm_uni$conf.int[1, "exp(coef)"]
lo_uni  <- sm_uni$conf.int[1, "lower .95"]
hi_uni  <- sm_uni$conf.int[1, "upper .95"]
p_uni   <- sm_uni$coefficients[1, "Pr(>|z|)"]
message(sprintf("[%s] Cox uni (%s) HR=%.3f (%.3f-%.3f) p=%.4g",
                SCRIPT_NAME, ENDPOINT, hr_uni, lo_uni, hi_uni, p_uni))

# --------------------------------------------------------------------------
# 3) Multivariate Cox CORE-A + score_z
# --------------------------------------------------------------------------
corea_vars  <- c("age", "er_status")
corea_avail <- intersect(corea_vars, names(df))
# Only include covariates with ≥80% non-NA (METABRIC er_status only ~22% populated)
corea_avail <- corea_avail[sapply(corea_avail, function(v) mean(!is.na(df[[v]])) >= 0.8)]
corea_used  <- if (length(corea_avail) > 0) paste(corea_avail, collapse = "+") else "none"

hr_multi <- lo_multi <- hi_multi <- p_multi <- NA_real_

if (length(corea_avail) > 0 && corea_used != "none") {
  df_m <- df[complete.cases(df[, c(corea_avail, "score_z", time_col, event_col)]), ]
  message(sprintf("[%s] CORE-A vars: %s | n complete: %d", SCRIPT_NAME, corea_used, nrow(df_m)))
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
    message(sprintf("[%s] Cox multi CORE-A HR=%.3f (%.3f-%.3f) p=%.4g",
                    SCRIPT_NAME, hr_multi, lo_multi, hi_multi, p_multi))
  }
}

# --------------------------------------------------------------------------
# 4) C-index with bootstrap CI
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
message(sprintf("[%s] Median follow-up (Reverse KM): %.1f months", SCRIPT_NAME, fu_median))

# --------------------------------------------------------------------------
# 5) KM: intra-cohort median (primary, DSS) + quartiles (sensitivity)
# --------------------------------------------------------------------------
# figures/{section}/{lang}/{ext}/ — dirs already created by 00_setup.R

km_cutpoint <- median(df$score_z, na.rm = TRUE)
df$risk_group_median <- ifelse(df$score_z >= km_cutpoint, "High", "Low")

old_warn <- getOption("warn"); options(warn = 0)
km_fit_med <- survfit(
  Surv(df[[time_col]], df[[event_col]]) ~ risk_group_median, data = df
)
options(warn = old_warn)
xlim_val  <- x_cutoff(km_fit_med, min_nrisk = 10)
break_val <- max(12L, as.integer(round(xlim_val / 5 / 12) * 12))
message(sprintf("[%s] KM xlim: %.0f months | break: %d months", SCRIPT_NAME, xlim_val, break_val))

# Bilingual KM main — EN and PT
old_warn <- getOption("warn"); options(warn = 0)
for (lang in c("EN", "PT")) {
  lang_lc <- tolower(lang)
  xlb <- if (lang == "EN") "Time (months)"                   else "Tempo (meses)"
  ylb <- if (lang == "EN") "Disease-specific survival"       else "Sobrevida específica"
  lbs <- if (lang == "EN") c("High risk", "Low risk")        else c("Alto risco", "Baixo risco")
  km_p <- ggsurvplot(
    km_fit_med, data = df, risk.table = TRUE, pval = TRUE, conf.int = TRUE,
    palette = c(COL$km_high, COL$km_low),
    xlab = xlb, ylab = ylb, legend.labs = lbs, legend.title = NULL,
    xlim = c(0, xlim_val), break.time.by = break_val,
    pval.coord = c(xlim_val * 0.5, 0.15),
    ggtheme = theme_classic()
  )
  # Add HR/CI annotation and pre-specified cutoff label
  hr_anno_lbl <- if (lang == "EN")
    sprintf("HR = %.2f (%.2f\u2013%.2f), p = %s",
            hr_uni, lo_uni, hi_uni, formatC(p_uni, format = "e", digits = 1))
  else
    sprintf("HR = %.2f (IC: %.2f\u2013%.2f), p = %s",
            hr_uni, lo_uni, hi_uni, formatC(p_uni, format = "e", digits = 1))
  subt_lbl <- if (lang == "EN")
    "Cutoff: intra-cohort median (pre-specified, not optimized)"
  else
    "Ponto de corte: mediana intracoorte (pr\u00e9-especificado, n\u00e3o otimizado)"
  km_p$plot <- km_p$plot +
    ggplot2::labs(subtitle = subt_lbl) +
    ggplot2::annotate("text",
                      x = xlim_val * 0.02, y = 0.10,
                      label = hr_anno_lbl,
                      hjust = 0, size = 3.0, colour = "grey20")
  pdf_path <- file.path(PATHS$figures[[paste0("main_", lang_lc, "_pdf")]],
                        sprintf("Fig2_KM_%s_%s_%s.pdf", COHORT, ENDPOINT, lang))
  png_path <- file.path(PATHS$figures[[paste0("main_", lang_lc, "_png")]],
                        sprintf("Fig2_KM_%s_%s_%s.png", COHORT, ENDPOINT, lang))
  cairo_pdf(pdf_path, width = 8, height = 6); print(km_p); dev.off()
  png(png_path, width = 8, height = 6, units = "in", res = 600); print(km_p); dev.off()
  registry_append(COHORT, sprintf("figure_km_main_%s", lang), pdf_path,
                  sha256_file(pdf_path), "ok", SCRIPT_NAME, file.info(pdf_path)$size / 1e6)
  registry_append(COHORT, sprintf("figure_km_main_png_%s", lang), png_path,
                  sha256_file(png_path), "ok", SCRIPT_NAME, file.info(png_path)$size / 1e6)
  cat(sprintf("[%s] [%s] %s | %s\n", SCRIPT_NAME, lang, pdf_path, png_path))
}
gc(); options(warn = old_warn)

# Quartiles
quart_cuts <- quantile(df$score_z, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
df$risk_group_quartile <- cut(df$score_z,
                               breaks = c(-Inf, quart_cuts[1], quart_cuts[2], quart_cuts[3], Inf),
                               labels = c("Q1", "Q2", "Q3", "Q4"))

old_warn <- getOption("warn"); options(warn = 0)
km_fit_q <- survfit(
  Surv(df[[time_col]], df[[event_col]]) ~ risk_group_quartile, data = df
)
options(warn = old_warn)

old_warn <- getOption("warn"); options(warn = 0)
for (lang in c("EN", "PT")) {
  lang_lc <- tolower(lang)
  xlb <- if (lang == "EN") "Time (months)"                  else "Tempo (meses)"
  ylb <- if (lang == "EN") "Disease-specific survival"      else "Sobrevida específica por doença"
  ttl <- if (lang == "EN")
    sprintf("KM CorePAM — %s | %s | Quartiles (Q1=lowest score/best prognosis, Q4=highest/worst)", COHORT, ENDPOINT)
  else
    sprintf("KM CorePAM — %s | %s | Quartis (Q1=menor escore/melhor prognóstico, Q4=maior/pior)", COHORT, ENDPOINT)
  lbs <- c("Q1", "Q2", "Q3", "Q4")
  km_plot_q <- ggsurvplot(
    km_fit_q, data = df, risk.table = TRUE, pval = TRUE, conf.int = FALSE,
    title = ttl, xlab = xlb, ylab = ylb, legend.labs = lbs,
    ggtheme = theme_classic()
  )
  km_q_pdf <- file.path(PATHS$figures[[paste0("supp_", lang_lc, "_pdf")]],
                        sprintf("FigS_KM_%s_%s_Quartis_CorePAM_%s.pdf", COHORT, ENDPOINT, lang))
  km_q_png <- file.path(PATHS$figures[[paste0("supp_", lang_lc, "_png")]],
                        sprintf("FigS_KM_%s_%s_Quartis_CorePAM_%s.png", COHORT, ENDPOINT, lang))
  cairo_pdf(km_q_pdf, width = 8, height = 6); print(km_plot_q); dev.off()
  png(km_q_png, width = 8, height = 6, units = "in", res = 600); print(km_plot_q); dev.off()
  registry_append(COHORT, sprintf("figure_km_quartile_%s", lang), km_q_pdf, sha256_file(km_q_pdf), "ok",
                  SCRIPT_NAME, file.info(km_q_pdf)$size / 1e6)
  cat(sprintf("[%s] [%s] %s | %s\n", SCRIPT_NAME, lang, km_q_pdf, km_q_png))
}
gc(); options(warn = old_warn)

# --------------------------------------------------------------------------
# 6) Sensitivity OS (METABRIC)
# --------------------------------------------------------------------------
hr_os <- lo_os <- hi_os <- p_os <- NA_real_
c_os  <- ci_os_lo <- ci_os_hi <- NA_real_

if (has_dss && all(c("os_time_months", "os_event") %in% names(df))) {
  df_os <- df[!is.na(df$os_time_months) & df$os_time_months > 0 &
                !is.na(df$os_event) & !is.na(df$score_z), ]
  old_warn <- getOption("warn"); options(warn = 0)
  cox_os  <- tryCatch(
    coxph(Surv(os_time_months, os_event) ~ score_z, data = df_os),
    error = function(e) NULL
  )
  options(warn = old_warn)
  if (!is.null(cox_os)) {
    sm_os  <- summary(cox_os)
    hr_os  <- sm_os$conf.int[1, "exp(coef)"]
    lo_os  <- sm_os$conf.int[1, "lower .95"]
    hi_os  <- sm_os$conf.int[1, "upper .95"]
    p_os   <- sm_os$coefficients[1, "Pr(>|z|)"]
    boot_os <- bootstrap_cindex(
      df_os$os_time_months, df_os$os_event, df_os$score_z,
      n_boot = as.integer(FREEZE$bootstrap_n),
      seed   = as.integer(FREEZE$seed_folds)
    )
    c_os      <- boot_os$c_index
    ci_os_lo  <- boot_os$ci_low
    ci_os_hi  <- boot_os$ci_high
    message(sprintf("[%s] Sensitivity OS: HR=%.3f (%.3f-%.3f) C=%.4f", SCRIPT_NAME, hr_os, lo_os, hi_os, c_os))
  }

  # KM OS sensitivity — bilingual (tryCatch: FigS5 may already exist if S7 GC error)
  old_warn <- getOption("warn"); options(warn = 0)
  km_fit_os <- survfit(Surv(os_time_months, os_event) ~ risk_group_median, data = df_os)
  options(warn = old_warn)
  xlim_os_val  <- x_cutoff(km_fit_os, min_nrisk = 10)
  break_os_val <- max(24L, as.integer(round(xlim_os_val / 5 / 24) * 24))
  message(sprintf("[%s] FigS5 OS xlim: %.0f mo | break: %d mo", SCRIPT_NAME, xlim_os_val, break_os_val))
  old_warn <- getOption("warn"); options(warn = 0)
  tryCatch({
    for (lang in c("EN", "PT")) {
      xlb <- if (lang == "EN") "Time (months)"              else "Tempo (meses)"
      ylb <- if (lang == "EN") "Overall survival (OS)"      else "Sobrevida global (OS)"
      ttl <- if (lang == "EN") "METABRIC — Sensitivity analysis using OS (primary endpoint: DSS)"
             else "METABRIC — Análise de sensibilidade usando OS (desfecho primário: DSS)"
      lbs <- if (lang == "EN") c("High risk", "Low risk") else c("Alto risco", "Baixo risco")
      lang_lc   <- tolower(lang)
      km_os_pdf <- file.path(PATHS$figures[[paste0("supp_", lang_lc, "_pdf")]],
                             sprintf("FigS5_METABRIC_Sensitivity_%s.pdf", lang))
      km_os_png <- file.path(PATHS$figures[[paste0("supp_", lang_lc, "_png")]],
                             sprintf("FigS5_METABRIC_Sensitivity_%s.png", lang))
      km_os_p <- ggsurvplot(
        km_fit_os, data = df_os, risk.table = TRUE, pval = TRUE, conf.int = TRUE,
        palette = c(COL$km_high, COL$km_low),
        xlab = xlb, ylab = ylb, legend.labs = lbs, legend.title = NULL,
        xlim = c(0, xlim_os_val), break.time.by = break_os_val,
        pval.coord = c(xlim_os_val * 0.5, 0.15),
        title = ttl,
        ggtheme = theme_classic()
      )
      cairo_pdf(km_os_pdf, width = 8, height = 6); print(km_os_p); dev.off()
      png(km_os_png, width = 8, height = 6, units = "in", res = 600); print(km_os_p); dev.off()
      if (file.exists(km_os_pdf))
        registry_append(COHORT, sprintf("figure_km_sens_os_%s", lang), km_os_pdf,
                        sha256_file(km_os_pdf), "ok", SCRIPT_NAME, file.info(km_os_pdf)$size / 1e6)
      cat(sprintf("[%s] [%s] %s | %s\n", SCRIPT_NAME, lang, km_os_pdf, km_os_png))
    }
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    message(sprintf("[%s] FigS5 ggsurvplot S7 error (using existing files): %s", SCRIPT_NAME, e$message))
  })
  gc(); options(warn = old_warn)
}

# --------------------------------------------------------------------------
# 7) Sensitivity Fine-Gray (competing risks)
# --------------------------------------------------------------------------
fg_coef <- fg_p <- NA_real_

if (has_cmprsk && has_dss && all(c("os_time_months", "os_event") %in% names(df))) {
  message(sprintf("[%s] Fine-Gray (sensitivity)...", SCRIPT_NAME))
  # Event: 1=cancer death, 2=other cause death, 0=alive/censored
  df$fg_status <- ifelse(df[[event_col]] == 1, 1,
                    ifelse(!is.na(df$os_event) & df$os_event == 1 & df[[event_col]] == 0, 2, 0))
  old_warn <- getOption("warn"); options(warn = 0)
  fg_res <- tryCatch(
    cmprsk::crr(
      ftime   = df[[time_col]],
      fstatus = df$fg_status,
      cov1    = as.matrix(df$score_z),
      failcode = 1,
      cencode  = 0
    ),
    error = function(e) { message("[", SCRIPT_NAME, "] Fine-Gray error:", e$message); NULL }
  )
  options(warn = old_warn)
  if (!is.null(fg_res)) {
    fg_coef <- fg_res$coef[1]
    # p-value: use p.values.2 if available (older cmprsk), else compute Wald p
    fg_p <- tryCatch(
      {
        pv <- fg_res$p.values.2
        if (!is.null(pv) && length(pv) >= 1 && is.numeric(pv)) pv[1]
        else {
          fg_se <- sqrt(fg_res$var[1, 1])
          2 * pnorm(-abs(fg_coef / fg_se))
        }
      },
      error = function(e) {
        fg_se <- sqrt(fg_res$var[1, 1])
        2 * pnorm(-abs(fg_coef / fg_se))
      }
    )
    message(sprintf("[%s] Fine-Gray: coef=%.4f (exp=%.3f) p=%.4g",
                    SCRIPT_NAME, fg_coef, exp(fg_coef), fg_p))
  }
}

# --------------------------------------------------------------------------
# 8) Output results
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
  hr_sens_os          = round(hr_os, 4),
  hr_sens_os_lo95     = round(lo_os, 4),
  hr_sens_os_hi95     = round(hi_os, 4),
  p_sens_os           = signif(p_os, 4),
  c_sens_os           = round(c_os, 4),
  fg_coef             = round(fg_coef, 6),
  fg_p                = signif(fg_p, 4)
)

supp_path <- file.path(supp_dir, sprintf("survival_results_%s.csv", COHORT))
readr::write_csv(res_df, supp_path)
h_res <- sha256_file(supp_path)
registry_append(COHORT, "survival_results", supp_path, h_res, "ok", SCRIPT_NAME,
                file.info(supp_path)$size / 1e6)

message(sprintf("[%s] COMPLETED for %s | %s", SCRIPT_NAME, COHORT, ENDPOINT))
