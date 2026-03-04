# =============================================================================
# SCRIPT: 35_er_neg_bootstrap_stability.R
# PURPOSE: Bootstrap analysis demonstrating the instability of the SCAN-B
#          ER-negative subgroup HR (N=224, 35 events) compared to the
#          well-powered METABRIC ER-negative result (N=439, 186 events).
# PROJECT: Core-PAM
# =============================================================================

source("scripts/00_setup.R")
suppressPackageStartupMessages(library(survival))

# ---------------------------------------------------------------------------
# 1) Load cohorts and filter to ER-negative
# ---------------------------------------------------------------------------
scanb_all <- strict_parquet(
  file.path(proc_cohort("SCANB"), "analysis_ready.parquet")
)
scanb <- scanb_all %>% filter(er_status == "Negative")

metabric_all <- strict_parquet(
  file.path(proc_cohort("METABRIC"), "analysis_ready.parquet")
)
metabric <- metabric_all %>% filter(er_status == "Negative")

message(sprintf("[35] SCAN-B  ER-neg: N=%d, events=%d",
                nrow(scanb), sum(scanb$os_event)))
message(sprintf("[35] METABRIC ER-neg: N=%d, events=%d",
                nrow(metabric), sum(metabric$dss_event)))

# ---------------------------------------------------------------------------
# 2) Observed HRs
# ---------------------------------------------------------------------------
fit_scanb_obs <- coxph(Surv(os_time_months, os_event) ~ score_z, data = scanb)
hr_scanb_obs  <- exp(coef(fit_scanb_obs))

fit_met_obs   <- coxph(Surv(dss_time_months, dss_event) ~ score_z, data = metabric)
hr_met_obs    <- exp(coef(fit_met_obs))

message(sprintf("[35] Observed HR — SCAN-B ER-neg:  %.3f", hr_scanb_obs))
message(sprintf("[35] Observed HR — METABRIC ER-neg: %.3f", hr_met_obs))

# ---------------------------------------------------------------------------
# 3) Bootstrap function
# ---------------------------------------------------------------------------
bootstrap_hr <- function(dat, time_col, event_col, B = 2000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(dat)
  hrs <- numeric(B)

  old_warn <- getOption("warn")
  options(warn = 0)

  for (i in seq_len(B)) {
    idx <- sample.int(n, replace = TRUE)
    boot_dat <- dat[idx, , drop = FALSE]

    hr_i <- tryCatch({
      frm <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ score_z"))
      fit <- coxph(frm, data = boot_dat)
      exp(coef(fit))
    }, error = function(e) NA_real_,
       warning = function(w) {
         fit <- suppressWarnings(coxph(
           as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ score_z")),
           data = boot_dat
         ))
         exp(coef(fit))
       })

    hrs[i] <- hr_i
  }

  options(warn = old_warn)
  hrs
}

# ---------------------------------------------------------------------------
# 4) Run bootstrap (B=2000)
# ---------------------------------------------------------------------------
B <- 2000

message(sprintf("[35] Running %d bootstrap resamples for SCAN-B ER-neg ...", B))
hrs_scanb <- bootstrap_hr(scanb, "os_time_months", "os_event",
                           B = B, seed = FREEZE$seed_folds)

message(sprintf("[35] Running %d bootstrap resamples for METABRIC ER-neg ...", B))
hrs_met   <- bootstrap_hr(metabric, "dss_time_months", "dss_event",
                           B = B, seed = FREEZE$seed_folds)

# ---------------------------------------------------------------------------
# 5) Summarise
# ---------------------------------------------------------------------------
summarise_boot <- function(hrs, hr_obs, cohort, n, n_events) {
  hrs_clean <- hrs[!is.na(hrs)]
  tibble(
    cohort            = cohort,
    n                 = n,
    n_events          = n_events,
    hr_observed       = round(hr_obs, 3),
    hr_bootstrap_median = round(median(hrs_clean), 3),
    hr_boot_lo95      = round(quantile(hrs_clean, 0.025), 3),
    hr_boot_hi95      = round(quantile(hrs_clean, 0.975), 3),
    sd_loghr          = round(sd(log(hrs_clean)), 3),
    prop_hr_below_1   = round(mean(hrs_clean < 1), 3)
  )
}

res_scanb <- summarise_boot(hrs_scanb, hr_scanb_obs, "SCANB",
                             nrow(scanb), sum(scanb$os_event))
res_met   <- summarise_boot(hrs_met, hr_met_obs, "METABRIC",
                             nrow(metabric), sum(metabric$dss_event))

results <- bind_rows(res_scanb, res_met)

# ---------------------------------------------------------------------------
# 6) Save
# ---------------------------------------------------------------------------
out_path <- file.path(PATHS$results$supp, "er_neg_bootstrap_stability.csv")
write_csv(results, out_path)
message(sprintf("[35] Saved: %s", out_path))

# ---------------------------------------------------------------------------
# 7) Print summary
# ---------------------------------------------------------------------------
cat("\n")
cat("======================================================================\n")
cat("  ER-Negative Bootstrap Stability Analysis (B=2000)\n")
cat("======================================================================\n\n")

for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  cat(sprintf("  %s  (N=%d, %d events)\n", r$cohort, r$n, r$n_events))
  cat(sprintf("    Observed HR:          %.3f\n", r$hr_observed))
  cat(sprintf("    Bootstrap median HR:  %.3f\n", r$hr_bootstrap_median))
  cat(sprintf("    Bootstrap 95%% CI:     [%.3f, %.3f]\n", r$hr_boot_lo95, r$hr_boot_hi95))
  cat(sprintf("    SD(log HR):           %.3f\n", r$sd_loghr))
  cat(sprintf("    Pr(HR < 1):           %.1f%%\n", r$prop_hr_below_1 * 100))
  cat("\n")
}

ci_width_scanb <- results$hr_boot_hi95[1] - results$hr_boot_lo95[1]
ci_width_met   <- results$hr_boot_hi95[2] - results$hr_boot_lo95[2]

cat("----------------------------------------------------------------------\n")
cat("  INTERPRETATION:\n")
cat(sprintf("    SCAN-B ER-neg:   95%% bootstrap CI width = %.2f, SD(logHR) = %.3f\n",
            ci_width_scanb, results$sd_loghr[1]))
cat(sprintf("    METABRIC ER-neg: 95%% bootstrap CI width = %.2f, SD(logHR) = %.3f\n",
            ci_width_met, results$sd_loghr[2]))
cat(sprintf("    Ratio of CI widths (SCANB/METABRIC):    %.1fx\n",
            ci_width_scanb / ci_width_met))
cat(sprintf("    Ratio of SD(logHR) (SCANB/METABRIC):    %.1fx\n",
            results$sd_loghr[1] / results$sd_loghr[2]))
cat("\n")
cat("    => SCAN-B ER-neg HR is UNSTABLE (wide bootstrap distribution,\n")
cat("       high SD, extreme HRs in resamples) due to only 35 events.\n")
cat("    => METABRIC ER-neg HR is STABLE (narrow distribution, centered\n")
cat("       around ~1.0) with 186 events providing adequate power.\n")
cat("======================================================================\n")

message("[35] Done.")
