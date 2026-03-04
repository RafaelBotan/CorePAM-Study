# =============================================================================
# SCRIPT: 32_calibration_formal.R
# PURPOSE: Formal calibration assessment via logistic recalibration at fixed
#          time horizon. Reports calibration intercept (alpha, ideal=0) and
#          slope (beta, ideal=1). Addresses reviewer request for quantitative
#          calibration beyond visual plots.
# OUTPUT:
#   results/supp/calibration_slope_intercept.csv
# PROJECT: Core-PAM (Memorial CorePAM — Major Revision)
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "32_calibration_formal.R"

suppressPackageStartupMessages({
  library(survival)
  library(arrow)
})

set.seed(FREEZE$seed_folds)
message(sprintf("[%s] ========== FORMAL CALIBRATION ==========", SCRIPT_NAME))

# --------------------------------------------------------------------------
# 1) Define cohorts + endpoints + horizons
# --------------------------------------------------------------------------
COHORTS <- list(
  SCANB     = list(time = "os_time_months",  event = "os_event",  horizon = 60),
  TCGA_BRCA = list(time = "os_time_months",  event = "os_event",  horizon = 24),
  METABRIC  = list(time = "dss_time_months", event = "dss_event", horizon = 60),
  GSE20685  = list(time = "os_time_months",  event = "os_event",  horizon = 60)
)

# --------------------------------------------------------------------------
# 2) Loop over cohorts
# --------------------------------------------------------------------------
results <- list()

for (coh_name in names(COHORTS)) {
  message(sprintf("[%s] Processing %s...", SCRIPT_NAME, coh_name))
  ep <- COHORTS[[coh_name]]

  ar_path <- file.path(proc_cohort(coh_name), "analysis_ready.parquet")
  if (!file.exists(ar_path)) {
    message(sprintf("[%s] %s: analysis_ready.parquet not found, skipping", SCRIPT_NAME, coh_name))
    next
  }

  df <- as.data.frame(strict_parquet(ar_path))
  time_col  <- ep$time
  event_col <- ep$event
  horizon   <- ep$horizon

  if (!all(c(time_col, event_col, "score_z") %in% names(df))) {
    message(sprintf("[%s] %s: required columns missing, skipping", SCRIPT_NAME, coh_name))
    next
  }

  df <- df[!is.na(df[[time_col]]) & df[[time_col]] > 0 &
             !is.na(df[[event_col]]) & !is.na(df$score_z), ]
  n_analysis <- nrow(df)

  if (n_analysis < 50) {
    message(sprintf("[%s] %s: n=%d too small for calibration, skipping", SCRIPT_NAME, coh_name, n_analysis))
    next
  }

  # Check max follow-up is sufficient
  if (max(df[[time_col]], na.rm = TRUE) < horizon) {
    message(sprintf("[%s] %s: max follow-up (%.1f) < horizon (%d), skipping",
                    SCRIPT_NAME, coh_name, max(df[[time_col]], na.rm = TRUE), horizon))
    next
  }

  message(sprintf("[%s] %s: N=%d, Events=%d, Horizon=%dm",
                  SCRIPT_NAME, coh_name, n_analysis, sum(df[[event_col]]), horizon))

  # --- Fit Cox model to get predicted risk at horizon ---
  old_warn <- getOption("warn"); options(warn = 0)
  cox_fit <- tryCatch(
    coxph(Surv(df[[time_col]], df[[event_col]]) ~ score_z, data = df, x = TRUE),
    error = function(e) NULL
  )
  options(warn = old_warn)

  if (is.null(cox_fit)) {
    message(sprintf("[%s] %s: Cox fit failed, skipping", SCRIPT_NAME, coh_name))
    next
  }

  # Predicted survival probability at horizon for each individual
  old_warn <- getOption("warn"); options(warn = 0)
  sf <- tryCatch(survfit(cox_fit, newdata = df), error = function(e) NULL)
  options(warn = old_warn)

  if (is.null(sf)) {
    message(sprintf("[%s] %s: survfit prediction failed, skipping", SCRIPT_NAME, coh_name))
    next
  }

  # S(t=horizon) for each sample
  t_idx <- max(which(sf$time <= horizon))
  pred_surv <- sf$surv[t_idx, ]
  pred_risk <- 1 - pred_surv  # predicted event probability

  # --- Observed binary outcome at horizon ---
  # Y_H = 1 if event before horizon, 0 otherwise (censor at horizon)
  y_h <- ifelse(df[[time_col]] <= horizon & df[[event_col]] == 1, 1, 0)
  # Exclude those censored before horizon (cannot determine outcome)
  usable <- df[[time_col]] > horizon | (df[[time_col]] <= horizon & df[[event_col]] == 1)
  n_usable <- sum(usable)

  if (n_usable < 30) {
    message(sprintf("[%s] %s: n_usable=%d too small for logistic recalibration, skipping",
                    SCRIPT_NAME, coh_name, n_usable))
    next
  }

  y_sub <- y_h[usable]
  pred_sub <- pred_risk[usable]

  # --- Logistic recalibration: logit(Y_H) = alpha + beta * logit(pred_risk) ---
  # Transform predicted risk to log-odds (clamped to avoid Inf)
  pred_clamped <- pmin(pmax(pred_sub, 0.001), 0.999)
  lp_pred <- log(pred_clamped / (1 - pred_clamped))

  old_warn <- getOption("warn"); options(warn = 0)
  logit_fit <- tryCatch(
    glm(y_sub ~ lp_pred, family = binomial(link = "logit")),
    error = function(e) NULL
  )
  options(warn = old_warn)

  alpha <- beta <- alpha_se <- beta_se <- NA_real_
  if (!is.null(logit_fit)) {
    alpha    <- coef(logit_fit)[1]       # intercept: ideal = 0
    beta     <- coef(logit_fit)[2]       # slope: ideal = 1
    alpha_se <- summary(logit_fit)$coefficients[1, "Std. Error"]
    beta_se  <- summary(logit_fit)$coefficients[2, "Std. Error"]

    message(sprintf("[%s] %s: alpha=%.4f (SE=%.4f), beta=%.4f (SE=%.4f)",
                    SCRIPT_NAME, coh_name, alpha, alpha_se, beta, beta_se))
  }

  # --- Calibration-in-the-large (CITL): mean(predicted) vs mean(observed) ---
  citl_pred <- mean(pred_sub, na.rm = TRUE)
  citl_obs  <- mean(y_sub, na.rm = TRUE)

  # --- E/O ratio (expected/observed events) ---
  eo_ratio <- sum(pred_sub, na.rm = TRUE) / sum(y_sub)

  # --- Brier score ---
  brier <- mean((pred_sub - y_sub)^2, na.rm = TRUE)

  results[[coh_name]] <- data.frame(
    cohort       = coh_name,
    horizon_m    = horizon,
    n            = n_analysis,
    n_usable     = n_usable,
    n_events_h   = sum(y_sub),
    alpha        = round(alpha, 4),
    alpha_se     = round(alpha_se, 4),
    beta         = round(beta, 4),
    beta_se      = round(beta_se, 4),
    citl_pred    = round(citl_pred, 4),
    citl_obs     = round(citl_obs, 4),
    eo_ratio     = round(eo_ratio, 4),
    brier_score  = round(brier, 4),
    stringsAsFactors = FALSE
  )
}

# --------------------------------------------------------------------------
# 3) Save results
# --------------------------------------------------------------------------
if (length(results) == 0) {
  stop(sprintf("[%s] No calibration results generated.", SCRIPT_NAME))
}

results_df <- do.call(rbind, results)
out_path   <- file.path(PATHS$results$supp, "calibration_slope_intercept.csv")
write.csv(results_df, out_path, row.names = FALSE)
h <- sha256_file(out_path)
registry_append("ALL", "calibration_slope_intercept", out_path, h, "ok", SCRIPT_NAME,
                file.info(out_path)$size / 1e6)

message(sprintf("[%s] Results saved: %s", SCRIPT_NAME, out_path))
message(sprintf("[%s] ========== DONE ==========", SCRIPT_NAME))
