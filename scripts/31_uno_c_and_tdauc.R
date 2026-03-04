# =============================================================================
# SCRIPT: 31_uno_c_and_tdauc.R
# PURPOSE: Compute Uno's C-statistic (IPCW) and time-dependent AUC (tdAUC)
#          as complementary discrimination metrics that handle heterogeneous
#          censoring. Addresses reviewer request for IPCW-based metrics.
# PACKAGES: survC1, timeROC
# OUTPUT:
#   results/supp/unos_c_and_tdauc.csv
# PROJECT: Core-PAM (Memorial CorePAM — Major Revision)
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "31_uno_c_and_tdauc.R"

suppressPackageStartupMessages({
  library(survival)
  library(arrow)
  library(survC1)
  library(timeROC)
})

set.seed(FREEZE$seed_folds)
message(sprintf("[%s] ========== UNO'S C + tdAUC ==========", SCRIPT_NAME))

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

  df <- strict_parquet(ar_path)
  time_col  <- ep$time
  event_col <- ep$event
  horizon   <- ep$horizon

  if (!all(c(time_col, event_col, "score_z") %in% names(df))) {
    message(sprintf("[%s] %s: required columns missing, skipping", SCRIPT_NAME, coh_name))
    next
  }

  df <- df[!is.na(df[[time_col]]) & df[[time_col]] > 0 &
             !is.na(df[[event_col]]) & !is.na(df$score_z), ]
  df <- as.data.frame(df)
  n_analysis <- nrow(df)

  if (n_analysis < 30) {
    message(sprintf("[%s] %s: n=%d too small, skipping", SCRIPT_NAME, coh_name, n_analysis))
    next
  }

  message(sprintf("[%s] %s: N=%d, Events=%d, Horizon=%dm",
                  SCRIPT_NAME, coh_name, n_analysis, sum(df[[event_col]]), horizon))

  # --- Uno's C-statistic via survC1::Est.Cval ---
  uno_c <- uno_c_lo <- uno_c_hi <- uno_c_se <- NA_real_

  old_warn <- getOption("warn"); options(warn = 0)
  uno_res <- tryCatch({
    Est.Cval(
      mydata  = data.frame(
        time  = df[[time_col]],
        event = df[[event_col]],
        score = df$score_z
      ),
      tau     = horizon,
      nofit   = TRUE
    )
  }, error = function(e) {
    message(sprintf("[%s] %s: Uno's C failed: %s", SCRIPT_NAME, coh_name, e$message))
    NULL
  })
  options(warn = old_warn)

  if (!is.null(uno_res)) {
    uno_c <- uno_res$Dhat
    # Est.Cval with nofit=TRUE may not return SE; extract safely
    if (!is.null(uno_res$se) && is.numeric(uno_res$se) && length(uno_res$se) == 1) {
      uno_c_se <- uno_res$se
    } else {
      # Use Jackknife/perturbation SE if available
      uno_c_se <- NA_real_
    }
    if (!is.na(uno_c_se)) {
      uno_c_lo <- uno_c - 1.96 * uno_c_se
      uno_c_hi <- uno_c + 1.96 * uno_c_se
    }
    message(sprintf("[%s] %s: Uno's C = %.4f (SE=%s)",
                    SCRIPT_NAME, coh_name, uno_c,
                    ifelse(is.na(uno_c_se), "NA", sprintf("%.4f", uno_c_se))))
  }

  # --- Time-dependent AUC via timeROC ---
  tdauc <- tdauc_lo <- tdauc_hi <- tdauc_se <- NA_real_

  old_warn <- getOption("warn"); options(warn = 0)
  troc_res <- tryCatch({
    timeROC(
      T         = df[[time_col]],
      delta     = df[[event_col]],
      marker    = df$score_z,
      cause     = 1,
      times     = horizon,
      iid       = TRUE,
      weighting = "marginal"
    )
  }, error = function(e) {
    message(sprintf("[%s] %s: timeROC failed: %s", SCRIPT_NAME, coh_name, e$message))
    NULL
  })
  options(warn = old_warn)

  if (!is.null(troc_res)) {
    idx_h <- which(troc_res$times == horizon)
    if (length(idx_h) > 0) {
      tdauc_raw <- troc_res$AUC[idx_h]
      # timeROC returns AUC on 0-1 scale (not percentage)
      tdauc <- tdauc_raw

      # Confidence interval via confint method
      old_warn <- getOption("warn"); options(warn = 0)
      ci_res <- tryCatch({
        confint(troc_res, level = 0.95)
      }, error = function(e) NULL)
      options(warn = old_warn)

      if (!is.null(ci_res) && !is.null(ci_res$CI_AUC)) {
        ci_mat <- ci_res$CI_AUC
        # confint.timeROC returns CI as PERCENTAGES (0-100 scale); convert to 0-1
        if (is.matrix(ci_mat) && nrow(ci_mat) >= 1) {
          ci_row <- min(idx_h, nrow(ci_mat))
          tdauc_lo <- ci_mat[ci_row, 1] / 100
          tdauc_hi <- ci_mat[ci_row, 2] / 100
        } else if (is.numeric(ci_mat) && length(ci_mat) >= 2) {
          tdauc_lo <- ci_mat[1] / 100
          tdauc_hi <- ci_mat[2] / 100
        }
      }

      # Fallback: compute SE from iid and build CI manually
      if (is.na(tdauc_lo) && !is.null(troc_res$inference) && !is.null(troc_res$inference$vect_iid_comp_time)) {
        iid_vec <- troc_res$inference$vect_iid_comp_time[, idx_h]
        tdauc_se <- sd(iid_vec)
        tdauc_lo <- tdauc - 1.96 * tdauc_se
        tdauc_hi <- tdauc + 1.96 * tdauc_se
      }

      message(sprintf("[%s] %s: tdAUC(%dm) = %.4f (%.4f-%.4f)",
                      SCRIPT_NAME, coh_name, horizon, tdauc,
                      ifelse(is.na(tdauc_lo), NA, tdauc_lo),
                      ifelse(is.na(tdauc_hi), NA, tdauc_hi)))
    }
  }

  # --- Harrell's C (for comparison) ---
  old_warn <- getOption("warn"); options(warn = 0)
  harrell_c <- tryCatch({
    cx <- concordance(Surv(df[[time_col]], df[[event_col]]) ~ df$score_z)
    max(cx$concordance, 1 - cx$concordance)
  }, error = function(e) NA_real_)
  options(warn = old_warn)

  # --- Store results ---
  results[[coh_name]] <- data.frame(
    cohort      = coh_name,
    n           = n_analysis,
    n_events    = sum(df[[event_col]]),
    horizon_m   = horizon,
    # Uno's C
    uno_c       = round(uno_c, 4),
    uno_c_lo95  = round(uno_c_lo, 4),
    uno_c_hi95  = round(uno_c_hi, 4),
    uno_c_se    = if (is.na(uno_c_se)) NA_real_ else round(uno_c_se, 4),
    # tdAUC
    tdauc       = round(tdauc, 4),
    tdauc_lo95  = round(ifelse(is.na(tdauc_lo), NA, tdauc_lo), 4),
    tdauc_hi95  = round(ifelse(is.na(tdauc_hi), NA, tdauc_hi), 4),
    # Harrell's C (reference)
    harrell_c   = round(harrell_c, 4),
    stringsAsFactors = FALSE
  )
}

# --------------------------------------------------------------------------
# 3) Save results
# --------------------------------------------------------------------------
if (length(results) == 0) {
  stop(sprintf("[%s] No results generated.", SCRIPT_NAME))
}

results_df <- do.call(rbind, results)
out_path   <- file.path(PATHS$results$supp, "unos_c_and_tdauc.csv")
write.csv(results_df, out_path, row.names = FALSE)
h <- sha256_file(out_path)
registry_append("ALL", "unos_c_and_tdauc", out_path, h, "ok", SCRIPT_NAME,
                file.info(out_path)$size / 1e6)

message(sprintf("[%s] Results saved: %s", SCRIPT_NAME, out_path))
message(sprintf("[%s] ========== DONE ==========", SCRIPT_NAME))
