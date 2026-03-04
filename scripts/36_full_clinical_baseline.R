# =============================================================================
# SCRIPT: 36_full_clinical_baseline.R
# PURPOSE: Test CorePAM's incremental prognostic value beyond a comprehensive
#          clinical baseline model: age + ER + T-stage + N-status + grade + HER2.
#          For SCAN-B (OS) and METABRIC (DSS) — the two cohorts with all six
#          clinical covariates available.
#          Also repeats the analysis restricted to ER+ patients.
#          Reports Delta C-index (bootstrap), LRT, multivariate HR for
#          CorePAM score, and calibration slope for each model.
# OUTPUT:
#   results/supp/full_clinical_baseline_incremental.csv
# PROJECT: Core-PAM (Memorial CorePAM — Major Revision)
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "36_full_clinical_baseline.R"

suppressPackageStartupMessages({
  library(survival)
  library(arrow)
})

set.seed(FREEZE$seed_folds)
message(sprintf("[%s] ========== FULL CLINICAL BASELINE ANALYSIS ==========",
                SCRIPT_NAME))

# ==========================================================================
# 1) Prepare SCAN-B data: merge analysis_ready with staging from RAW
# ==========================================================================
message(sprintf("[%s] --- Preparing SCAN-B ---", SCRIPT_NAME))

scanb_raw <- strict_rds(file.path(PATHS$raw, "GSE96058/GSE96058_clinical_raw.rds"))
scanb_ar  <- as.data.frame(strict_parquet(
  file.path(PATHS$processed, "SCANB/analysis_ready.parquet")
))

# Derive T-stage from tumor size (mm) per AJCC criteria
old_warn <- getOption("warn"); options(warn = 0)
scanb_size <- as.numeric(scanb_raw[["tumor size:ch1"]])
options(warn = old_warn)

scanb_t <- ifelse(scanb_size <= 20, "T1",
           ifelse(scanb_size <= 50, "T2",
           ifelse(scanb_size >  50, "T3", NA_character_)))

# Derive N-status from lymph node status
scanb_node_raw <- scanb_raw[["lymph node status:ch1"]]
scanb_n <- ifelse(scanb_node_raw == "NodeNegative", "N0",
           ifelse(scanb_node_raw == "NodePositive", "N+", NA_character_))

scanb_staging <- data.frame(
  sample_id = normalize_id(scanb_raw[["geo_accession"]]),
  t_stage   = scanb_t,
  n_status  = scanb_n,
  stringsAsFactors = FALSE
)
message(sprintf("[%s] SCAN-B staging: T available %d/%d (%.1f%%), N available %d/%d (%.1f%%)",
                SCRIPT_NAME,
                sum(!is.na(scanb_t)), nrow(scanb_raw),
                100 * mean(!is.na(scanb_t)),
                sum(!is.na(scanb_n)), nrow(scanb_raw),
                100 * mean(!is.na(scanb_n))))

# Merge staging into analysis_ready
scanb_df <- merge(scanb_ar, scanb_staging, by = "sample_id", all.x = TRUE)
scanb_df$cohort_name <- "SCAN-B"
scanb_df$time_col    <- scanb_df$os_time_months
scanb_df$event_col   <- scanb_df$os_event
message(sprintf("[%s] SCAN-B merged: N=%d", SCRIPT_NAME, nrow(scanb_df)))

# ==========================================================================
# 2) Prepare METABRIC data: merge analysis_ready with staging from RAW
# ==========================================================================
message(sprintf("[%s] --- Preparing METABRIC ---", SCRIPT_NAME))

metabric_dir <- file.path(PATHS$raw, "METABRIC/brca_metabric")

.read_cbio <- function(fname) {
  fpath <- file.path(metabric_dir, fname)
  if (!file.exists(fpath)) stop("File not found: ", fpath)
  lines <- readLines(fpath, warn = FALSE)
  skip_n <- sum(startsWith(lines, "#"))
  old_warn <- getOption("warn"); options(warn = 0)
  out <- readr::read_tsv(fpath, skip = skip_n, show_col_types = FALSE,
                          name_repair = "unique")
  options(warn = old_warn)
  out
}

met_pat <- .read_cbio("data_clinical_patient.txt")
met_sam <- .read_cbio("data_clinical_sample.txt")
met_full <- merge(as.data.frame(met_pat), as.data.frame(met_sam),
                  by = "PATIENT_ID", all = FALSE)

# Derive T-stage from TUMOR_SIZE (mm)
old_warn <- getOption("warn"); options(warn = 0)
met_size <- as.numeric(met_full[["TUMOR_SIZE"]])
options(warn = old_warn)

met_t <- ifelse(met_size <= 20, "T1",
         ifelse(met_size <= 50, "T2",
         ifelse(met_size >  50, "T3", NA_character_)))

# Derive N-status from LYMPH_NODES_EXAMINED_POSITIVE (0 = N0, >0 = N+)
old_warn <- getOption("warn"); options(warn = 0)
met_nodes <- as.numeric(met_full[["LYMPH_NODES_EXAMINED_POSITIVE"]])
options(warn = old_warn)

met_n <- ifelse(!is.na(met_nodes) & met_nodes == 0, "N0",
         ifelse(!is.na(met_nodes) & met_nodes > 0,  "N+", NA_character_))

met_staging <- data.frame(
  sample_id = normalize_id(met_full[["SAMPLE_ID"]]),
  t_stage   = met_t,
  n_status  = met_n,
  stringsAsFactors = FALSE
)
message(sprintf("[%s] METABRIC staging: T available %d/%d (%.1f%%), N available %d/%d (%.1f%%)",
                SCRIPT_NAME,
                sum(!is.na(met_t)), nrow(met_full),
                100 * mean(!is.na(met_t)),
                sum(!is.na(met_n)), nrow(met_full),
                100 * mean(!is.na(met_n))))

met_ar <- as.data.frame(strict_parquet(
  file.path(PATHS$processed, "METABRIC/analysis_ready.parquet")
))
met_df <- merge(met_ar, met_staging, by = "sample_id", all.x = TRUE)
met_df$cohort_name <- "METABRIC"
met_df$time_col    <- met_df$dss_time_months
met_df$event_col   <- met_df$dss_event
message(sprintf("[%s] METABRIC merged: N=%d", SCRIPT_NAME, nrow(met_df)))

# ==========================================================================
# 3) Analysis function: incremental value of CorePAM over full clinical model
# ==========================================================================
run_incremental <- function(df, cohort, endpoint, population, time_var,
                            event_var, horizon = 60) {
  message(sprintf("[%s] --- %s | %s | %s ---",
                  SCRIPT_NAME, cohort, endpoint, population))

  # Prepare clinical variables
  df$age_           <- df$age
  df$t_stage_       <- factor(df$t_stage, levels = c("T1", "T2", "T3"))
  df$n_status_      <- factor(df$n_status, levels = c("N0", "N+"))
  df$grade_         <- factor(df$grade, levels = c(1, 2, 3), ordered = TRUE)
  df$score_z_       <- df$score_z
  df$surv_time_     <- df[[time_var]]
  df$surv_event_    <- df[[event_var]]

  # ER status: convert to factor
  if (is.numeric(df$er_status)) {
    df$er_status_ <- factor(ifelse(df$er_status == 1, "Positive", "Negative"),
                            levels = c("Negative", "Positive"))
  } else {
    df$er_status_ <- factor(df$er_status, levels = c("Negative", "Positive"))
  }

  # HER2 status: convert to factor
  df$her2_status_ <- factor(df$her2_status, levels = c(0, 1),
                            labels = c("Negative", "Positive"))

  # Filter to complete cases
  cov_vars <- c("surv_time_", "surv_event_", "age_", "er_status_",
                "t_stage_", "n_status_", "grade_", "her2_status_", "score_z_")
  df_cc <- df[complete.cases(df[, cov_vars]), ]
  df_cc <- df_cc[df_cc$surv_time_ > 0, ]

  n_total <- nrow(df_cc)
  n_events <- sum(df_cc$surv_event_)
  n_cov_complete <- n_total
  message(sprintf("[%s] Complete cases: N=%d, Events=%d", SCRIPT_NAME,
                  n_total, n_events))

  if (n_total < 50 || n_events < 20) {
    message(sprintf("[%s] Insufficient data for %s %s %s (N=%d, ev=%d). Skipping.",
                    SCRIPT_NAME, cohort, endpoint, population, n_total, n_events))
    return(NULL)
  }

  # ------------------------------------------------------------------
  # 3a) Fit Cox models
  # ------------------------------------------------------------------
  fm_clinical <- Surv(surv_time_, surv_event_) ~ age_ + er_status_ +
    t_stage_ + n_status_ + grade_ + her2_status_
  fm_full     <- Surv(surv_time_, surv_event_) ~ age_ + er_status_ +
    t_stage_ + n_status_ + grade_ + her2_status_ + score_z_

  old_warn <- getOption("warn"); options(warn = 0)
  cox_clinical <- tryCatch(coxph(fm_clinical, data = df_cc), error = function(e) NULL)
  cox_full     <- tryCatch(coxph(fm_full,     data = df_cc), error = function(e) NULL)
  options(warn = old_warn)

  if (is.null(cox_clinical) || is.null(cox_full)) {
    message(sprintf("[%s] Cox fitting failed for %s %s %s",
                    SCRIPT_NAME, cohort, endpoint, population))
    return(NULL)
  }

  # ------------------------------------------------------------------
  # 3b) C-indices
  # ------------------------------------------------------------------
  old_warn <- getOption("warn"); options(warn = 0)
  c_clinical <- concordance(cox_clinical)$concordance
  c_full     <- concordance(cox_full)$concordance
  options(warn = old_warn)

  delta_c_obs <- c_full - c_clinical

  # ------------------------------------------------------------------
  # 3c) Bootstrap Delta C-index (1000 resamples)
  # ------------------------------------------------------------------
  n_boot <- as.integer(FREEZE$bootstrap_n)
  boot_seed <- as.integer(FREEZE$seed_folds)
  set.seed(boot_seed)
  n <- nrow(df_cc)
  delta_boot <- numeric(n_boot)

  old_warn <- getOption("warn"); options(warn = 0)
  for (b in seq_len(n_boot)) {
    idx <- sample(n, n, replace = TRUE)
    df_b <- df_cc[idx, ]
    cox_b_clin <- tryCatch(coxph(fm_clinical, data = df_b), error = function(e) NULL)
    cox_b_full <- tryCatch(coxph(fm_full,     data = df_b), error = function(e) NULL)
    if (!is.null(cox_b_clin) && !is.null(cox_b_full)) {
      c_b_clin <- tryCatch(concordance(cox_b_clin)$concordance, error = function(e) NA_real_)
      c_b_full <- tryCatch(concordance(cox_b_full)$concordance, error = function(e) NA_real_)
      delta_boot[b] <- c_b_full - c_b_clin
    } else {
      delta_boot[b] <- NA_real_
    }
  }
  options(warn = old_warn)

  delta_lo95 <- quantile(delta_boot, 0.025, na.rm = TRUE)
  delta_hi95 <- quantile(delta_boot, 0.975, na.rm = TRUE)

  message(sprintf("[%s] C_clinical=%.4f | C_clinical+CorePAM=%.4f | DeltaC=%.4f (%.4f, %.4f)",
                  SCRIPT_NAME, c_clinical, c_full, delta_c_obs,
                  delta_lo95, delta_hi95))

  # ------------------------------------------------------------------
  # 3d) Likelihood Ratio Test (nested models)
  # ------------------------------------------------------------------
  old_warn <- getOption("warn"); options(warn = 0)
  lrt <- tryCatch(anova(cox_clinical, cox_full, test = "Chisq"),
                  error = function(e) NULL)
  options(warn = old_warn)

  lrt_chisq <- lrt_df <- lrt_p <- NA_real_
  if (!is.null(lrt)) {
    # anova returns a data.frame; row 2 is the comparison
    lrt_chisq <- lrt[["Chisq"]][2]
    lrt_df    <- lrt[["Df"]][2]
    lrt_p     <- lrt[["Pr(>|Chi|)"]][2]
    message(sprintf("[%s] LRT: chisq=%.3f, df=%d, p=%g",
                    SCRIPT_NAME, lrt_chisq, lrt_df, lrt_p))
  }

  # ------------------------------------------------------------------
  # 3e) Multivariate HR for score_z in full model
  # ------------------------------------------------------------------
  sm_full <- summary(cox_full)
  coef_idx <- grep("score_z_", rownames(sm_full$coefficients))
  hr_adj <- hr_adj_lo <- hr_adj_hi <- p_adj <- NA_real_
  if (length(coef_idx) > 0) {
    hr_adj    <- sm_full$conf.int[coef_idx, "exp(coef)"]
    hr_adj_lo <- sm_full$conf.int[coef_idx, "lower .95"]
    hr_adj_hi <- sm_full$conf.int[coef_idx, "upper .95"]
    p_adj     <- sm_full$coefficients[coef_idx, "Pr(>|z|)"]
    message(sprintf("[%s] Adjusted HR=%.3f (%.3f-%.3f), p=%g",
                    SCRIPT_NAME, hr_adj, hr_adj_lo, hr_adj_hi, p_adj))
  }

  # ------------------------------------------------------------------
  # 3f) Calibration slope at horizon
  # ------------------------------------------------------------------
  cal_slope_clinical <- cal_slope_full <- NA_real_

  old_warn <- getOption("warn"); options(warn = 0)

  # Clinical-only calibration slope
  cal_clin <- tryCatch({
    lp_clin <- predict(cox_clinical, type = "lp")
    # Binary outcome at horizon: event observed before horizon
    obs_event_h <- as.integer(df_cc$surv_time_ <= horizon & df_cc$surv_event_ == 1)
    # Exclude those censored before horizon without event
    include <- df_cc$surv_time_ > horizon | df_cc$surv_event_ == 1
    if (sum(include) >= 30) {
      glm_clin <- glm(obs_event_h[include] ~ lp_clin[include], family = binomial)
      coef(glm_clin)[2]  # slope (ideal = 1)
    } else NA_real_
  }, error = function(e) NA_real_)
  cal_slope_clinical <- cal_clin

  # Clinical + CorePAM calibration slope
  cal_full <- tryCatch({
    lp_full <- predict(cox_full, type = "lp")
    obs_event_h <- as.integer(df_cc$surv_time_ <= horizon & df_cc$surv_event_ == 1)
    include <- df_cc$surv_time_ > horizon | df_cc$surv_event_ == 1
    if (sum(include) >= 30) {
      glm_full <- glm(obs_event_h[include] ~ lp_full[include], family = binomial)
      coef(glm_full)[2]
    } else NA_real_
  }, error = function(e) NA_real_)
  cal_slope_full <- cal_full

  options(warn = old_warn)

  message(sprintf("[%s] Calibration slope (clinical): %.3f | (clinical+CorePAM): %.3f",
                  SCRIPT_NAME,
                  ifelse(is.na(cal_slope_clinical), NA, cal_slope_clinical),
                  ifelse(is.na(cal_slope_full), NA, cal_slope_full)))

  # ------------------------------------------------------------------
  # Return row
  # ------------------------------------------------------------------
  tibble(
    cohort                = cohort,
    endpoint              = endpoint,
    population            = population,
    n                     = n_total,
    n_events              = n_events,
    n_covariates_complete = n_cov_complete,
    c_clinical            = round(c_clinical, 4),
    c_clinical_corepam    = round(c_full, 4),
    delta_c               = round(delta_c_obs, 4),
    delta_c_lo95          = round(as.numeric(delta_lo95), 4),
    delta_c_hi95          = round(as.numeric(delta_hi95), 4),
    hr_corepam_adjusted   = round(hr_adj, 4),
    hr_adj_lo95           = round(hr_adj_lo, 4),
    hr_adj_hi95           = round(hr_adj_hi, 4),
    p_corepam_adjusted    = signif(p_adj, 4),
    lrt_chisq             = round(lrt_chisq, 3),
    lrt_df                = lrt_df,
    lrt_p                 = signif(lrt_p, 4),
    cal_slope_clinical    = round(cal_slope_clinical, 4),
    cal_slope_full        = round(cal_slope_full, 4),
    horizon_months        = horizon
  )
}

# ==========================================================================
# 4) Run analyses for both cohorts: All patients + ER+ restricted
# ==========================================================================
results_all <- list()

# --- SCAN-B: All ---
results_all[["SCANB_All"]] <- run_incremental(
  df = scanb_df, cohort = "SCAN-B", endpoint = "OS",
  population = "All", time_var = "os_time_months",
  event_var = "os_event", horizon = 60
)

# --- SCAN-B: ER+ ---
scanb_erp <- scanb_df[!is.na(scanb_df$er_status) &
                         tolower(scanb_df$er_status) %in% c("positive", "1"), ]
message(sprintf("[%s] SCAN-B ER+ subset: N=%d", SCRIPT_NAME, nrow(scanb_erp)))
results_all[["SCANB_ERpos"]] <- run_incremental(
  df = scanb_erp, cohort = "SCAN-B", endpoint = "OS",
  population = "ER+", time_var = "os_time_months",
  event_var = "os_event", horizon = 60
)

# --- METABRIC: All ---
results_all[["METABRIC_All"]] <- run_incremental(
  df = met_df, cohort = "METABRIC", endpoint = "DSS",
  population = "All", time_var = "dss_time_months",
  event_var = "dss_event", horizon = 60
)

# --- METABRIC: ER+ ---
met_erp <- met_df[!is.na(met_df$er_status) &
                    tolower(met_df$er_status) %in% c("positive", "1"), ]
message(sprintf("[%s] METABRIC ER+ subset: N=%d", SCRIPT_NAME, nrow(met_erp)))
results_all[["METABRIC_ERpos"]] <- run_incremental(
  df = met_erp, cohort = "METABRIC", endpoint = "DSS",
  population = "ER+", time_var = "dss_time_months",
  event_var = "dss_event", horizon = 60
)

# ==========================================================================
# 5) Combine and save
# ==========================================================================
results_df <- bind_rows(Filter(Negate(is.null), results_all))

if (nrow(results_df) == 0) {
  stop(sprintf("[%s] No results generated. Check data availability.", SCRIPT_NAME))
}

out_dir <- PATHS$results$supp
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(out_dir, "full_clinical_baseline_incremental.csv")
write_csv(results_df, out_path)

h <- sha256_file(out_path)
registry_append("SCANB_METABRIC", "full_clinical_baseline_incremental",
                out_path, h, "ok", SCRIPT_NAME,
                file.info(out_path)$size / 1e6)

message(sprintf("\n[%s] Results saved: %s", SCRIPT_NAME, out_path))
message(sprintf("[%s] %d analysis rows generated", SCRIPT_NAME, nrow(results_df)))

# Print summary table
message(sprintf("\n[%s] ===== SUMMARY =====", SCRIPT_NAME))
for (i in seq_len(nrow(results_df))) {
  r <- results_df[i, ]
  message(sprintf(
    "  %s | %s | %s | N=%d (ev=%d) | C_clin=%.4f | C_full=%.4f | DC=%.4f [%.4f, %.4f] | HR=%.3f (%.3f-%.3f) p=%g | LRT p=%g",
    r$cohort, r$endpoint, r$population,
    r$n, r$n_events,
    r$c_clinical, r$c_clinical_corepam,
    r$delta_c, r$delta_c_lo95, r$delta_c_hi95,
    r$hr_corepam_adjusted, r$hr_adj_lo95, r$hr_adj_hi95,
    r$p_corepam_adjusted,
    r$lrt_p
  ))
}

message(sprintf("\n[%s] ========== DONE ==========", SCRIPT_NAME))
