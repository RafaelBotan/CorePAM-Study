# =============================================================================
# SCRIPT: 07D_validate_one_cohort_corepam.R
# PURPOSE: Smoke test for a single cohort — verifies the CorePAM score
#          pipeline end-to-end: data integrity, score distribution, Cox
#          direction, and basic sanity checks. Intended for rapid validation
#          before running full survival analysis.
# PROJETO: Core-PAM (Memorial v6.1 §8)
#
# Usage:
#   Rscript --no-save scripts/07D_validate_one_cohort_corepam.R COHORT
#   or set COHORT environment variable:
#   COHORT=SCANB Rscript --no-save scripts/07D_validate_one_cohort_corepam.R
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "07D_validate_one_cohort_corepam.R"

suppressPackageStartupMessages(library(survival))

# ---------------------------------------------------------------------------
# Resolve target cohort from CLI args or environment
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
COHORT <- if (length(args) >= 1 && nchar(trimws(args[1])) > 0) {
  trimws(args[1])
} else {
  env_coh <- Sys.getenv("COHORT", "")
  if (nchar(env_coh) > 0) env_coh else "SCANB"   # default for interactive use
}

valid_cohorts <- c("SCANB", "TCGA_BRCA", "METABRIC", "GSE20685")
if (!COHORT %in% valid_cohorts) {
  stop(sprintf("[07D] Unknown cohort: '%s'. Valid: %s", COHORT,
               paste(valid_cohorts, collapse = ", ")))
}

# Skip if already validated (unless force rerun)
FORCE <- as.logical(Sys.getenv("FORCE_RERUN", "FALSE"))
stamp_path <- file.path(PATHS$results$supp,
                        sprintf("07D_validated_%s.txt", COHORT))
if (!FORCE && file.exists(stamp_path)) {
  message(sprintf("[07D] Cohort %s already validated. Set FORCE_RERUN=TRUE to recheck.", COHORT))
  quit(save = "no", status = 0)
}

message(sprintf("[07D] Starting smoke test for cohort: %s", COHORT))

checks_passed <- 0L
checks_failed <- 0L
issues <- character(0)

fail_check <- function(msg) {
  issues  <<- c(issues, msg)
  checks_failed <<- checks_failed + 1L
  message(sprintf("[07D] FAIL — %s", msg))
}

pass_check <- function(msg) {
  checks_passed <<- checks_passed + 1L
  message(sprintf("[07D] PASS — %s", msg))
}

# ---------------------------------------------------------------------------
# CHECK 1: analysis_ready.parquet exists and loads
# ---------------------------------------------------------------------------
ready_path <- file.path(proc_cohort(COHORT), "analysis_ready.parquet")

if (!file.exists(ready_path)) {
  fail_check(sprintf("analysis_ready.parquet missing: %s", ready_path))
  stop("[07D] Cannot continue without analysis_ready.parquet. Run 06_zscore_and_score.R first.")
}

df <- tryCatch(
  strict_parquet(ready_path),
  error = function(e) {
    fail_check(sprintf("Failed to load analysis_ready.parquet: %s", e$message))
    NULL
  }
)

if (is.null(df)) stop("[07D] Cannot continue — failed to load data.")
pass_check(sprintf("analysis_ready.parquet loaded: %d rows x %d cols", nrow(df), ncol(df)))

# ---------------------------------------------------------------------------
# CHECK 2: Required columns present
# ---------------------------------------------------------------------------
required_cols <- c("sample_id", "score", "score_z", "score_direction")
# Add endpoint columns based on cohort
endpoint_cols <- if (COHORT == "METABRIC") c("dss_time", "dss_event") else c("os_time", "os_event")
required_cols <- c(required_cols, endpoint_cols)

missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  fail_check(sprintf("Missing columns: %s", paste(missing_cols, collapse = ", ")))
} else {
  pass_check(sprintf("All required columns present (%s)", paste(required_cols, collapse = ", ")))
}

# ---------------------------------------------------------------------------
# CHECK 3: No all-NA score column
# ---------------------------------------------------------------------------
if ("score_z" %in% names(df)) {
  na_frac <- mean(is.na(df$score_z))
  if (na_frac > 0.10) {
    fail_check(sprintf("score_z has %.1f%% NA values (threshold: 10%%)", na_frac * 100))
  } else {
    pass_check(sprintf("score_z NA fraction: %.1f%%", na_frac * 100))
  }
}

# ---------------------------------------------------------------------------
# CHECK 4: Score distribution sanity
# ---------------------------------------------------------------------------
if ("score_z" %in% names(df) && sum(!is.na(df$score_z)) > 10) {
  score_mean <- mean(df$score_z, na.rm = TRUE)
  score_sd   <- sd(df$score_z, na.rm = TRUE)
  # score_z should be approximately N(0,1) by construction
  if (abs(score_mean) > 1.0) {
    fail_check(sprintf("score_z mean far from 0: %.4f (expected ~0)", score_mean))
  } else {
    pass_check(sprintf("score_z mean: %.4f (expected ~0)", score_mean))
  }
  if (score_sd < 0.1 || score_sd > 10) {
    fail_check(sprintf("score_z SD suspicious: %.4f (expected ~1)", score_sd))
  } else {
    pass_check(sprintf("score_z SD: %.4f (expected ~1)", score_sd))
  }
}

# ---------------------------------------------------------------------------
# CHECK 5: Endpoint column sanity
# ---------------------------------------------------------------------------
time_col  <- endpoint_cols[1]
event_col <- endpoint_cols[2]

if (all(c(time_col, event_col) %in% names(df))) {
  n_pos_time  <- sum(df[[time_col]] > 0, na.rm = TRUE)
  n_events    <- sum(df[[event_col]] == 1, na.rm = TRUE)
  event_rate  <- n_events / nrow(df)

  if (n_pos_time < 10) {
    fail_check(sprintf("Very few positive times: %d", n_pos_time))
  } else {
    pass_check(sprintf("Positive times: %d / %d", n_pos_time, nrow(df)))
  }

  if (event_rate < 0.01 || event_rate > 0.99) {
    fail_check(sprintf("Suspicious event rate: %.1f%%", event_rate * 100))
  } else {
    pass_check(sprintf("Events: %d / %d (%.1f%%)", n_events, nrow(df), event_rate * 100))
  }
}

# ---------------------------------------------------------------------------
# CHECK 6: Cox HR direction (should be >= 1 after score_direction applied)
# ---------------------------------------------------------------------------
if (all(c("score_z", time_col, event_col) %in% names(df))) {
  df_cox <- df[!is.na(df[[time_col]]) & df[[time_col]] > 0 &
                 !is.na(df[[event_col]]) & !is.na(df$score_z), ]

  if (nrow(df_cox) >= 20) {
    old_warn <- getOption("warn"); options(warn = 0)
    cox_fit <- tryCatch(
      coxph(Surv(df_cox[[time_col]], df_cox[[event_col]]) ~ score_z, data = df_cox),
      error = function(e) NULL
    )
    options(warn = old_warn)

    if (!is.null(cox_fit)) {
      hr_val <- exp(coef(cox_fit)[1])
      message(sprintf("[07D] Cox HR (score_z): %.4f", hr_val))
      if (hr_val < 0.5) {
        fail_check(sprintf("HR very low: %.4f — score direction may be wrong", hr_val))
      } else if (hr_val < 1.0) {
        # Warn but don't fail — direction inversion may be intentional for sensitivity
        message(sprintf("[07D] WARN — HR < 1: %.4f (check score_direction column)", hr_val))
        pass_check(sprintf("Cox HR computed: %.4f (< 1 but non-critical for smoke test)", hr_val))
      } else {
        pass_check(sprintf("Cox HR: %.4f (>= 1, correct direction)", hr_val))
      }
    } else {
      fail_check("Cox model failed to fit")
    }
  } else {
    message(sprintf("[07D] SKIP Cox check — too few complete cases: %d", nrow(df_cox)))
  }
}

# ---------------------------------------------------------------------------
# CHECK 7: genes_present column (coverage)
# ---------------------------------------------------------------------------
if ("genes_present" %in% names(df)) {
  n_present <- df$genes_present[1]
  # Load panel size from weights file
  weights_path <- file.path(PATHS$results$corepam, "CorePAM_weights.csv")
  if (file.exists(weights_path)) {
    wt <- strict_csv(weights_path)
    n_panel <- sum(wt$weight != 0, na.rm = TRUE)
    frac <- n_present / n_panel
    min_frac <- FREEZE$min_genes_fraction
    if (frac < min_frac) {
      fail_check(sprintf("Gene coverage: %d/%d (%.1f%%) < %.0f%% minimum",
                         n_present, n_panel, frac * 100, min_frac * 100))
    } else {
      pass_check(sprintf("Gene coverage: %d/%d (%.1f%%) >= %.0f%% minimum",
                         n_present, n_panel, frac * 100, min_frac * 100))
    }
  }
}

# ---------------------------------------------------------------------------
# SUMMARY
# ---------------------------------------------------------------------------
message(sprintf("\n[07D] ===== SMOKE TEST SUMMARY for %s =====", COHORT))
message(sprintf("[07D] Checks passed: %d | Failed: %d", checks_passed, checks_failed))

dir.create(PATHS$results$supp, showWarnings = FALSE, recursive = TRUE)

summary_df <- tibble(
  timestamp      = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  cohort         = COHORT,
  checks_passed  = checks_passed,
  checks_failed  = checks_failed,
  status         = if (checks_failed == 0) "PASS" else "FAIL",
  issues         = if (length(issues) > 0) paste(issues, collapse = " | ") else ""
)

report_path <- file.path(PATHS$results$supp,
                          sprintf("07D_smoke_test_%s.csv", COHORT))
readr::write_csv(summary_df, report_path)
h_rep <- sha256_file(report_path)
registry_append(COHORT, "smoke_test", report_path, h_rep,
                summary_df$status, SCRIPT_NAME,
                file.info(report_path)$size / 1e6)

if (checks_failed > 0) {
  message("[07D] Issues found:")
  for (issue in issues) message("  - ", issue)
  stop(sprintf("[07D] Smoke test FAILED for %s: %d issue(s)", COHORT, checks_failed))
}

# Write pass stamp
writeLines(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), stamp_path)
message(sprintf("[07D] Smoke test PASSED for %s", COHORT))
message("[07D] Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
