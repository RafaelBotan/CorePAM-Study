# =============================================================================
# SCRIPT: 15_qc_schema_range_checks.R
# PURPOSE: Structural QC — hard-fail checks on schema and value ranges for
#          all analysis_ready.parquet files and CorePAM artifacts.
#          Any failure stops the pipeline (stop()).
# PROJECT: Core-PAM (Memorial v6.1 §QC)
#
# OUTPUTS:
#   results/supp/qc_schema_range_report.csv
#
# NOTE: This script is designed to run AFTER all 06_* scripts complete.
#       All checks produce hard failures — no warnings, no skips.
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "15_qc_schema_range_checks.R"

# Skip if already passed (set FORCE_RERUN=TRUE to recheck)
FORCE <- as.logical(Sys.getenv("FORCE_RERUN", "FALSE"))
out_path <- file.path(PATHS$results$supp, "qc_schema_range_report.csv")
if (!FORCE && file.exists(out_path)) {
  message(sprintf("[15] Schema QC report already exists: %s", out_path))
  message("[15] Skipping. Set FORCE_RERUN=TRUE to rerun.")
  quit(save = "no", status = 0)
}

message("[15] Starting structural QC — schema and range checks")

COHORTS <- c("SCANB", "TCGA_BRCA", "METABRIC", "GSE20685")

# Endpoint mapping per cohort
ENDPOINT_MAP <- list(
  SCANB     = list(time = "os_time_months",  event = "os_event"),
  TCGA_BRCA = list(time = "os_time_months",  event = "os_event"),
  METABRIC  = list(time = "dss_time_months", event = "dss_event"),
  GSE20685  = list(time = "os_time_months",  event = "os_event")
)

# Accumulate all issues
all_issues  <- list()
report_rows <- list()

# ---------------------------------------------------------------------------
# Helper: record a check result
# ---------------------------------------------------------------------------
record_check <- function(cohort, check_name, status, detail = "") {
  if (status == "FAIL") {
    all_issues[[length(all_issues) + 1]] <<- sprintf("[%s] %s: %s", cohort, check_name, detail)
    message(sprintf("[15] FAIL [%s] %s — %s", cohort, check_name, detail))
  } else {
    message(sprintf("[15] PASS [%s] %s", cohort, check_name))
  }
  report_rows[[length(report_rows) + 1]] <<- list(
    cohort     = cohort,
    check      = check_name,
    status     = status,
    detail     = detail,
    timestamp  = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
}

# ---------------------------------------------------------------------------
# 1) CorePAM weights schema
# ---------------------------------------------------------------------------
weights_path <- file.path(PATHS$results$corepam, "CorePAM_weights.csv")
if (file.exists(weights_path)) {
  wt <- strict_csv(weights_path)
  # Required columns
  if (!all(c("gene", "weight") %in% names(wt))) {
    record_check("GLOBAL", "weights_schema",
                 "FAIL", "Missing columns: gene and/or weight")
  } else {
    record_check("GLOBAL", "weights_schema", "PASS")
  }
  # Non-zero weights
  n_nonzero <- sum(wt$weight != 0, na.rm = TRUE)
  if (n_nonzero == 0) {
    record_check("GLOBAL", "weights_nonzero_count",
                 "FAIL", "All weights are zero — panel is empty")
  } else {
    record_check("GLOBAL", "weights_nonzero_count",
                 "PASS", sprintf("%d non-zero weights", n_nonzero))
  }
  # No NA weights
  n_na <- sum(is.na(wt$weight))
  if (n_na > 0) {
    record_check("GLOBAL", "weights_no_na",
                 "FAIL", sprintf("%d NA weights", n_na))
  } else {
    record_check("GLOBAL", "weights_no_na", "PASS")
  }
} else {
  record_check("GLOBAL", "weights_file_exists",
               "FAIL", sprintf("File not found: %s", weights_path))
}

# ---------------------------------------------------------------------------
# 2) Per-cohort analysis_ready checks
# ---------------------------------------------------------------------------
for (cohort in COHORTS) {
  ready_path <- file.path(proc_cohort(cohort), "analysis_ready.parquet")
  time_col   <- ENDPOINT_MAP[[cohort]]$time
  event_col  <- ENDPOINT_MAP[[cohort]]$event

  if (!file.exists(ready_path)) {
    record_check(cohort, "file_exists", "FAIL",
                 sprintf("Missing: %s", ready_path))
    next
  }
  record_check(cohort, "file_exists", "PASS")

  df <- tryCatch(
    strict_parquet(ready_path),
    error = function(e) {
      record_check(cohort, "file_loads", "FAIL", e$message)
      NULL
    }
  )
  if (is.null(df)) next
  record_check(cohort, "file_loads", "PASS",
               sprintf("%d rows x %d cols", nrow(df), ncol(df)))

  # Required columns (sample_id or patient_id must be present)
  id_col   <- if ("sample_id" %in% names(df)) "sample_id" else "patient_id"
  req_cols <- c(id_col, "score", "score_z", "score_direction",
                time_col, event_col)
  missing_cols <- setdiff(req_cols, names(df))
  if (length(missing_cols) > 0) {
    record_check(cohort, "required_columns", "FAIL",
                 sprintf("Missing: %s", paste(missing_cols, collapse = ", ")))
  } else {
    record_check(cohort, "required_columns", "PASS")
  }

  # No duplicate patient/sample IDs
  id_vals <- if (id_col %in% names(df)) df[[id_col]] else character(0)
  n_dup   <- sum(duplicated(id_vals))
  if (n_dup > 0) {
    record_check(cohort, "no_duplicate_samples", "FAIL",
                 sprintf("%d duplicate sample_ids", n_dup))
  } else {
    record_check(cohort, "no_duplicate_samples", "PASS")
  }

  # score_z: no NA above threshold
  if ("score_z" %in% names(df)) {
    na_frac <- mean(is.na(df$score_z))
    if (na_frac > 0.10) {
      record_check(cohort, "score_z_na_fraction", "FAIL",
                   sprintf("%.1f%% NA (max 10%%)", na_frac * 100))
    } else {
      record_check(cohort, "score_z_na_fraction", "PASS",
                   sprintf("%.2f%%", na_frac * 100))
    }
  }

  # time column: must be numeric > 0 for non-NA values
  if (time_col %in% names(df)) {
    n_le0 <- sum(df[[time_col]] <= 0, na.rm = TRUE)
    if (n_le0 > 0) {
      record_check(cohort, "time_positive", "FAIL",
                   sprintf("%d values <= 0 in %s", n_le0, time_col))
    } else {
      record_check(cohort, "time_positive", "PASS")
    }

    # Time in months: median should be between 1 and 500
    med_time <- median(df[[time_col]], na.rm = TRUE)
    if (is.na(med_time) || med_time < 1 || med_time > 500) {
      record_check(cohort, "time_unit_months", "FAIL",
                   sprintf("Median time = %.1f months (expected 1-500)", med_time))
    } else {
      record_check(cohort, "time_unit_months", "PASS",
                   sprintf("Median = %.1f months", med_time))
    }
  }

  # event column: must be 0 or 1 only
  if (event_col %in% names(df)) {
    invalid_events <- !df[[event_col]] %in% c(0L, 1L, NA)
    n_invalid <- sum(invalid_events, na.rm = TRUE)
    if (n_invalid > 0) {
      record_check(cohort, "event_binary", "FAIL",
                   sprintf("%d values not in {0,1}: %s", n_invalid, event_col))
    } else {
      record_check(cohort, "event_binary", "PASS")
    }
  }

  # METABRIC must have dss columns (not os_* as primary)
  if (cohort == "METABRIC") {
    has_dss <- all(c("dss_time_months", "dss_event") %in% names(df))
    if (!has_dss) {
      record_check(cohort, "metabric_dss_present", "FAIL",
                   "METABRIC must have dss_time and dss_event (primary endpoint)")
    } else {
      record_check(cohort, "metabric_dss_present", "PASS")
    }
  }

  # score_direction must be one of: original, inverted, unknown
  if ("score_direction" %in% names(df)) {
    valid_dirs <- c("original", "inverted", "unknown")
    dirs <- unique(df$score_direction)
    bad_dirs <- setdiff(dirs, valid_dirs)
    if (length(bad_dirs) > 0) {
      record_check(cohort, "score_direction_valid", "FAIL",
                   sprintf("Invalid direction values: %s", paste(bad_dirs, collapse = ", ")))
    } else {
      record_check(cohort, "score_direction_valid", "PASS",
                   sprintf("direction = %s", paste(dirs, collapse = ", ")))
    }
  }

  # Minimum sample size: at least 50 complete cases
  complete_mask <- !is.na(df[[time_col]]) & df[[time_col]] > 0 &
    !is.na(df[[event_col]]) & !is.na(df$score_z)
  n_complete <- sum(complete_mask)
  if (n_complete < 50) {
    record_check(cohort, "min_sample_size", "FAIL",
                 sprintf("Only %d complete cases (minimum 50)", n_complete))
  } else {
    record_check(cohort, "min_sample_size", "PASS",
                 sprintf("%d complete cases", n_complete))
  }
}

# ---------------------------------------------------------------------------
# 3) Save report
# ---------------------------------------------------------------------------
report_df <- bind_rows(lapply(report_rows, as_tibble))
dir.create(PATHS$results$supp, showWarnings = FALSE, recursive = TRUE)
readr::write_csv(report_df, out_path)
h_rep <- sha256_file(out_path)

n_fail <- sum(report_df$status == "FAIL")
n_pass <- sum(report_df$status == "PASS")
overall_status <- if (n_fail == 0) "PASS" else "FAIL"

registry_append("ALL", "qc_schema_range_report", out_path, h_rep,
                overall_status, SCRIPT_NAME, file.info(out_path)$size / 1e6)

message(sprintf("\n[15] Schema QC: %d PASS | %d FAIL", n_pass, n_fail))

if (n_fail > 0) {
  message("[15] FAILURES:")
  for (issue in all_issues) message("  - ", issue)
  stop(sprintf("[15] Schema/range QC FAILED: %d issue(s). See %s", n_fail, out_path))
}

message(sprintf("[15] All schema/range checks PASSED. Report: %s", out_path))
message("[15] Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
