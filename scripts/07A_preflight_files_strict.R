# =============================================================================
# SCRIPT: 07A_preflight_files_strict.R
# PURPOSE: Strict pre-flight check — verifies that all required input files
#          exist and are non-empty before survival analysis scripts (07_*) run.
#          Fails hard (stop()) if any required file is missing.
# PROJETO: Core-PAM (Memorial v6.1 §8)
#
# Run this before any 07_survival_analysis_<COHORT>.R script.
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "07A_preflight_files_strict.R"

# Skip if all outputs already verified (unless force rerun)
FORCE <- as.logical(Sys.getenv("FORCE_RERUN", "FALSE"))
stamp_path <- file.path(PATHS$results$supp, "07A_preflight_passed.txt")
if (!FORCE && file.exists(stamp_path)) {
  message("[07A] Pre-flight stamp found — already passed. Set FORCE_RERUN=TRUE to recheck.")
  quit(save = "no", status = 0)
}

message("[07A] Starting strict pre-flight check for all cohorts...")

COHORTS <- c("SCANB", "TCGA_BRCA", "METABRIC", "GSE20685")

# Required files per cohort
REQUIRED_FILES <- list(
  analysis_ready   = "analysis_ready.parquet",
  expression_preZ  = "expression_genelevel_preZ.parquet",
  clinical_final   = "clinical_FINAL.parquet"
)

# Global model artifacts
GLOBAL_REQUIRED <- list(
  weights  = file.path(PATHS$results$corepam, "CorePAM_weights.csv"),
  model    = file.path(PATHS$results$corepam, "CorePAM_model.rds"),
  card     = file.path(PATHS$results$corepam, "CorePAM_training_card.json")
)

all_ok  <- TRUE
issues  <- character(0)
checked <- character(0)

# ---------------------------------------------------------------------------
# 1) Global artifacts
# ---------------------------------------------------------------------------
message("[07A] Checking global CorePAM artifacts...")
for (nm in names(GLOBAL_REQUIRED)) {
  path <- GLOBAL_REQUIRED[[nm]]
  if (!file.exists(path)) {
    msg <- sprintf("MISSING global artifact [%s]: %s", nm, path)
    issues  <- c(issues, msg)
    all_ok  <- FALSE
    message(sprintf("[07A] FAIL — %s", msg))
  } else {
    sz <- file.info(path)$size
    if (sz == 0) {
      msg <- sprintf("EMPTY global artifact [%s]: %s", nm, path)
      issues  <- c(issues, msg)
      all_ok  <- FALSE
      message(sprintf("[07A] FAIL — %s", msg))
    } else {
      message(sprintf("[07A] OK   — %s (%.1f KB)", nm, sz / 1024))
      checked <- c(checked, path)
    }
  }
}

# ---------------------------------------------------------------------------
# 2) Per-cohort files
# ---------------------------------------------------------------------------
message("[07A] Checking per-cohort files...")
for (cohort in COHORTS) {
  proc_dir <- proc_cohort(cohort)
  for (nm in names(REQUIRED_FILES)) {
    path <- file.path(proc_dir, REQUIRED_FILES[[nm]])
    if (!file.exists(path)) {
      msg <- sprintf("MISSING [%s] %s: %s", cohort, nm, path)
      issues  <- c(issues, msg)
      all_ok  <- FALSE
      message(sprintf("[07A] FAIL — %s", msg))
    } else {
      sz <- file.info(path)$size
      if (sz == 0) {
        msg <- sprintf("EMPTY [%s] %s: %s", cohort, nm, path)
        issues  <- c(issues, msg)
        all_ok  <- FALSE
        message(sprintf("[07A] FAIL — %s", msg))
      } else {
        message(sprintf("[07A] OK   — [%s] %s (%.2f MB)", cohort, nm, sz / 1e6))
        checked <- c(checked, path)
      }
    }
  }
}

# ---------------------------------------------------------------------------
# 3) Validate CorePAM_weights.csv structure
# ---------------------------------------------------------------------------
weights_path <- GLOBAL_REQUIRED$weights
if (file.exists(weights_path) && file.info(weights_path)$size > 0) {
  wt <- strict_csv(weights_path)
  if (!all(c("gene", "weight") %in% names(wt))) {
    msg <- "CorePAM_weights.csv missing required columns: gene, weight"
    issues <- c(issues, msg)
    all_ok <- FALSE
    message(sprintf("[07A] FAIL — %s", msg))
  } else {
    n_nonzero <- sum(wt$weight != 0, na.rm = TRUE)
    message(sprintf("[07A] CorePAM panel: %d non-zero weights", n_nonzero))
    if (n_nonzero == 0) {
      msg <- "CorePAM_weights.csv has 0 non-zero weights — panel is empty!"
      issues <- c(issues, msg)
      all_ok <- FALSE
      message(sprintf("[07A] FAIL — %s", msg))
    }
  }
}

# ---------------------------------------------------------------------------
# 4) Save pre-flight report
# ---------------------------------------------------------------------------
dir.create(PATHS$results$supp, showWarnings = FALSE, recursive = TRUE)

report_df <- tibble(
  timestamp   = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  n_checked   = length(checked),
  n_issues    = length(issues),
  status      = if (all_ok) "PASS" else "FAIL",
  issues_text = if (length(issues) > 0) paste(issues, collapse = " | ") else ""
)

report_path <- file.path(PATHS$results$supp, "07A_preflight_report.csv")
readr::write_csv(report_df, report_path)
h_rep <- sha256_file(report_path)
registry_append("ALL", "preflight_report", report_path, h_rep, report_df$status,
                SCRIPT_NAME, file.info(report_path)$size / 1e6)

# ---------------------------------------------------------------------------
# 5) Hard-fail or stamp
# ---------------------------------------------------------------------------
if (!all_ok) {
  message("\n[07A] ===== PRE-FLIGHT FAILED =====")
  for (issue in issues) message("  - ", issue)
  message("Fix the issues above before running 07_survival_analysis_*.R scripts.")
  stop(sprintf("[07A] Pre-flight failed: %d issue(s). See %s", length(issues), report_path))
}

# Write pass stamp
writeLines(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), stamp_path)
message(sprintf("[07A] ===== PRE-FLIGHT PASSED ===== (%d files checked)", length(checked)))
message(sprintf("[07A] Stamp written: %s", stamp_path))
message("[07A] Ready to run 07_survival_analysis_*.R scripts.")
message("[07A] Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
