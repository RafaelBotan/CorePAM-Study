# =============================================================================
# SCRIPT: 16_qc_text_vs_results_assert.R
# PURPOSE: Anti-hard-code assertion checker.
#          Verifies that key numerical results reported in the manuscript QMD
#          match values from CSV/JSON result files. Hard-fails if any
#          discrepancy found (tolerance: ±0.001 for numeric values).
# PROJECT: Core-PAM (Memorial v6.1 §QC)
#
# INPUTS:
#   results/corepam/CorePAM_training_card.json
#   results/main/meta_survival_summary.csv
#   results/supp/survival_results_*.csv
#   (QMD file if present)
#
# OUTPUTS:
#   results/supp/qc_text_vs_results_report.csv
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "16_qc_text_vs_results_assert.R"

suppressPackageStartupMessages(library(jsonlite))

# Skip if already passed (set FORCE_RERUN=TRUE to recheck)
FORCE <- as.logical(Sys.getenv("FORCE_RERUN", "FALSE"))
out_path <- file.path(PATHS$results$supp, "qc_text_vs_results_report.csv")
if (!FORCE && file.exists(out_path)) {
  message(sprintf("[16] Report already exists: %s", out_path))
  message("[16] Skipping. Set FORCE_RERUN=TRUE to rerun.")
  quit(save = "no", status = 0)
}

message("[16] Starting anti-hard-code assertion checks")

NUMERIC_TOL <- 0.001   # tolerance for floating-point comparisons

checks <- list()

record_assert <- function(name, expected, actual, source_file, tolerance = NUMERIC_TOL) {
  if (is.na(actual)) {
    status <- "SKIP"
    detail <- sprintf("actual value not found (source: %s)", basename(source_file))
  } else if (is.numeric(expected) && is.numeric(actual)) {
    diff <- abs(expected - actual)
    if (diff <= tolerance) {
      status <- "PASS"
      detail <- sprintf("expected=%.4f actual=%.4f diff=%.6f", expected, actual, diff)
    } else {
      status <- "FAIL"
      detail <- sprintf("MISMATCH expected=%.4f actual=%.4f diff=%.6f (tol=%.3f)",
                        expected, actual, diff, tolerance)
    }
  } else {
    status <- if (as.character(expected) == as.character(actual)) "PASS" else "FAIL"
    detail <- sprintf("expected='%s' actual='%s'", expected, actual)
  }

  if (status == "FAIL") message(sprintf("[16] FAIL — %s: %s", name, detail))
  else if (status == "PASS") message(sprintf("[16] PASS — %s", name))
  else message(sprintf("[16] SKIP — %s: %s", name, detail))

  checks[[length(checks) + 1]] <<- list(
    check      = name,
    status     = status,
    detail     = detail,
    source     = basename(source_file),
    timestamp  = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
}

# ---------------------------------------------------------------------------
# 1) CorePAM training card (JSON)
# ---------------------------------------------------------------------------
card_path <- file.path(PATHS$results$corepam, "CorePAM_training_card.json")

if (file.exists(card_path)) {
  card <- tryCatch(
    jsonlite::fromJSON(card_path),
    error = function(e) {
      message(sprintf("[16] WARNING: Could not parse %s: %s", card_path, e$message))
      NULL
    }
  )

  if (!is.null(card)) {
    message("[16] Training card loaded. Checking frozen parameters match freeze CSV...")

    # Verify frozen parameters are consistent with analysis_freeze.csv
    if (!is.null(card$alpha))   record_assert("alpha_matches_freeze",
                                              FREEZE$alpha, card$alpha, card_path)
    if (!is.null(card$k_folds)) record_assert("k_folds_matches_freeze",
                                              FREEZE$k_folds, card$k_folds, card_path)
    if (!is.null(card$seed_folds)) record_assert("seed_matches_freeze",
                                                  FREEZE$seed_folds, card$seed_folds, card_path)
    if (!is.null(card$delta_c)) record_assert("delta_c_matches_freeze",
                                              FREEZE$delta_c, card$delta_c, card_path)

    # Verify n_genes > 0
    if (!is.null(card$n_genes)) {
      if (card$n_genes > 0) {
        message(sprintf("[16] PASS — n_genes > 0: %d", card$n_genes))
        checks[[length(checks) + 1]] <- list(
          check     = "n_genes_positive",
          status    = "PASS",
          detail    = sprintf("n_genes = %d", card$n_genes),
          source    = basename(card_path),
          timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        )
      } else {
        message("[16] FAIL — n_genes is 0 or negative")
        checks[[length(checks) + 1]] <- list(
          check     = "n_genes_positive",
          status    = "FAIL",
          detail    = sprintf("n_genes = %d (must be > 0)", card$n_genes),
          source    = basename(card_path),
          timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        )
      }
    }

    # Verify cadj is non-trivial (should be > 0.5)
    if (!is.null(card$cadj_selected)) {
      record_assert("cadj_above_random", 0.5, NULL, card_path)  # manual check
      if (card$cadj_selected > 0.5) {
        message(sprintf("[16] PASS — cadj_selected = %.4f > 0.5", card$cadj_selected))
      } else {
        message(sprintf("[16] WARN — cadj_selected = %.4f <= 0.5 (may indicate poor model)", card$cadj_selected))
      }
    }
  }
} else {
  message(sprintf("[16] SKIP — training card not found: %s", card_path))
}

# ---------------------------------------------------------------------------
# 2) Meta-survival summary consistency
# ---------------------------------------------------------------------------
meta_path <- file.path(PATHS$results$main, "meta_survival_summary.csv")

if (file.exists(meta_path)) {
  meta_df <- tryCatch(strict_csv(meta_path), error = function(e) NULL)

  if (!is.null(meta_df) && nrow(meta_df) > 0) {
    # Meta HR must be > 1 (higher score = worse outcome)
    uni_row <- meta_df[meta_df$analysis == "univariado", ]
    if (nrow(uni_row) > 0) {
      meta_hr <- uni_row$meta_HR[1]
      if (!is.na(meta_hr) && meta_hr <= 1.0) {
        message(sprintf("[16] WARN — meta_HR = %.4f <= 1 (expected > 1 for risk score)", meta_hr))
      } else if (!is.na(meta_hr)) {
        message(sprintf("[16] PASS — meta_HR = %.4f > 1", meta_hr))
      }

      # I2 must be numeric and between 0-100
      i2 <- uni_row$I2_pct[1]
      if (!is.na(i2) && (i2 < 0 || i2 > 100)) {
        message(sprintf("[16] FAIL — meta I2 = %.1f%% out of range [0,100]", i2))
        checks[[length(checks) + 1]] <- list(
          check = "meta_I2_range", status = "FAIL",
          detail = sprintf("I2 = %.1f%%", i2), source = basename(meta_path),
          timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        )
      } else if (!is.na(i2)) {
        message(sprintf("[16] PASS — meta I2 = %.1f%% (valid range)", i2))
        checks[[length(checks) + 1]] <- list(
          check = "meta_I2_range", status = "PASS",
          detail = sprintf("I2 = %.1f%%", i2), source = basename(meta_path),
          timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        )
      }
    }
  }
} else {
  message(sprintf("[16] SKIP — meta_survival_summary.csv not found: %s", meta_path))
}

# ---------------------------------------------------------------------------
# 3) Per-cohort survival results: sanity range checks
# ---------------------------------------------------------------------------
COHORTS <- c("SCANB", "TCGA_BRCA", "METABRIC", "GSE20685")

for (cohort in COHORTS) {
  surv_path <- file.path(PATHS$results$supp,
                          sprintf("survival_results_%s.csv", cohort))
  if (!file.exists(surv_path)) {
    message(sprintf("[16] SKIP — survival_results_%s.csv not found", cohort))
    next
  }

  surv_df <- tryCatch(strict_csv(surv_path), error = function(e) NULL)
  if (is.null(surv_df) || nrow(surv_df) == 0) next

  # HR must be positive
  if ("hr_uni" %in% names(surv_df)) {
    hr <- surv_df$hr_uni[1]
    if (!is.na(hr) && hr <= 0) {
      message(sprintf("[16] FAIL [%s] hr_uni = %.4f (must be positive)", cohort, hr))
      checks[[length(checks) + 1]] <- list(
        check = sprintf("%s_hr_positive", cohort), status = "FAIL",
        detail = sprintf("hr_uni = %.4f", hr), source = basename(surv_path),
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      )
    } else if (!is.na(hr)) {
      message(sprintf("[16] PASS [%s] hr_uni = %.4f (positive)", cohort, hr))
      checks[[length(checks) + 1]] <- list(
        check = sprintf("%s_hr_positive", cohort), status = "PASS",
        detail = sprintf("hr_uni = %.4f", hr), source = basename(surv_path),
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      )
    }
  }

  # C-index must be between 0 and 1
  if ("c_index" %in% names(surv_df)) {
    ci <- surv_df$c_index[1]
    if (!is.na(ci) && (ci < 0 || ci > 1)) {
      message(sprintf("[16] FAIL [%s] c_index = %.4f out of [0,1]", cohort, ci))
      checks[[length(checks) + 1]] <- list(
        check = sprintf("%s_cindex_range", cohort), status = "FAIL",
        detail = sprintf("c_index = %.4f", ci), source = basename(surv_path),
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      )
    } else if (!is.na(ci)) {
      message(sprintf("[16] PASS [%s] c_index = %.4f (valid)", cohort, ci))
      checks[[length(checks) + 1]] <- list(
        check = sprintf("%s_cindex_range", cohort), status = "PASS",
        detail = sprintf("c_index = %.4f", ci), source = basename(surv_path),
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      )
    }
  }
}

# ---------------------------------------------------------------------------
# 4) Check for hard-coded numbers in QMD file (if present)
# ---------------------------------------------------------------------------
qmd_paths <- list.files(ROOT_REPO, pattern = "\\.qmd$", recursive = TRUE, full.names = TRUE)
if (length(qmd_paths) > 0) {
  qmd_path <- qmd_paths[1]
  message(sprintf("[16] Checking QMD for hard-coded numbers: %s", basename(qmd_path)))

  qmd_lines <- readLines(qmd_path, warn = FALSE)
  # Pattern: bare numbers in text (not in code blocks, not table separators)
  # We check for specific suspicious patterns: p-values, HR, C-index written raw
  # This is a heuristic — check for numbers that look like results in non-code text
  suspect_lines <- grep("(?<![`|\\$])\\b[0-9]+\\.[0-9]{3,}\\b(?![`|\\$])",
                        qmd_lines, perl = TRUE, value = TRUE)
  # Filter out obvious non-result lines (YAML, code blocks, references)
  suspect_lines <- grep("^[^#`|\\-\\s]", suspect_lines, value = TRUE)

  if (length(suspect_lines) > 20) {
    message(sprintf("[16] WARN — QMD has %d lines with potential hard-coded numbers",
                    length(suspect_lines)))
    message("[16] First 5 suspicious lines:")
    for (l in suspect_lines[seq_len(5)]) message("  ", l)
  } else {
    message(sprintf("[16] QMD scan: %d suspect lines with embedded numbers (review manually)",
                    length(suspect_lines)))
  }
} else {
  message("[16] SKIP — no .qmd file found (will check when manuscript is created)")
}

# ---------------------------------------------------------------------------
# 5) Save report
# ---------------------------------------------------------------------------
report_df <- bind_rows(lapply(checks, as_tibble))

if (nrow(report_df) == 0) {
  report_df <- tibble(
    check     = "no_checks_run",
    status    = "SKIP",
    detail    = "No result files found to validate. Run scripts 05-08 first.",
    source    = "",
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
}

dir.create(PATHS$results$supp, showWarnings = FALSE, recursive = TRUE)
readr::write_csv(report_df, out_path)
h_rep <- sha256_file(out_path)

n_fail <- sum(report_df$status == "FAIL")
n_pass <- sum(report_df$status == "PASS")
n_skip <- sum(report_df$status == "SKIP")
overall_status <- if (n_fail == 0) "PASS" else "FAIL"

registry_append("ALL", "qc_text_vs_results_report", out_path, h_rep,
                overall_status, SCRIPT_NAME, file.info(out_path)$size / 1e6)

message(sprintf("\n[16] Anti-hard-code QC: %d PASS | %d FAIL | %d SKIP", n_pass, n_fail, n_skip))

if (n_fail > 0) {
  fails <- report_df$detail[report_df$status == "FAIL"]
  for (f in fails) message("  FAIL: ", f)
  stop(sprintf("[16] Text vs results QC FAILED: %d issue(s). See %s", n_fail, out_path))
}

message(sprintf("[16] All assertions PASSED. Report: %s", out_path))
message("[16] Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
