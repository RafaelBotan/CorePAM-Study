# =============================================================================
# SCRIPT: 00_validate_all.R
# PURPOSE: Master validation runner.
#          1. Syntax-checks all scripts/*.R (parse()).
#          2. Checks that expected output files exist (pipeline progress).
#          3. Reports PASS/FAIL per script in a summary table.
#          Does NOT run the scripts (data may not be present).
#          To run the full pipeline, use 18_make_submission_bundle.R runbook.
#
# Usage:
#   Rscript scripts/00_validate_all.R
#   FORCE_RERUN=TRUE Rscript scripts/00_validate_all.R
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "00_validate_all.R"

FORCE <- as.logical(Sys.getenv("FORCE_RERUN", "FALSE"))

message("[validate] Starting full pipeline validation check")
message("[validate] Syntax check + output file presence check")
message(strrep("=", 60))

# ---------------------------------------------------------------------------
# 1) SYNTAX CHECK: parse all scripts
# ---------------------------------------------------------------------------
script_dir   <- file.path(ROOT_REPO, "scripts")
script_files <- sort(list.files(script_dir, pattern = "\\.R$", full.names = TRUE))
# Exclude self (avoid circular check triggering recursion confusion)
script_files <- script_files[basename(script_files) != "00_validate_all.R"]

syntax_results <- lapply(script_files, function(f) {
  result <- tryCatch({
    parse(file = f)
    list(file = basename(f), status = "PASS", detail = "syntax OK")
  }, error = function(e) {
    list(file = basename(f), status = "FAIL",
         detail = sprintf("syntax error: %s", conditionMessage(e)))
  })
  if (result$status == "FAIL") {
    message(sprintf("[validate] SYNTAX FAIL: %s — %s", result$file, result$detail))
  } else {
    message(sprintf("[validate] SYNTAX PASS: %s", result$file))
  }
  result
})

syntax_df <- bind_rows(lapply(syntax_results, as_tibble))

# ---------------------------------------------------------------------------
# 2) OUTPUT FILE PRESENCE CHECK
# ---------------------------------------------------------------------------
COHORTS <- c("SCANB", "GSE96058", "TCGA_BRCA", "METABRIC", "GSE20685")

# Define expected outputs for each pipeline stage
expected_outputs <- list(
  # Stage 01: downloads (check at least one raw dir exists)
  list(stage = "01_download",
       label = "RAW dirs exist",
       paths = lapply(COHORTS, function(c) raw_cohort(c))),

  # Stage 02: clinical harmonization
  list(stage = "02_clinical",
       label = "clinical_FINAL.parquet",
       paths = lapply(COHORTS, function(c)
         file.path(proc_cohort(c), "clinical_FINAL.parquet"))),

  # Stage 03: expression preprocessing
  list(stage = "03_expression",
       label = "expression_genelevel_preZ.parquet",
       paths = lapply(COHORTS, function(c)
         file.path(proc_cohort(c), "expression_genelevel_preZ.parquet"))),

  # Stage 04: gene audit
  list(stage = "04_audit",
       label = "gene_audit_by_cohort.csv",
       paths = list(file.path(PATHS$results$supp, "gene_audit_by_cohort.csv"))),

  # Stage 05: CorePAM derivation
  list(stage = "05_corepam",
       label = "CorePAM_weights.csv",
       paths = list(file.path(PATHS$results$corepam, "CorePAM_weights.csv"))),
  list(stage = "05_corepam",
       label = "CorePAM_training_card.json",
       paths = list(file.path(PATHS$results$corepam, "CorePAM_training_card.json"))),

  # Stage 06: scoring
  list(stage = "06_score",
       label = "analysis_ready.parquet",
       paths = lapply(COHORTS, function(c)
         file.path(proc_cohort(c), "analysis_ready.parquet"))),

  # Stage 07: survival analysis
  list(stage = "07_survival",
       label = "survival_results_*.csv",
       paths = lapply(COHORTS, function(c)
         file.path(PATHS$results$supp, sprintf("survival_results_%s.csv", c)))),

  # Stage 08: meta-analysis
  list(stage = "08_meta",
       label = "meta_survival_summary.csv",
       paths = list(file.path(PATHS$results$main, "meta_survival_summary.csv"))),

  # Stage 11: incremental value
  list(stage = "11_incremental",
       label = "incremental_value_results.csv",
       paths = list(file.path(PATHS$results$supp, "incremental_value_results.csv"))),

  # Stage 13: QC correlations
  list(stage = "13_qc_corr",
       label = "qc_offdiag_correlations.csv",
       paths = list(file.path(PATHS$results$supp, "qc_offdiag_correlations.csv"))),

  # Stage 14: METABRIC PCA
  list(stage = "14_pca",
       label = "metabric_pca_forensics.csv",
       paths = list(file.path(PATHS$results$supp, "metabric_pca_forensics.csv"))),

  # Stage 15: schema QC
  list(stage = "15_schema",
       label = "qc_schema_range_report.csv",
       paths = list(file.path(PATHS$results$supp, "qc_schema_range_report.csv"))),

  # Stage 16: text vs results
  list(stage = "16_text_vs",
       label = "qc_text_vs_results_report.csv",
       paths = list(file.path(PATHS$results$supp, "qc_text_vs_results_report.csv")))
)

message("\n[validate] Checking expected output files...")
output_results <- list()

for (item in expected_outputs) {
  paths_to_check <- unlist(item$paths)
  n_total   <- length(paths_to_check)
  n_present <- sum(file.exists(paths_to_check))
  status    <- if (n_present == n_total) "PASS" else
               if (n_present > 0) "PARTIAL" else "MISSING"

  detail <- sprintf("%d/%d files present", n_present, n_total)
  if (status != "PASS") {
    missing_files <- paths_to_check[!file.exists(paths_to_check)]
    detail <- paste0(detail, sprintf("; missing: %s",
                                      paste(basename(missing_files), collapse = ", ")))
  }
  message(sprintf("[validate] %-20s %-40s %s",
                  item$stage, item$label, status))

  output_results[[length(output_results) + 1]] <- list(
    check  = paste(item$stage, item$label, sep = ": "),
    status = status,
    detail = detail
  )
}

output_df <- bind_rows(lapply(output_results, as_tibble))

# ---------------------------------------------------------------------------
# 3) Summary report
# ---------------------------------------------------------------------------
message(strrep("=", 60))
message("[validate] SUMMARY")
message(strrep("=", 60))

n_syntax_fail <- sum(syntax_df$status == "FAIL")
n_syntax_pass <- sum(syntax_df$status == "PASS")
n_output_pass <- sum(output_df$status == "PASS")
n_output_part <- sum(output_df$status == "PARTIAL")
n_output_miss <- sum(output_df$status == "MISSING")

message(sprintf("[validate] Syntax: %d PASS | %d FAIL", n_syntax_pass, n_syntax_fail))
message(sprintf("[validate] Outputs: %d PASS | %d PARTIAL | %d MISSING",
                n_output_pass, n_output_part, n_output_miss))

if (n_syntax_fail > 0) {
  message("\n[validate] SYNTAX FAILURES:")
  fails <- syntax_df[syntax_df$status == "FAIL", ]
  for (i in seq_len(nrow(fails))) {
    message(sprintf("  %s: %s", fails$file[i], fails$detail[i]))
  }
}

# Pipeline stage completion
total_stages  <- nrow(output_df)
stages_done   <- n_output_pass
pct_done      <- round(100 * stages_done / total_stages)
message(sprintf("\n[validate] Pipeline completion: %d/%d stages (%d%%)",
                stages_done, total_stages, pct_done))

# ---------------------------------------------------------------------------
# 4) Save full report
# ---------------------------------------------------------------------------
dir.create(PATHS$results$supp, showWarnings = FALSE, recursive = TRUE)

report_path <- file.path(PATHS$results$supp, "validation_report.csv")
full_report <- bind_rows(
  mutate(syntax_df, check_type = "syntax"),
  mutate(output_df, check_type = "output_presence", file = NA_character_)
)
readr::write_csv(full_report, report_path)

message(sprintf("\n[validate] Full report saved: %s", report_path))

if (n_syntax_fail > 0) {
  stop(sprintf("[validate] %d syntax error(s) found. Fix before running pipeline.",
               n_syntax_fail))
}

message("[validate] All scripts pass syntax check.")
message("[validate] Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
