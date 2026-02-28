options(warn = 0)
library(nanoparquet)

COHORTS <- c("SCANB", "TCGA_BRCA", "METABRIC", "GSE20685")
BASE    <- "Y:/Phd-Genomic-claude/01_Base_Pura_CorePAM/PROCESSED"

for (cohort in COHORTS) {
  cat(sprintf("\n=== %s ===\n", cohort))

  # Expression parquet
  expr_path <- file.path(BASE, cohort, "expression_genelevel_preZ.parquet")
  if (file.exists(expr_path)) {
    e <- read_parquet(expr_path)
    cat(sprintf("  Expression dims: %d rows x %d cols\n", nrow(e), ncol(e)))
    cat(sprintf("  First 4 col names: %s\n", paste(head(names(e), 4), collapse=", ")))
    cat(sprintf("  First col values (first 2): %s\n", paste(head(e[[1]], 2), collapse=", ")))
  } else {
    cat("  Expression: NOT FOUND\n")
  }

  # Clinical parquet — try both names
  clin_final <- file.path(BASE, cohort, "clinical_FINAL.parquet")
  clin_harm  <- file.path(BASE, cohort, "clinical_harmonized.parquet")
  clin_path  <- if (file.exists(clin_final)) clin_final
                else if (file.exists(clin_harm)) clin_harm
                else NA

  if (!is.na(clin_path)) {
    c <- read_parquet(clin_path)
    cat(sprintf("  Clinical file: %s\n", basename(clin_path)))
    cat(sprintf("  Clinical cols: %s\n", paste(names(c), collapse=", ")))
    time_cols  <- grep("time|os|dss|surv|months", names(c), ignore.case=TRUE, value=TRUE)
    event_cols <- grep("event|status|dead|vital", names(c), ignore.case=TRUE, value=TRUE)
    cat(sprintf("  Time-like cols: %s\n", paste(time_cols, collapse=", ")))
    cat(sprintf("  Event-like cols: %s\n", paste(event_cols, collapse=", ")))
  } else {
    cat("  Clinical: NOT FOUND (tried clinical_FINAL.parquet and clinical_harmonized.parquet)\n")
  }
}
cat("\nDone.\n")
