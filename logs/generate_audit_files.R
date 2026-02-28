options(warn = 0)
library(nanoparquet)

# =============================================================================
# generate_audit_files.R
# Generates:
#   1) results/supp/pam50_gene_list_audit.csv
#   2) results/supp/train_inclusion_log.csv
# =============================================================================

source("scripts/00_setup.R")

PAM50_CANONICAL_50 <- c(
  "ACTR3B","ANLN","BAG1","BCL2","BIRC5","BLVRA","CCNB1","CCNE1","CDC20","CDC6",
  "CDH3","CENPF","CEP55","CXXC5","EGFR","ERBB2","ESR1","EXO1","FGFR4","FOXA1",
  "FOXC1","GPR160","GRB7","KIF2C","KRT14","KRT17","KRT5","MAPT","MDM2",
  "MELK","MIA","MKI67","MLPH","MMP11","MYBL2","MYC","NAT1","NDC80","NUF2",
  "ORC6","PGR","PHGDH","PTTG1","RRM2","SFRP1","SLC39A6","TMEM45B","TYMS","UBE2C","UBE2T"
)
stopifnot(length(PAM50_CANONICAL_50) == 50)

COHORTS <- c("SCANB","TCGA_BRCA","METABRIC","GSE20685")

# ── 1) PAM50 gene audit across cohorts ───────────────────────────────────────
cat("=== PAM50 Gene List Audit ===\n")
cat("Canonical list: Parker et al. 2009 (J Clin Oncol), 50 genes\n\n")

audit_rows <- list()
for (cohort in COHORTS) {
  path_expr <- file.path(proc_cohort(cohort), "expression_genelevel_preZ.parquet")
  if (!file.exists(path_expr)) {
    cat(cohort, ": expression file not found — skipping\n")
    next
  }
  expr <- nanoparquet::read_parquet(path_expr)
  genes_in_matrix <- expr$gene
  for (g in PAM50_CANONICAL_50) {
    present <- g %in% genes_in_matrix
    note <- if (!present) paste0("absent in ", cohort) else ""
    audit_rows[[length(audit_rows) + 1]] <- data.frame(
      cohort         = cohort,
      gene           = g,
      in_canonical_50 = TRUE,
      in_matrix      = present,
      note           = note,
      stringsAsFactors = FALSE
    )
  }
  n_present <- sum(PAM50_CANONICAL_50 %in% genes_in_matrix)
  cat(sprintf("  %-12s: %d/50 PAM50 genes present\n", cohort, n_present))
}

audit_df <- do.call(rbind, audit_rows)

# Also flag KRT8 explicitly as "not in canonical 50"
krt8_row <- data.frame(
  cohort          = COHORTS,
  gene            = "KRT8",
  in_canonical_50 = FALSE,
  in_matrix       = vapply(COHORTS, function(cohort) {
    path_expr <- file.path(proc_cohort(cohort), "expression_genelevel_preZ.parquet")
    if (!file.exists(path_expr)) return(NA)
    expr <- nanoparquet::read_parquet(path_expr)
    "KRT8" %in% expr$gene
  }, logical(1)),
  note = "NOT in canonical PAM50 (Parker 2009) — excluded from analysis",
  stringsAsFactors = FALSE
)
audit_df <- rbind(audit_df, krt8_row)

dir.create("results/supp", showWarnings = FALSE, recursive = TRUE)
write.csv(audit_df, "results/supp/pam50_gene_list_audit.csv", row.names = FALSE)
cat("\nAudit saved: results/supp/pam50_gene_list_audit.csv\n")
cat("Total rows:", nrow(audit_df), "\n")

# Wide summary: gene × cohort
audit_wide <- reshape(
  audit_df[audit_df$in_canonical_50 == TRUE, c("cohort","gene","in_matrix")],
  idvar = "gene", timevar = "cohort", direction = "wide"
)
names(audit_wide) <- sub("in_matrix\\.", "", names(audit_wide))
audit_wide$n_cohorts_present <- rowSums(audit_wide[, COHORTS[COHORTS %in% names(audit_wide)]], na.rm=TRUE)
write.csv(audit_wide, "results/supp/pam50_gene_list_audit_wide.csv", row.names = FALSE)
cat("Wide audit saved: results/supp/pam50_gene_list_audit_wide.csv\n\n")

# ── 2) Train inclusion log (SCAN-B) ──────────────────────────────────────────
cat("=== Train Inclusion Log (SCAN-B) ===\n")

# Raw clinical from parquet
path_clin_raw <- file.path(proc_cohort("SCANB"), "clinical_FINAL.parquet")
if (!file.exists(path_clin_raw)) {
  cat("clinical_FINAL.parquet not found — cannot build inclusion log\n")
} else {
  clin <- nanoparquet::read_parquet(path_clin_raw)

  # Expression parquet for matching
  path_expr_scanb <- file.path(proc_cohort("SCANB"), "expression_genelevel_preZ.parquet")
  expr_scanb <- nanoparquet::read_parquet(path_expr_scanb)
  # Get sample IDs from expression matrix (column names except 'gene')
  expr_ids <- setdiff(names(expr_scanb), "gene")

  # Step-by-step inclusion
  n_raw           <- nrow(clin)
  # Detect time/event column names (lowercase in SCANB parquet)
  time_col  <- if ("os_time_months" %in% names(clin)) "os_time_months" else "OS_Time"
  event_col <- if ("os_event"       %in% names(clin)) "os_event"       else "OS_Event"
  cat(sprintf("  Detected: time=%s | event=%s\n", time_col, event_col))

  n_has_os_time   <- sum(!is.na(clin[[time_col]])  & clin[[time_col]] > 0, na.rm = TRUE)
  n_has_os_event  <- sum(!is.na(clin[[event_col]]), na.rm = TRUE)
  n_complete_surv <- sum(!is.na(clin[[time_col]])  & clin[[time_col]] > 0 &
                           !is.na(clin[[event_col]]), na.rm = TRUE)
  n_in_expression <- sum(clin$patient_id %in% expr_ids, na.rm = TRUE)

  # Dedup: check for duplicate patient_ids
  dup_ids <- clin$patient_id[duplicated(clin$patient_id)]
  n_dedup <- n_raw - length(dup_ids)

  # Final: complete survival + in expression
  clin_complete <- clin[!is.na(clin[[time_col]]) & clin[[time_col]] > 0 &
                           !is.na(clin[[event_col]]), ]
  clin_matched  <- clin_complete[clin_complete$patient_id %in% expr_ids, ]
  n_final       <- nrow(clin_matched)

  cat(sprintf("  N raw (clinical):              %d\n", n_raw))
  cat(sprintf("  N duplicate patient_ids:       %d\n", length(dup_ids)))
  cat(sprintf("  N after dedup:                 %d\n", n_dedup))
  cat(sprintf("  N with complete OS (time>0):   %d\n", n_complete_surv))
  cat(sprintf("  N with expression data:        %d\n", n_in_expression))
  cat(sprintf("  N final (clin+expr matched):   %d\n", n_final))

  log_df <- data.frame(
    step = c(
      "N_raw_clinical",
      "N_duplicate_patient_ids_removed",
      "N_after_dedup",
      "N_complete_OS_time_and_event",
      "N_with_expression_data",
      "N_final_training_set"
    ),
    n = c(
      n_raw,
      length(dup_ids),
      n_dedup,
      n_complete_surv,
      n_in_expression,
      n_final
    ),
    note = c(
      "All SCAN-B samples in clinical_FINAL.parquet (GSE96058, Option A)",
      "Technical replicates or duplicate IDs removed",
      "After removing duplicate patient_ids",
      "OS_Time > 0 AND OS_Event not NA",
      "patient_id present in expression matrix columns",
      "Used in script 05 model training"
    ),
    stringsAsFactors = FALSE
  )
  write.csv(log_df, "results/supp/train_inclusion_log.csv", row.names = FALSE)
  cat("\nInclusion log saved: results/supp/train_inclusion_log.csv\n")

  # Dedup log
  if (length(dup_ids) > 0) {
    dup_df <- clin[clin$patient_id %in% dup_ids, c("patient_id","OS_Time","OS_Event")]
    write.csv(dup_df, "results/supp/dedup_techrep_scanb.csv", row.names = FALSE)
    cat("Dedup log saved: results/supp/dedup_techrep_scanb.csv\n")
  } else {
    cat("No duplicate patient_ids found — dedup_techrep_scanb.csv not needed (N_dedup=0)\n")
    write.csv(data.frame(note="No duplicate patient_ids in SCAN-B clinical_FINAL.parquet"),
              "results/supp/dedup_techrep_scanb.csv", row.names = FALSE)
  }
}

cat("\nDone.\n")
