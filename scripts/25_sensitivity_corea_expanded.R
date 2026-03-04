# =============================================================================
# SCRIPT: 25_sensitivity_corea_expanded.R
# PURPOSE: Sensitivity analysis — CORE-A Expanded (age + ER + T-stage + N-status)
#          Tests whether CorePAM score remains significant after adjusting for
#          anatomical staging variables (pT, pN) beyond CORE-A (age + ER).
# PROJECT: Core-PAM (Memorial CorePAM)
#
# OUTPUT:
#   results/supp/EXP_corea_expanded_results.csv
# =============================================================================
source("scripts/00_setup.R")
SCRIPT_NAME <- "25_sensitivity_corea_expanded.R"

suppressPackageStartupMessages({
  library(survival)
  library(arrow)
})

OUT_DIR <- file.path(PATHS$results$supp)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

message("[25_SENS] ========== CORE-A EXPANDED ==========")

# --------------------------------------------------------------------------
# 1) Extract T-stage + N-status from RAW clinical files
# --------------------------------------------------------------------------

# --- SCAN-B ---
scanb_raw <- readRDS(file.path(PATHS$raw, "GSE96058/GSE96058_clinical_raw.rds"))
scanb_ar  <- read_parquet(file.path(PATHS$processed, "SCANB/analysis_ready.parquet"))

# Tumor size (mm) -> T-stage (suppress NA coercion warnings)
old_warn <- getOption("warn"); options(warn = 0)
scanb_size <- as.numeric(scanb_raw$`tumor size:ch1`)
scanb_t <- ifelse(scanb_size <= 20, "T1",
           ifelse(scanb_size <= 50, "T2",
           ifelse(scanb_size > 50,  "T3", NA_character_)))

# Node status (binary)
scanb_node_raw <- scanb_raw$`lymph node status:ch1`
scanb_n <- ifelse(scanb_node_raw == "NodeNegative", "N0",
           ifelse(scanb_node_raw == "NodePositive", "N+", NA_character_))

options(warn = old_warn)

scanb_staging <- data.frame(
  sample_id = scanb_raw$geo_accession,
  t_stage = scanb_t,
  n_status = scanb_n,
  stringsAsFactors = FALSE
)
message(sprintf("[25_SENS] SCAN-B staging: T=%d/%d (%.1f%%), N=%d/%d (%.1f%%)",
                sum(!is.na(scanb_t)), nrow(scanb_raw), 100*mean(!is.na(scanb_t)),
                sum(!is.na(scanb_n)), nrow(scanb_raw), 100*mean(!is.na(scanb_n))))

# --- TCGA-BRCA ---
tcga_raw <- readRDS(file.path(PATHS$raw, "TCGA_BRCA/TCGA_BRCA_clinical_raw.rds"))
tcga_ar  <- read_parquet(file.path(PATHS$processed, "TCGA_BRCA/analysis_ready.parquet"))

# T-stage from ajcc_pathologic_t
tcga_t_raw <- tcga_raw$ajcc_pathologic_t
tcga_t <- ifelse(grepl("^T1", tcga_t_raw), "T1",
          ifelse(grepl("^T2", tcga_t_raw), "T2",
          ifelse(grepl("^T3", tcga_t_raw), "T3",
          ifelse(grepl("^T4", tcga_t_raw), "T4",
          NA_character_))))

# N-status from ajcc_pathologic_n
tcga_n_raw <- tcga_raw$ajcc_pathologic_n
tcga_n <- ifelse(grepl("^N0", tcga_n_raw), "N0",
          ifelse(grepl("^N[1-3]", tcga_n_raw), "N+",
          NA_character_))

# Find patient ID column (TCGA uses submitter_id, not patient_id)
tcga_id_col <- intersect(names(tcga_raw), c("patient_id", "submitter_id", "bcr_patient_barcode", "case_id"))
tcga_patient_ids <- tcga_raw[[tcga_id_col[1]]]
message(sprintf("[25_SENS] TCGA ID column: %s (N=%d)", tcga_id_col[1], length(tcga_patient_ids)))

tcga_staging <- data.frame(
  patient_id = tcga_patient_ids,
  t_stage = tcga_t,
  n_status = tcga_n,
  stringsAsFactors = FALSE
)
message(sprintf("[25_SENS] TCGA-BRCA staging: T=%d/%d (%.1f%%), N=%d/%d (%.1f%%)",
                sum(!is.na(tcga_t)), nrow(tcga_raw), 100*mean(!is.na(tcga_t)),
                sum(!is.na(tcga_n)), nrow(tcga_raw), 100*mean(!is.na(tcga_n))))

# --- METABRIC ---
metabric_clin_p <- file.path(PATHS$raw, "METABRIC/brca_metabric/data_clinical_patient.txt")
metabric_clin_s <- file.path(PATHS$raw, "METABRIC/brca_metabric/data_clinical_sample.txt")
metabric_ar <- read_parquet(file.path(PATHS$processed, "METABRIC/analysis_ready.parquet"))

# Read cBioPortal clinical files (skip comment lines starting with #)
old_warn <- getOption("warn"); options(warn = 0)
mb_pat <- read.delim(metabric_clin_p, comment.char = "#", stringsAsFactors = FALSE)
mb_sam <- read.delim(metabric_clin_s, comment.char = "#", stringsAsFactors = FALSE)
options(warn = old_warn)

# Tumor size (mm) -> T-stage (suppress NA coercion warnings)
old_warn <- getOption("warn"); options(warn = 0)
mb_size <- as.numeric(mb_sam$TUMOR_SIZE)
mb_t <- ifelse(mb_size <= 20, "T1",
        ifelse(mb_size <= 50, "T2",
        ifelse(mb_size > 50,  "T3", NA_character_)))

# Node status (binary) from lymph node count
mb_nodes <- as.numeric(mb_pat$LYMPH_NODES_EXAMINED_POSITIVE)
options(warn = old_warn)
mb_n <- ifelse(!is.na(mb_nodes) & mb_nodes == 0, "N0",
        ifelse(!is.na(mb_nodes) & mb_nodes > 0,  "N+", NA_character_))

metabric_staging <- data.frame(
  patient_id = mb_pat$PATIENT_ID,
  t_stage = mb_t[match(mb_pat$PATIENT_ID, mb_sam$PATIENT_ID)],
  n_status = mb_n,
  stringsAsFactors = FALSE
)
message(sprintf("[25_SENS] METABRIC staging: T=%d/%d (%.1f%%), N=%d/%d (%.1f%%)",
                sum(!is.na(metabric_staging$t_stage)), nrow(metabric_staging),
                100*mean(!is.na(metabric_staging$t_stage)),
                sum(!is.na(metabric_staging$n_status)), nrow(metabric_staging),
                100*mean(!is.na(metabric_staging$n_status))))

# --- GSE20685 ---
gse_raw <- readRDS(file.path(PATHS$raw, "GSE20685/GSE20685_clinical_raw.rds"))
gse_ar  <- read_parquet(file.path(PATHS$processed, "GSE20685/analysis_ready.parquet"))

# T-stage from t_stage:ch1
gse_t_raw <- gse_raw$`t_stage:ch1`
gse_t <- ifelse(gse_t_raw %in% c("1", "1c"), "T1",
         ifelse(gse_t_raw == "2", "T2",
         ifelse(gse_t_raw == "3", "T3",
         ifelse(gse_t_raw == "4", "T4", NA_character_))))

# N-status from n_stage:ch1
gse_n_raw <- gse_raw$`n_stage:ch1`
gse_n <- ifelse(gse_n_raw == "0", "N0",
         ifelse(gse_n_raw %in% c("1","2","3"), "N+", NA_character_))

gse_staging <- data.frame(
  sample_id = gse_raw$geo_accession,
  t_stage = gse_t,
  n_status = gse_n,
  stringsAsFactors = FALSE
)
message(sprintf("[25_SENS] GSE20685 staging: T=%d/%d (%.1f%%), N=%d/%d (%.1f%%)",
                sum(!is.na(gse_t)), nrow(gse_raw), 100*mean(!is.na(gse_t)),
                sum(!is.na(gse_n)), nrow(gse_raw), 100*mean(!is.na(gse_n))))

# --------------------------------------------------------------------------
# 2) Merge staging with analysis_ready and run survival models
# --------------------------------------------------------------------------
COHORTS <- list(
  SCANB    = list(ar = scanb_ar,    staging = scanb_staging,    join_by = "sample_id",
                  time = "os_time_months", event = "os_event"),
  TCGA_BRCA = list(ar = tcga_ar,   staging = tcga_staging,     join_by = "patient_id",
                  time = "os_time_months", event = "os_event"),
  METABRIC = list(ar = metabric_ar, staging = metabric_staging, join_by = "patient_id",
                  time = "dss_time_months", event = "dss_event"),
  GSE20685 = list(ar = gse_ar,     staging = gse_staging,      join_by = "sample_id",
                  time = "os_time_months", event = "os_event")
)

results <- list()

for (coh_name in names(COHORTS)) {
  coh <- COHORTS[[coh_name]]
  ar <- as.data.frame(coh$ar)
  stg <- coh$staging

  # Merge
  jby <- coh$join_by
  if (jby %in% names(ar) && jby %in% names(stg)) {
    merged <- merge(ar, stg, by = jby, all.x = TRUE)
  } else if ("patient_id" %in% names(ar) && "patient_id" %in% names(stg)) {
    merged <- merge(ar, stg, by = "patient_id", all.x = TRUE)
  } else {
    message(sprintf("[25_SENS] %s: cannot merge, skipping", coh_name))
    next
  }

  time_col  <- coh$time
  event_col <- coh$event

  # Filter valid survival data
  merged <- merged[!is.na(merged[[time_col]]) & merged[[time_col]] > 0 &
                   !is.na(merged[[event_col]]), ]

  message(sprintf("[25_SENS] %s: N=%d, T available=%d, N available=%d",
                  coh_name, nrow(merged),
                  sum(!is.na(merged$t_stage)),
                  sum(!is.na(merged$n_status))))

  # Factorize
  merged$t_stage  <- factor(merged$t_stage, levels = c("T1","T2","T3","T4"))
  merged$n_status <- factor(merged$n_status, levels = c("N0","N+"))

  old_warn <- getOption("warn"); options(warn = 0)

  surv_obj <- Surv(merged[[time_col]], merged[[event_col]])

  # Model 1: Current CORE-A (age + ER)
  has_er <- "er_status" %in% names(merged) && sum(!is.na(merged$er_status)) > 10
  if (has_er) {
    fit_current <- coxph(surv_obj ~ score_z + age + er_status, data = merged)
  } else {
    fit_current <- coxph(surv_obj ~ score_z + age, data = merged)
  }

  # Model 2: Expanded CORE-A (age + ER + T + N)
  # Only use complete cases for T and N
  merged_tn <- merged[!is.na(merged$t_stage) & !is.na(merged$n_status), ]
  surv_tn <- Surv(merged_tn[[time_col]], merged_tn[[event_col]])

  if (has_er) {
    fit_expanded <- coxph(surv_tn ~ score_z + age + er_status + t_stage + n_status, data = merged_tn)
    fit_clinical_only <- coxph(surv_tn ~ age + er_status + t_stage + n_status, data = merged_tn)
    fit_clinical_plus <- coxph(surv_tn ~ age + er_status + t_stage + n_status + score_z, data = merged_tn)
  } else {
    fit_expanded <- coxph(surv_tn ~ score_z + age + t_stage + n_status, data = merged_tn)
    fit_clinical_only <- coxph(surv_tn ~ age + t_stage + n_status, data = merged_tn)
    fit_clinical_plus <- coxph(surv_tn ~ age + t_stage + n_status + score_z, data = merged_tn)
  }

  # Extract CorePAM HR from expanded model
  coef_exp <- summary(fit_expanded)$coefficients
  hr_exp <- exp(coef_exp["score_z", "coef"])
  p_exp  <- coef_exp["score_z", "Pr(>|z|)"]
  ci_lo <- exp(coef_exp["score_z", "coef"] - 1.96 * coef_exp["score_z", "se(coef)"])
  ci_hi <- exp(coef_exp["score_z", "coef"] + 1.96 * coef_exp["score_z", "se(coef)"])

  # C-index: clinical-only vs clinical+score (expanded)
  c_clin <- concordance(fit_clinical_only)$concordance
  c_plus <- concordance(fit_clinical_plus)$concordance
  delta_c_expanded <- c_plus - c_clin

  # Current CORE-A comparison (for reference)
  coef_cur <- summary(fit_current)$coefficients
  hr_cur <- exp(coef_cur["score_z", "coef"])
  p_cur  <- coef_cur["score_z", "Pr(>|z|)"]

  options(warn = old_warn)

  results[[coh_name]] <- data.frame(
    cohort = coh_name,
    n_total = nrow(merged),
    n_tn_complete = nrow(merged_tn),
    hr_current_corea = round(hr_cur, 4),
    p_current_corea = signif(p_cur, 4),
    hr_expanded_corea = round(hr_exp, 4),
    hr_exp_lo95 = round(ci_lo, 4),
    hr_exp_hi95 = round(ci_hi, 4),
    p_expanded_corea = signif(p_exp, 4),
    c_clinical_expanded = round(c_clin, 4),
    c_clinical_plus_score = round(c_plus, 4),
    delta_c_expanded = round(delta_c_expanded, 4),
    corea_vars_current = if(has_er) "age+ER" else "age",
    corea_vars_expanded = if(has_er) "age+ER+T+N" else "age+T+N",
    stringsAsFactors = FALSE
  )

  message(sprintf("[25_SENS] %s: HR_current=%.3f (p=%.4g) → HR_expanded=%.3f (p=%.4g) | ΔC_expanded=%.4f",
                  coh_name, hr_cur, p_cur, hr_exp, p_exp, delta_c_expanded))
}

# --------------------------------------------------------------------------
# 3) Save results
# --------------------------------------------------------------------------
results_df <- do.call(rbind, results)
out_path <- file.path(OUT_DIR, "EXP_corea_expanded_results.csv")
write.csv(results_df, out_path, row.names = FALSE)
message(sprintf("[25_SENS] Results saved to %s", out_path))
message("[25_SENS] ========== DONE ==========")
