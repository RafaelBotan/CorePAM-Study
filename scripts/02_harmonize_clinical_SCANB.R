# =============================================================================
# SCRIPT: 02_harmonize_clinical_SCANB.R
# PURPOSE: Clinical harmonization of SCAN-B (GSE96058) and deterministic split
#          training (SCAN-B) vs validation (GSE96058).
# PROJETO: Core-PAM (Memorial v6.1)
#
# INPUT:
#   01_Base_Pura_CorePAM/RAW/GSE96058/GSE96058_clinical_raw.rds  (pData GEO)
#
# OUTPUTS:
#   01_Base_Pura_CorePAM/PROCESSED/SCANB/clinical_FINAL.parquet
#   01_Base_Pura_CorePAM/PROCESSED/GSE96058/clinical_FINAL.parquet
#   01_docs/endpoint_mapping_templates/endpoint_mapping_SCANB.csv
#   01_docs/endpoint_mapping_templates/endpoint_mapping_GSE96058.csv
#
# RULES (Memorial v6.1):
#   - Time in months (days / 30.4375).
#   - time <= 0: DROP (report count).
#   - event OS: 1 = death any cause; 0 = censored.
#   - No imputation: report effective N per variable.
#   - Leakage-proof: training and validation IDs are disjoint sets.
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "02_harmonize_clinical_SCANB.R"

# =============================================================================
# 1) STRICT READ
# =============================================================================
raw_path <- file.path(raw_cohort("GSE96058"), "GSE96058_clinical_raw.rds")
pheno    <- strict_rds(raw_path)

message(sprintf("[02_SCANB] pData loaded: %d samples x %d columns",
                nrow(pheno), ncol(pheno)))

# Inspect available columns (log for audit)
message("[02_SCANB] Available columns:")
print(names(pheno))

# =============================================================================
# 2) IDENTIFY TRAINING / VALIDATION SPLIT COLUMN
#    GEO GSE96058 contains a column indicating "training" vs "validation".
#    Exact name may vary — identify automatically and register.
# =============================================================================
.split_candidates <- grep(
  "train|valid|group|cohort|set",
  tolower(names(pheno)), value = TRUE
)
message("[02_SCANB] Split column candidates: ",
        paste(names(pheno)[tolower(names(pheno)) %in% .split_candidates],
              collapse = ", "))

# Adjust SPLIT_COL if the automatic name is not correct
SPLIT_COL <- names(pheno)[
  tolower(names(pheno)) %in% .split_candidates
][1]

if (is.na(SPLIT_COL) || length(SPLIT_COL) == 0) {
  stop(
    "[02_SCANB] Split column not detected automatically.\n",
    "Inspect names(pheno) and define SPLIT_COL manually."
  )
}
message("[02_SCANB] Using split column: ", SPLIT_COL)
message("[02_SCANB] Unique values: ", paste(unique(pheno[[SPLIT_COL]]), collapse = " | "))

# =============================================================================
# 3) COLUMN MAPPING (adjust names according to actual pData)
#    These are the EXPECTED names in the GSE96058 pData.
#    If they differ, adjust the map below.
# =============================================================================
COL_MAP <- list(
  sample_id   = "geo_accession",          # GSM ID
  patient_id  = "title",                  # patient ID (SCAN-B barcode)
  os_time_raw = "overall survival (months):ch1",   # time in months (verify unit)
  os_event    = "overall survival event:ch1",       # 0/1 ou dead/alive
  age         = "age at diagnosis:ch1",
  er_status   = "er status:ch1",
  stage       = "tumor stage:ch1"
)

# =============================================================================
# 4) HARMONIZATION FUNCTION (applied to training and validation)
# =============================================================================
harmonize_scanb <- function(df, cohort_label) {
  n_raw <- nrow(df)

  out <- tibble(
    sample_id  = normalize_id(df[[COL_MAP$sample_id]]),
    patient_id = normalize_id(df[[COL_MAP$patient_id]])
  )

  # Survival time — check if already in months or in days
  raw_time <- suppressWarnings(as.numeric(df[[COL_MAP$os_time_raw]]))

  # Heuristic: if median > 200 assumes days → convert
  if (median(raw_time, na.rm = TRUE) > 200) {
    message(sprintf(
      "[02_%s] Time in days detected (median=%.1f). Converting to months.",
      cohort_label, median(raw_time, na.rm = TRUE)
    ))
    out$os_time_months <- raw_time / FREEZE$time_unit_divisor
  } else {
    message(sprintf(
      "[02_%s] Time already in months (median=%.1f).",
      cohort_label, median(raw_time, na.rm = TRUE)
    ))
    out$os_time_months <- raw_time
  }

  # OS event: standardize to 0/1
  raw_event <- tolower(trimws(df[[COL_MAP$os_event]]))
  out$os_event <- dplyr::case_when(
    raw_event %in% c("1", "dead", "deceased", "yes", "true")  ~ 1L,
    raw_event %in% c("0", "alive", "living",  "no",  "false") ~ 0L,
    TRUE ~ NA_integer_
  )

  # CORE-A covariates
  out$age      <- suppressWarnings(as.numeric(df[[COL_MAP$age]]))
  out$er_status <- dplyr::case_when(
    tolower(trimws(df[[COL_MAP$er_status]])) %in% c("positive", "pos", "1", "+") ~ "Positive",
    tolower(trimws(df[[COL_MAP$er_status]])) %in% c("negative", "neg", "0", "-") ~ "Negative",
    TRUE ~ NA_character_
  )
  out$stage <- trimws(df[[COL_MAP$stage]])

  n_le0 <- sum(out$os_time_months <= 0, na.rm = TRUE)
  out   <- out |> filter(os_time_months > 0)
  n_na  <- sum(is.na(out$os_time_months))
  out   <- out |> filter(!is.na(os_time_months))

  message(sprintf(
    "[02_%s] N raw=%d | Removed time<=0: %d | Removed time=NA: %d | N final=%d",
    cohort_label, n_raw, n_le0, n_na, nrow(out)
  ))
  message(sprintf(
    "[02_%s] Eventos OS: %d (%.1f%%) | ER+ %d | ER- %d | ER NA %d",
    cohort_label,
    sum(out$os_event, na.rm = TRUE),
    100 * mean(out$os_event, na.rm = TRUE),
    sum(out$er_status == "Positive", na.rm = TRUE),
    sum(out$er_status == "Negative", na.rm = TRUE),
    sum(is.na(out$er_status))
  ))

  out
}

# =============================================================================
# 5) SPLIT AND HARMONIZATION
# =============================================================================
train_mask <- grepl("train", tolower(pheno[[SPLIT_COL]]))
valid_mask <- !train_mask

pheno_train <- pheno[train_mask, ]
pheno_valid <- pheno[valid_mask, ]

message(sprintf("[02_SCANB] Split: training=%d | validation=%d",
                nrow(pheno_train), nrow(pheno_valid)))

# Check disjunction (leakage-proof)
ids_train <- normalize_id(pheno_train[[COL_MAP$sample_id]])
ids_valid <- normalize_id(pheno_valid[[COL_MAP$sample_id]])
overlap   <- intersect(ids_train, ids_valid)
if (length(overlap) > 0) {
  stop(sprintf("[02_SCANB] LEAKAGE: %d IDs overlapping between training and validation: %s",
               length(overlap), paste(overlap[1:min(5, length(overlap))], collapse = ",")))
}
message("[02_SCANB] Leakage check: PASS (overlap = 0)")

clin_train <- harmonize_scanb(pheno_train, "SCANB")
clin_valid <- harmonize_scanb(pheno_valid, "GSE96058")

# =============================================================================
# 6) SAVE FINAL PARQUETS
# =============================================================================
.save_clinical <- function(df, cohort) {
  dest_dir <- proc_cohort(cohort)
  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
  out_path <- file.path(dest_dir, "clinical_FINAL.parquet")
  arrow::write_parquet(df, out_path)

  h    <- sha256_file(out_path)
  size <- file.info(out_path)$size / 1024^2
  registry_append(cohort, "Clinical_FINAL", out_path, h,
                  "INTEGRO", SCRIPT_NAME, size)
  message(sprintf("[02_%s] Saved: %s | %d samples | %.2f MB | SHA256: %s",
                  cohort, out_path, nrow(df), size, h))
}

.save_clinical(clin_train, "SCANB")
.save_clinical(clin_valid, "GSE96058")

# =============================================================================
# 7) ENDPOINT MAPPING (audit dictionary)
# =============================================================================
.save_endpoint_map <- function(cohort, n_samples, n_events,
                               time_col_origin, event_col_origin,
                               time_unit_origin) {
  out <- tibble(
    cohort           = cohort,
    endpoint         = "OS",
    col_time_origin  = time_col_origin,
    col_event_origin = event_col_origin,
    unit_origin      = time_unit_origin,
    unit_final       = "months",
    conversion       = ifelse(time_unit_origin == "days",
                              "/ 30.4375", "none"),
    event_1_means    = "death_any_cause",
    n_samples        = n_samples,
    n_events         = n_events,
    script           = SCRIPT_NAME,
    timestamp        = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  dest <- file.path(PATHS$docs, "endpoint_mapping_templates",
                    sprintf("endpoint_mapping_%s.csv", cohort))
  write_csv(out, dest)
  message(sprintf("[02_%s] Endpoint map saved: %s", cohort, dest))
}

.save_endpoint_map("SCANB", nrow(clin_train),
                   sum(clin_train$os_event, na.rm = TRUE),
                   COL_MAP$os_time_raw, COL_MAP$os_event, "months_or_days")
.save_endpoint_map("GSE96058", nrow(clin_valid),
                   sum(clin_valid$os_event, na.rm = TRUE),
                   COL_MAP$os_time_raw, COL_MAP$os_event, "months_or_days")

message("\n[02_SCANB] Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message("Next step: scripts/03_expression_preprocess_SCANB.R")
