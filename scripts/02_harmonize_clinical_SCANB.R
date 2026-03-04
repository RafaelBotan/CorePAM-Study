# =============================================================================
# SCRIPT: 02_harmonize_clinical_SCANB.R
# PURPOSE: Clinical harmonization of SCAN-B (GSE96058).
#          ALL 3069 samples are the SCAN-B training cohort (no internal split).
#          Decision: Option A, confirmed 2026-02-28 (Rafael Botan).
# PROJECT: Core-PAM (Memorial v6.1)
#
# INPUT:
#   01_Base_Pura_CorePAM/RAW/GSE96058/GSE96058_clinical_raw.rds  (pData GEO)
#
# OUTPUTS:
#   01_Base_Pura_CorePAM/PROCESSED/SCANB/clinical_FINAL.parquet
#   01_docs/endpoint_mapping_templates/endpoint_mapping_SCANB.csv
#
# RULES (Memorial v6.1):
#   - Time in months (days / 30.4375).
#   - time <= 0: DROP (report count).
#   - event OS: 1 = death any cause; 0 = censored.
#   - No imputation: report effective N per variable.
#   - ALL GSE96058 samples = SCAN-B training (no train/val split in pData).
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "02_harmonize_clinical_SCANB.R"

# Skip if output already exists (unless FORCE_RERUN=TRUE)
FORCE    <- as.logical(Sys.getenv("FORCE_RERUN", "FALSE"))
out_path <- file.path(proc_cohort("SCANB"), "clinical_FINAL.parquet")
if (!FORCE && file.exists(out_path)) {
  message(sprintf("[02_SCANB] Output exists — skipping. Set FORCE_RERUN=TRUE to rerun.\n  %s", out_path))
  quit(save = "no", status = 0)
}

# =============================================================================
# 1) STRICT READ
# =============================================================================
raw_path <- file.path(raw_cohort("GSE96058"), "GSE96058_clinical_raw.rds")
pheno    <- strict_rds(raw_path)

message(sprintf("[02_SCANB] pData loaded: %d samples x %d columns",
                nrow(pheno), ncol(pheno)))
message("[02_SCANB] ALL samples assigned to SCAN-B training cohort (Option A)")

# =============================================================================
# 2) COLUMN MAPPING (verified against actual GSE96058 pData 2026-02-28)
#    - os_time_raw: "overall survival days:ch1" (values ~2000+ days)
#    - os_event:    "overall survival event:ch1" (0/1)
#    - age:         "age at diagnosis:ch1"
#    - er_status:   "er status:ch1"  (values: 0, 1, NA)
#    - NOTE: no tumor stage column exists in this pData
# =============================================================================
COL_MAP <- list(
  sample_id   = "geo_accession",
  patient_id  = "title",
  os_time_raw = "overall survival days:ch1",
  os_event    = "overall survival event:ch1",
  age         = "age at diagnosis:ch1",
  er_status   = "er status:ch1"
)

# Verify required columns exist
missing_cols <- setdiff(unlist(COL_MAP), names(pheno))
if (length(missing_cols) > 0) {
  stop(sprintf(
    "[02_SCANB] Missing expected columns in pData:\n  %s\nAvailable columns:\n  %s",
    paste(missing_cols, collapse = "\n  "),
    paste(names(pheno), collapse = "\n  ")
  ))
}
message("[02_SCANB] Column mapping verified — all required columns present")

# =============================================================================
# 3) HARMONIZATION
# =============================================================================
n_raw <- nrow(pheno)

clin <- tibble(
  sample_id  = normalize_id(pheno[[COL_MAP$sample_id]]),
  patient_id = normalize_id(pheno[[COL_MAP$patient_id]])
)

# Survival time — column is in days; convert to months
raw_time <- suppressWarnings(as.numeric(pheno[[COL_MAP$os_time_raw]]))
time_median <- median(raw_time, na.rm = TRUE)
message(sprintf("[02_SCANB] Time column: '%s' | median=%.1f days → converting to months",
                COL_MAP$os_time_raw, time_median))
clin$os_time_months <- raw_time / FREEZE$time_unit_divisor

# OS event: standardize to 0/1
raw_event <- tolower(trimws(pheno[[COL_MAP$os_event]]))
clin$os_event <- dplyr::case_when(
  raw_event %in% c("1", "dead", "deceased", "yes", "true")  ~ 1L,
  raw_event %in% c("0", "alive", "living",  "no",  "false") ~ 0L,
  TRUE ~ NA_integer_
)

# Age
clin$age <- suppressWarnings(as.numeric(pheno[[COL_MAP$age]]))

# ER status: values in pData are 0/1 (integer), no text labels
clin$er_status <- dplyr::case_when(
  tolower(trimws(pheno[[COL_MAP$er_status]])) %in% c("1", "positive", "pos", "+") ~ "Positive",
  tolower(trimws(pheno[[COL_MAP$er_status]])) %in% c("0", "negative", "neg", "-") ~ "Negative",
  TRUE ~ NA_character_
)

# Grade: from "nhg:ch1" (Nottingham Histologic Grade, values G1/G2/G3)
clin$grade <- suppressWarnings(as.integer(sub("^G", "", pheno[["nhg:ch1"]])))
message(sprintf("[02_SCANB] Grade: 1=%d, 2=%d, 3=%d, NA=%d",
                sum(clin$grade == 1L, na.rm = TRUE),
                sum(clin$grade == 2L, na.rm = TRUE),
                sum(clin$grade == 3L, na.rm = TRUE),
                sum(is.na(clin$grade))))

# HER2 status: from "her2 status:ch1" (values 0/1)
clin$her2_status <- suppressWarnings(as.integer(pheno[["her2 status:ch1"]]))
message(sprintf("[02_SCANB] HER2: 0=%d, 1=%d, NA=%d",
                sum(clin$her2_status == 0L, na.rm = TRUE),
                sum(clin$her2_status == 1L, na.rm = TRUE),
                sum(is.na(clin$her2_status))))

# =============================================================================
# 4) FILTER: remove time <= 0 and time = NA
# =============================================================================
n_le0 <- sum(!is.na(clin$os_time_months) & clin$os_time_months <= 0)
clin  <- clin |> filter(!is.na(os_time_months) & os_time_months > 0)

message(sprintf(
  "[02_SCANB] N raw=%d | Removed time<=0: %d | N final=%d",
  n_raw, n_le0, nrow(clin)
))
message(sprintf(
  "[02_SCANB] OS events: %d (%.1f%%) | ER+ %d | ER- %d | ER NA %d | Age NA %d",
  sum(clin$os_event, na.rm = TRUE),
  100 * mean(clin$os_event, na.rm = TRUE),
  sum(clin$er_status == "Positive", na.rm = TRUE),
  sum(clin$er_status == "Negative", na.rm = TRUE),
  sum(is.na(clin$er_status)),
  sum(is.na(clin$age))
))

# =============================================================================
# 5) SAVE PARQUET + REGISTRY
# =============================================================================
dest_dir <- proc_cohort("SCANB")
dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)

old_warn <- getOption("warn"); options(warn = 0)
arrow::write_parquet(clin, out_path)
options(warn = old_warn)

h    <- sha256_file(out_path)
size <- file.info(out_path)$size / 1024^2
registry_append("SCANB", "Clinical_FINAL", out_path, h,
                "INTEGRO", SCRIPT_NAME, size)
message(sprintf("[02_SCANB] Saved: %s | %d samples | %.2f MB | SHA256: %s",
                out_path, nrow(clin), size, h))

# =============================================================================
# 6) ENDPOINT MAPPING (audit dictionary)
# =============================================================================
endpt_map <- tibble(
  cohort           = "SCANB",
  endpoint         = "OS",
  col_time_origin  = COL_MAP$os_time_raw,
  col_event_origin = COL_MAP$os_event,
  unit_origin      = "days",
  unit_final       = "months",
  conversion       = "/ 30.4375",
  event_1_means    = "death_any_cause",
  n_samples        = nrow(clin),
  n_events         = sum(clin$os_event, na.rm = TRUE),
  script           = SCRIPT_NAME,
  timestamp        = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)
endpt_dir <- file.path(PATHS$docs, "endpoint_mapping_templates")
dir.create(endpt_dir, showWarnings = FALSE, recursive = TRUE)
endpt_path <- file.path(endpt_dir, "endpoint_mapping_SCANB.csv")
write_csv(endpt_map, endpt_path)
message(sprintf("[02_SCANB] Endpoint map saved: %s", endpt_path))

message("\n[02_SCANB] Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message("Next step: scripts/03_expression_preprocess_SCANB.R")
