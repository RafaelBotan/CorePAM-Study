# =============================================================================
# SCRIPT: 02_harmonize_clinical_METABRIC.R
# PURPOSE: Clinical harmonization of METABRIC — correctly creates OS and DSS.
# PROJECT: Core-PAM (Memorial v6.1)
#
# INPUTS:
#   01_Base_Pura_CorePAM/RAW/METABRIC/brca_metabric/data_clinical_patient.txt
#   01_Base_Pura_CorePAM/RAW/METABRIC/brca_metabric/data_clinical_sample.txt
#
# OUTPUTS:
#   01_Base_Pura_CorePAM/PROCESSED/METABRIC/clinical_FINAL.parquet
#   01_docs/endpoint_mapping_templates/endpoint_mapping_METABRIC.csv
#
# RULES (Memorial v6.1 §7.1-7.2):
#   - Primary endpoint: DSS (dss_event=1 only cancer cause).
#   - OS as sensitivity.
#   - dss_time == os_time (time to death, regardless of cause).
#   - Other causes of death → censored at time of death (competing events).
#   - Time in months (days / 30.4375).
#   - time <= 0: DROP.
#   - clinical_FINAL.parquet is the primary source for all downstream scripts.
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "02_harmonize_clinical_METABRIC.R"

# =============================================================================
# 1) READ (cBioPortal: header lines start with "#")
# =============================================================================
metabric_dir <- file.path(raw_cohort("METABRIC"), "brca_metabric")

.read_cbio <- function(fname) {
  fpath <- file.path(metabric_dir, fname)
  if (!file.exists(fpath)) stop("File not found: ", fpath)
  # cBioPortal uses "#" for metadata lines → skip
  lines    <- readLines(fpath, warn = FALSE)
  skip_n   <- sum(startsWith(lines, "#"))
  read_tsv(fpath, skip = skip_n, show_col_types = FALSE,
           name_repair = "unique")
}

message("[02_METABRIC] Reading data_clinical_patient.txt ...")
clin_pat <- .read_cbio("data_clinical_patient.txt")

message("[02_METABRIC] Reading data_clinical_sample.txt ...")
clin_sam <- .read_cbio("data_clinical_sample.txt")

message(sprintf("[02_METABRIC] Patients: %d | Samples: %d",
                nrow(clin_pat), nrow(clin_sam)))
message("[02_METABRIC] Patient columns: ", paste(names(clin_pat), collapse = ", "))
message("[02_METABRIC] Sample columns:  ", paste(names(clin_sam), collapse = ", "))

# =============================================================================
# 2) JOIN patient x sample
# =============================================================================
# Join key: PATIENT_ID (cBioPortal standard)
clin_full <- clin_pat |>
  left_join(clin_sam, by = "PATIENT_ID", suffix = c("", ".sample"))

message(sprintf("[02_METABRIC] After join: %d rows", nrow(clin_full)))

# =============================================================================
# 3) COLUMN MAPPING cBioPortal → standardized names
#    Standard cBioPortal METABRIC names (check if they differ in current version)
# =============================================================================
# OS
OS_TIME_COL   <- "OS_MONTHS"      # tempo em meses (já em meses no cBioPortal)
OS_STATUS_COL <- "OS_STATUS"      # "0:LIVING" / "1:DECEASED"

# DSS (Disease-Specific Survival)
# cBioPortal METABRIC provides VITAL_STATUS with cause of death
DSS_CAUSE_COL <- "VITAL_STATUS"   # "Died of Disease" / "Died of Other Causes" / "Living"

# Covariates
AGE_COL    <- "AGE_AT_DIAGNOSIS"
ER_COL     <- "ER_IHC"            # or "ER_STATUS" depending on version
STAGE_COL  <- "TUMOR_STAGE"
SAMPLE_COL <- "SAMPLE_ID"
PAT_COL    <- "PATIENT_ID"

# =============================================================================
# 4) HARMONIZATION
# =============================================================================
n_raw <- nrow(clin_full)

out <- tibble(
  sample_id  = normalize_id(clin_full[[SAMPLE_COL]]),
  patient_id = normalize_id(clin_full[[PAT_COL]])
)

# OS: time already in months in cBioPortal METABRIC
out$os_time_months <- suppressWarnings(as.numeric(clin_full[[OS_TIME_COL]]))

# OS event: 1 = death any cause
raw_os <- tolower(trimws(clin_full[[OS_STATUS_COL]]))
out$os_event <- dplyr::case_when(
  grepl("^1|deceased|dead|died", raw_os) ~ 1L,
  grepl("^0|living|alive",       raw_os) ~ 0L,
  TRUE ~ NA_integer_
)

# DSS: dss_time == os_time; dss_event=1 only cancer cause
out$dss_time_months <- out$os_time_months   # mesmo tempo
raw_cause <- tolower(trimws(clin_full[[DSS_CAUSE_COL]]))
out$dss_event <- dplyr::case_when(
  grepl("died of disease|disease", raw_cause)       ~ 1L,  # cancer death
  grepl("died of other|other cause|living|alive",
        raw_cause)                                  ~ 0L,  # censored (competing)
  TRUE ~ NA_integer_
)

# CORE-A covariates
out$age       <- suppressWarnings(as.numeric(clin_full[[AGE_COL]]))
out$er_status <- dplyr::case_when(
  tolower(trimws(clin_full[[ER_COL]])) %in%
    c("positive", "positve", "pos", "1", "+") ~ "Positive",  # "positve" = cBioPortal typo
  tolower(trimws(clin_full[[ER_COL]])) %in%
    c("negative", "neg", "0", "-") ~ "Negative",
  TRUE ~ NA_character_
)
out$stage <- trimws(clin_full[[STAGE_COL]])

n_le0 <- sum(out$os_time_months <= 0, na.rm = TRUE)
n_na  <- sum(is.na(out$os_time_months))
out   <- out |> filter(os_time_months > 0, !is.na(os_time_months))

message(sprintf(
  "[02_METABRIC] N raw=%d | Removed time<=0: %d | NA time: %d | N final=%d",
  n_raw, n_le0, n_na, nrow(out)
))
message(sprintf(
  "[02_METABRIC] OS eventos: %d (%.1f%%) | DSS eventos: %d (%.1f%%)",
  sum(out$os_event,  na.rm = TRUE), 100 * mean(out$os_event,  na.rm = TRUE),
  sum(out$dss_event, na.rm = TRUE), 100 * mean(out$dss_event, na.rm = TRUE)
))
message(sprintf(
  "[02_METABRIC] ER+: %d | ER-: %d | ER NA: %d",
  sum(out$er_status == "Positive", na.rm = TRUE),
  sum(out$er_status == "Negative", na.rm = TRUE),
  sum(is.na(out$er_status))
))

# =============================================================================
# 5) SAVE FINAL PARQUET (primary source for all downstream scripts)
# =============================================================================
dest_dir  <- proc_cohort("METABRIC")
dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
out_path  <- file.path(dest_dir, "clinical_FINAL.parquet")

arrow::write_parquet(out, out_path)
h    <- sha256_file(out_path)
size <- file.info(out_path)$size / 1024^2
registry_append("METABRIC", "Clinical_FINAL", out_path, h,
                "INTEGRO", SCRIPT_NAME, size)
message(sprintf("[02_METABRIC] Saved: %s | %d samples | %.2f MB | SHA256: %s",
                out_path, nrow(out), size, h))

# =============================================================================
# 6) ENDPOINT MAPPING
# =============================================================================
ep_map <- tibble(
  cohort            = "METABRIC",
  endpoint_primary  = "DSS",
  endpoint_sensitivity = "OS",
  col_time_origin   = OS_TIME_COL,
  col_os_event      = OS_STATUS_COL,
  col_cause         = DSS_CAUSE_COL,
  unit_origin       = "months",
  unit_final        = "months",
  conversion        = "none",
  dss_event_1_rule  = "VITAL_STATUS contains 'died of disease'",
  dss_event_0_rule  = "all other (died other cause = competing censored at death time)",
  n_samples         = nrow(out),
  n_os_events       = sum(out$os_event,  na.rm = TRUE),
  n_dss_events      = sum(out$dss_event, na.rm = TRUE),
  script            = SCRIPT_NAME,
  timestamp         = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)
ep_path <- file.path(PATHS$docs, "endpoint_mapping_templates",
                     "endpoint_mapping_METABRIC.csv")
write_csv(ep_map, ep_path)
message("[02_METABRIC] Endpoint map: ", ep_path)

message("\n[02_METABRIC] Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message("Next step: scripts/03_expression_preprocess_METABRIC.R")
