# =============================================================================
# SCRIPT: 02_harmonize_clinical_GSE20685.R
# PURPOSE: Clinical harmonization of GSE20685 (Taiwan; Affymetrix microarray).
# PROJETO: Core-PAM (Memorial v6.1)
#
# INPUT:
#   01_Base_Pura_CorePAM/RAW/GSE20685/GSE20685_clinical_raw.rds  (pData GEO)
#
# OUTPUT:
#   01_Base_Pura_CorePAM/PROCESSED/GSE20685/clinical_FINAL.parquet
#   01_docs/endpoint_mapping_templates/endpoint_mapping_GSE20685.csv
#
# REFERENCE: Lu et al. 2012 — 327 patients with BC (BC specific-survival)
#
# RULES (Memorial v6.1):
#   - Endpoint: OS (verify if dataset provides OS or RFS/BCSS).
#   - Time in months (days / 30.4375 if needed).
#   - time <= 0: DROP.
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "02_harmonize_clinical_GSE20685.R"

# =============================================================================
# 1) STRICT READ
# =============================================================================
raw_path <- file.path(raw_cohort("GSE20685"), "GSE20685_clinical_raw.rds")
pheno    <- strict_rds(raw_path)

message(sprintf("[02_GSE20685] pData: %d samples x %d columns",
                nrow(pheno), ncol(pheno)))
message("[02_GSE20685] Available columns:")
print(names(pheno))
message("[02_GSE20685] Sample of first rows (characteristics columns):")
char_cols <- grep("characteristics", names(pheno), value = TRUE)
if (length(char_cols) > 0) print(head(pheno[, char_cols, drop = FALSE], 3))

# =============================================================================
# 2) MAPPING (GSE20685 — adjust after inspecting actual pData)
#    GSE20685 typically has clinical data in "characteristics_ch1.*" columns
#    Inspect and adjust the names below.
# =============================================================================
# Automatic identification of time and event columns
.find_col <- function(pheno, patterns) {
  candidates <- names(pheno)[
    grepl(paste(patterns, collapse = "|"), tolower(names(pheno)))
  ]
  if (length(candidates) == 0) return(NA_character_)
  candidates[1]
}

TIME_COL  <- .find_col(pheno, c("survival", "time", "months", "days",
                                "os.time", "os_time", "fu.time"))
EVENT_COL <- .find_col(pheno, c("status", "event", "death", "dead",
                                "os.event", "os_event", "censor"))
AGE_COL   <- .find_col(pheno, c("age"))
ER_COL    <- .find_col(pheno, c("er ", "er_", "estrogen"))

message(sprintf(
  "[02_GSE20685] Detected columns — time: '%s' | event: '%s' | age: '%s' | ER: '%s'",
  TIME_COL, EVENT_COL, AGE_COL, ER_COL
))

if (is.na(TIME_COL) || is.na(EVENT_COL)) {
  message("[02_GSE20685] WARNING: time/event columns not detected automatically.")
  message("  Inspect pheno and define TIME_COL / EVENT_COL manually.")
  message("  Available columns: ", paste(names(pheno), collapse = " | "))
  stop("Define TIME_COL and EVENT_COL manually before proceeding.")
}

# =============================================================================
# 3) HARMONIZATION
# =============================================================================
n_raw <- nrow(pheno)

out <- tibble(
  sample_id  = normalize_id(pheno$geo_accession),
  patient_id = normalize_id(pheno$title)
)

# Survival time
raw_time <- suppressWarnings(as.numeric(pheno[[TIME_COL]]))

# Unit heuristic: if median > 200 → days; else → months
if (!all(is.na(raw_time)) && median(raw_time, na.rm = TRUE) > 200) {
  message(sprintf(
    "[02_GSE20685] Time detected in DAYS (median=%.1f). Converting to months.",
    median(raw_time, na.rm = TRUE)
  ))
  out$os_time_months <- raw_time / FREEZE$time_unit_divisor
  time_unit_origin   <- "days"
} else {
  message(sprintf(
    "[02_GSE20685] Time in MONTHS (median=%.1f).",
    median(raw_time, na.rm = TRUE)
  ))
  out$os_time_months <- raw_time
  time_unit_origin   <- "months"
}

# OS event
raw_event <- tolower(trimws(pheno[[EVENT_COL]]))
out$os_event <- dplyr::case_when(
  raw_event %in% c("1", "dead", "deceased", "yes", "true",  "death") ~ 1L,
  raw_event %in% c("0", "alive", "living",  "no",  "false", "censored") ~ 0L,
  TRUE ~ NA_integer_
)

# Age
if (!is.na(AGE_COL)) {
  out$age <- suppressWarnings(as.numeric(pheno[[AGE_COL]]))
} else {
  out$age <- NA_real_
}

# ER status
if (!is.na(ER_COL)) {
  out$er_status <- dplyr::case_when(
    tolower(trimws(pheno[[ER_COL]])) %in%
      c("positive", "pos", "1", "+", "er+") ~ "Positive",
    tolower(trimws(pheno[[ER_COL]])) %in%
      c("negative", "neg", "0", "-", "er-") ~ "Negative",
    TRUE ~ NA_character_
  )
} else {
  out$er_status <- NA_character_
  message("[02_GSE20685] WARNING: ER status not found — marked as NA.")
}

# DROP: time <= 0
n_le0 <- sum(out$os_time_months <= 0, na.rm = TRUE)
n_na  <- sum(is.na(out$os_time_months))
out   <- out |> filter(os_time_months > 0, !is.na(os_time_months))

message(sprintf(
  "[02_GSE20685] N raw=%d | Removed time<=0: %d | NA time: %d | N final=%d",
  n_raw, n_le0, n_na, nrow(out)
))
message(sprintf(
  "[02_GSE20685] OS events: %d (%.1f%%) | Median follow-up: %.1f months",
  sum(out$os_event, na.rm = TRUE),
  100 * mean(out$os_event, na.rm = TRUE),
  median(out$os_time_months, na.rm = TRUE)
))

# =============================================================================
# 4) SAVE
# =============================================================================
dest_dir <- proc_cohort("GSE20685")
dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(dest_dir, "clinical_FINAL.parquet")

arrow::write_parquet(out, out_path)
h    <- sha256_file(out_path)
size <- file.info(out_path)$size / 1024^2
registry_append("GSE20685", "Clinical_FINAL", out_path, h,
                "INTEGRO", SCRIPT_NAME, size)
message(sprintf("[02_GSE20685] Saved: %s | %d samples | %.2f MB | SHA256: %s",
                out_path, nrow(out), size, h))

# =============================================================================
# 5) ENDPOINT MAPPING
# =============================================================================
ep_map <- tibble(
  cohort           = "GSE20685",
  endpoint_primary = "OS",
  col_time_origin  = TIME_COL,
  col_event_origin = EVENT_COL,
  unit_origin      = time_unit_origin,
  unit_final       = "months",
  conversion       = ifelse(time_unit_origin == "days",
                            paste0("/ ", FREEZE$time_unit_divisor), "none"),
  event_1_means    = "death_any_cause",
  n_samples        = nrow(out),
  n_events         = sum(out$os_event, na.rm = TRUE),
  script           = SCRIPT_NAME,
  timestamp        = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)
ep_path <- file.path(PATHS$docs, "endpoint_mapping_templates",
                     "endpoint_mapping_GSE20685.csv")
write_csv(ep_map, ep_path)
message("[02_GSE20685] Endpoint map: ", ep_path)

message("\n[02_GSE20685] Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message("Next step: scripts/03_expression_preprocess_GSE20685.R")
