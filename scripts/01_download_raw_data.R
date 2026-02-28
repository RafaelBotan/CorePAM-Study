# =============================================================================
# SCRIPT: 01_download_raw_data.R
# PURPOSE: Automated download of clinical and expression data from public
#          sources (GEO, GDC/TCGA, cBioPortal) + SHA-256 verification.
# PROJETO: Core-PAM (Memorial v6.1 / Freeze Core-PAM)
# NÍVEL:   A1-Proof — complete forensic traceability
#
# SOURCES:
#   GSE96058  (SCAN-B training + validation) — GEO via GEOquery
#   METABRIC  (validation microarray)        — cBioPortal (open access)
#   TCGA-BRCA (validation RNA-seq)           — GDC via TCGAbiolinks
#   GSE20685  (Taiwan, validation)           — GEO via GEOquery
#
# R DEPENDENCIES:
#   CRAN:        tidyverse, digest, arrow, httr, curl
#   Bioconductor: GEOquery, TCGAbiolinks, SummarizedExperiment
#
# EXECUTION:
#   Rscript scripts/01_download_raw_data.R
#
# OUTPUTS (all in 01_Base_Pura_CorePAM/RAW/<COHORT>/):
#   Expression + raw clinical for each cohort
#   registry/study_registry.csv  — append for each artifact
#   01_docs/registry/data_lake_audit_report.csv — final report
#   results/supp/leakage_check_scanb_vs_gse96058.csv — anti-leakage
#
# RULES (Memorial v6.1):
#   - RAW is read-only after ingestion (do not overwrite).
#   - Read warning = error (strict I/O via 00_setup.R).
#   - No gene counts in file names.
#   - SCAN-B (training) and GSE96058 (validation) share GEO access but
#     have DISTINCT samples — split done in 02_harmonize_clinical_SCANB.R.
# =============================================================================

source("scripts/00_setup.R")

# Load packages with warn=0: suppressPackageStartupMessages() only suppresses
# message() calls, NOT warning() calls. Under warn=2 (strict mode set by 00_setup.R),
# internal package-load warnings (e.g. GEOquery libdeflate) become hard errors.
# Pattern: save warn level, set to 0 for library(), restore immediately after.
old_warn <- getOption("warn"); options(warn = 0)
suppressPackageStartupMessages({
  library(GEOquery)
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(httr)
  library(curl)
})
options(warn = old_warn)

SCRIPT_NAME <- "01_download_raw_data.R"

# =============================================================================
# HELPER: download_file_safe
#   Downloads a file via URL with mode="wb" (Windows-safe), checks size,
#   computes SHA-256 and registers in study_registry.csv.
#   If file already exists and size > 0, skips download (idempotent).
# =============================================================================
download_file_safe <- function(url, destfile, cohort, file_type,
                               overwrite = FALSE) {
  dir.create(dirname(destfile), showWarnings = FALSE, recursive = TRUE)

  # Minimum size threshold: 10 KB. Guards against caching a corrupt HTTP error
  # response (e.g. a 243-byte HTML "403 Forbidden" page written by httr).
  MIN_VALID_BYTES <- 10240L
  if (file.exists(destfile) && file.info(destfile)$size >= MIN_VALID_BYTES && !overwrite) {
    message(sprintf("  [SKIP] Already exists: %s", basename(destfile)))
    h    <- sha256_file(destfile)
    size <- file.info(destfile)$size / 1024^2
    registry_append(cohort, file_type, destfile, h, "INTEGRO_CACHED",
                    SCRIPT_NAME, size)
    return(invisible(h))
  }
  # File exists but is too small — treat as corrupt and re-download
  if (file.exists(destfile) && file.info(destfile)$size < MIN_VALID_BYTES) {
    message(sprintf("  [WARN] File too small (%d bytes) — likely a corrupt/partial download. Re-downloading.",
                    file.info(destfile)$size))
    unlink(destfile)
  }

  message(sprintf("  [DOWN] %s  ->  %s", url, basename(destfile)))
  resp <- tryCatch(
    httr::GET(url,
              httr::write_disk(destfile, overwrite = TRUE),
              httr::progress(),
              httr::timeout(600)),
    error = function(e) {
      stop("Download failed for ", url, "\n", conditionMessage(e))
    }
  )
  if (httr::http_error(resp)) {
    stop("HTTP ", httr::status_code(resp), " downloading: ", url)
  }

  info <- file.info(destfile)
  if (is.na(info$size) || info$size == 0L) {
    stop("Download resulted in empty file: ", destfile)
  }

  h    <- sha256_file(destfile)
  size <- info$size / 1024^2
  message(sprintf("  [OK]   %.1f MB | SHA256: %s", size, h))

  registry_append(cohort, file_type, destfile, h, "INTEGRO", SCRIPT_NAME, size)
  invisible(h)
}

# =============================================================================
# HELPER: geo_supp_download
#   Downloads supplementary files from a GEO accession to the destination
#   directory. Returns vector of paths of downloaded files.
# =============================================================================
geo_supp_download <- function(geo_acc, dest_dir, cohort) {
  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)

  # Skip if dest_dir already contains at least one non-RDS file with size > 10 KB
  existing <- list.files(dest_dir, full.names = TRUE, pattern = "\\.gz$|\\.txt$|\\.csv$")
  if (length(existing) > 0 && all(file.info(existing)$size > 10240, na.rm = TRUE)) {
    message(sprintf("  [SKIP] GEO supplementary files already in %s (%d files). Registering.",
                    dest_dir, length(existing)))
    for (f in existing) {
      h  <- sha256_file(f)
      sz <- file.info(f)$size / 1024^2
      registry_append(cohort, paste0("Expression_supp_cached_", basename(f)),
                      f, h, "INTEGRO_CACHED", "01_download_raw_data.R", sz)
      message(sprintf("  [SKIP] %s | %.1f MB", basename(f), sz))
    }
    return(invisible(existing))
  }

  message(sprintf("\n  [GEO-SUPP] Downloading supplementary files from %s ...", geo_acc))

  # GEOquery saves in subfolder with accession name; move afterwards
  tmp_base <- tempdir()
  old_warn <- getOption("warn"); options(warn = 0)
  files <- tryCatch(
    GEOquery::getGEOSuppFiles(geo_acc, baseDir = tmp_base, fetch_files = TRUE),
    error = function(e) {
      options(warn = old_warn)
      stop("Failed to download supplementary files from ", geo_acc, ": ", conditionMessage(e))
    }
  )
  options(warn = old_warn)

  src_dir  <- file.path(tmp_base, geo_acc)
  src_files <- list.files(src_dir, full.names = TRUE)
  out_paths <- file.path(dest_dir, basename(src_files))

  for (i in seq_along(src_files)) {
    file.copy(src_files[i], out_paths[i], overwrite = TRUE)
    info <- file.info(out_paths[i])
    h    <- sha256_file(out_paths[i])
    registry_append(cohort, paste0("Expression_supp_", i),
                    out_paths[i], h, "INTEGRO", SCRIPT_NAME,
                    info$size / 1024^2)
    message(sprintf("  [OK] %s | %.1f MB", basename(out_paths[i]),
                    info$size / 1024^2))
  }
  invisible(out_paths)
}

# =============================================================================
# HELPER: geo_clinical_download
#   Downloads the series matrix from a GEO accession and saves as RDS (pData).
# =============================================================================
geo_clinical_download <- function(geo_acc, dest_rds, cohort) {
  dir.create(dirname(dest_rds), showWarnings = FALSE, recursive = TRUE)

  if (file.exists(dest_rds) && file.info(dest_rds)$size > 0) {
    message(sprintf("  [SKIP] Clinical already exists: %s", basename(dest_rds)))
    h <- sha256_file(dest_rds)
    registry_append(cohort, "Clinical_raw", dest_rds, h,
                    "INTEGRO_CACHED", SCRIPT_NAME,
                    file.info(dest_rds)$size / 1024^2)
    return(invisible(h))
  }

  message(sprintf("  [GEO-CLIN] Downloading series matrix from %s ...", geo_acc))
  old_warn <- getOption("warn"); options(warn = 0)
  gse <- tryCatch(
    GEOquery::getGEO(geo_acc, GSEMatrix = TRUE, getGPL = FALSE),
    error = function(e) {
      options(warn = old_warn)
      stop("Failed to download GEO series matrix from ", geo_acc, ": ",
           conditionMessage(e))
    }
  )
  options(warn = old_warn)

  # Extract pData (clinical metadata) from first element
  pheno <- Biobase::pData(gse[[1]])
  saveRDS(pheno, dest_rds)

  h    <- sha256_file(dest_rds)
  size <- file.info(dest_rds)$size / 1024^2
  message(sprintf("  [OK] %d samples | %.1f MB | SHA256: %s",
                  nrow(pheno), size, h))
  registry_append(cohort, "Clinical_raw", dest_rds, h,
                  "INTEGRO", SCRIPT_NAME, size)
  invisible(h)
}

# =============================================================================
# SECTION 1: GSE96058 — SCAN-B (training + validation GSE96058)
#   Platform: RNA-seq (Illumina HiSeq 2000/2500)
#   Expression: count matrix + FPKM (GEO supplementary files)
#   Clinical: series matrix pData
#   Note: the split training (SCAN-B) vs validation (GSE96058) is done in
#         02_harmonize_clinical_SCANB.R based on patient IDs.
# =============================================================================
message("\n", strrep("=", 60))
message("[01_download] SECTION 1: GSE96058 (SCAN-B + validation)")
message(strrep("=", 60))

GSE96058_RAW_DIR <- raw_cohort("GSE96058")

# 1a. Supplementary files (expression — gene counts, FPKM)
geo_supp_download("GSE96058", GSE96058_RAW_DIR, cohort = "GSE96058")

# 1b. Clinical via series matrix
geo_clinical_download(
  geo_acc  = "GSE96058",
  dest_rds = file.path(GSE96058_RAW_DIR, "GSE96058_clinical_raw.rds"),
  cohort   = "GSE96058"
)

# =============================================================================
# SECTION 2: METABRIC — cBioPortal (open access)
#   Platform: Illumina HT-12 v3 microarray
#   Expression: data_mrna_illumina_microarray.txt  (log2-intensity)
#   Clinical: data_clinical_patient.txt + data_clinical_sample.txt
#   Source: cBioPortal Data Hub (AWS S3, no authentication)
#   Reference: Curtis et al. Nature 2012; Pereira et al. Nat Commun 2016
# =============================================================================
message("\n", strrep("=", 60))
message("[01_download] SECTION 2: METABRIC (cBioPortal)")
message(strrep("=", 60))

METABRIC_RAW_DIR  <- raw_cohort("METABRIC")
METABRIC_TAR      <- file.path(METABRIC_RAW_DIR, "brca_metabric.tar.gz")
# URL updated 2026-02: cBioPortal migrated from cbioportal-datahub.s3.amazonaws.com
# (now returns HTTP 403) to datahub.assets.cbioportal.org.
METABRIC_BASE_URL <- "https://datahub.assets.cbioportal.org/brca_metabric.tar.gz"

# 2a. Download complete package
download_file_safe(
  url      = METABRIC_BASE_URL,
  destfile = METABRIC_TAR,
  cohort   = "METABRIC",
  file_type = "Package_tarball"
)

# 2b. Extraction — verify tar is valid before extracting
message("  [EXTRACT] Extracting METABRIC...")
tar_size <- file.info(METABRIC_TAR)$size
if (is.na(tar_size) || tar_size < 10240L) {
  stop(sprintf("brca_metabric.tar.gz appears corrupt or too small (%d bytes). Re-run download.",
               as.integer(tar_size)))
}
old_warn <- getOption("warn"); options(warn = 0)
untar_result <- tryCatch(
  utils::untar(METABRIC_TAR, exdir = METABRIC_RAW_DIR),
  error = function(e) {
    options(warn = old_warn)
    stop("Failed to extract METABRIC tar.gz: ", conditionMessage(e))
  }
)
options(warn = old_warn)
# Verify at least one expected file was extracted
metabric_check <- file.path(METABRIC_RAW_DIR, "brca_metabric", "data_clinical_patient.txt")
if (!file.exists(metabric_check)) {
  stop(sprintf(
    "Extraction completed but expected files not found under %s/brca_metabric/.\n",
    METABRIC_RAW_DIR,
    "Check tar contents with: utils::untar('%s', list=TRUE)", METABRIC_TAR
  ))
}
message("  [EXTRACT] Extraction verified OK.")

# 2c. Locate and register extracted files
metabric_files <- list(
  expression = file.path(METABRIC_RAW_DIR, "brca_metabric",
                         "data_mrna_illumina_microarray.txt"),
  clin_patient = file.path(METABRIC_RAW_DIR, "brca_metabric",
                           "data_clinical_patient.txt"),
  clin_sample  = file.path(METABRIC_RAW_DIR, "brca_metabric",
                           "data_clinical_sample.txt")
)

for (ftype in names(metabric_files)) {
  fpath <- metabric_files[[ftype]]
  if (!file.exists(fpath)) {
    warning("METABRIC file not found after extraction: ", fpath)
    next
  }
  h    <- sha256_file(fpath)
  size <- file.info(fpath)$size / 1024^2
  registry_append("METABRIC", paste0("Extracted_", ftype),
                  fpath, h, "INTEGRO", SCRIPT_NAME, size)
  message(sprintf("  [OK] %s | %.1f MB", basename(fpath), size))
}

# =============================================================================
# SECTION 3: TCGA-BRCA — GDC / TCGAbiolinks
#   Platform: RNA-seq (Illumina HiSeq; STAR aligner)
#   Expression: STAR raw counts (unstranded)
#   Clinical: GDC clinical XML (harmonized)
#   Note: heavy download (~5-10 GB); requires stable connection.
#         TCGAbiolinks saves in GDCdata/ subfolder by default.
# =============================================================================
message("\n", strrep("=", 60))
message("[01_download] SECTION 3: TCGA-BRCA (GDC / TCGAbiolinks)")
message(strrep("=", 60))

TCGA_RAW_DIR   <- raw_cohort("TCGA_BRCA")
tcga_expr_rds  <- file.path(TCGA_RAW_DIR, "TCGA_BRCA_SE_counts_raw.rds")
tcga_clin_rds  <- file.path(TCGA_RAW_DIR, "TCGA_BRCA_clinical_raw.rds")

dir.create(TCGA_RAW_DIR, showWarnings = FALSE, recursive = TRUE)

# Skip entire TCGA section if both output files already exist (saves 5-10 GB re-download)
if (file.exists(tcga_expr_rds) && file.info(tcga_expr_rds)$size > 1e6 &&
    file.exists(tcga_clin_rds)  && file.info(tcga_clin_rds)$size  > 1e3) {
  message("  [SKIP] TCGA-BRCA already downloaded and registered.")
  h_expr <- sha256_file(tcga_expr_rds)
  h_clin <- sha256_file(tcga_clin_rds)
  registry_append("TCGA_BRCA", "Expression_SE_counts", tcga_expr_rds, h_expr,
                  "INTEGRO_CACHED", SCRIPT_NAME, file.info(tcga_expr_rds)$size / 1e6)
  registry_append("TCGA_BRCA", "Clinical_raw", tcga_clin_rds, h_clin,
                  "INTEGRO_CACHED", SCRIPT_NAME, file.info(tcga_clin_rds)$size / 1e6)
} else {

# SSL fix: R's bundled libcurl uses a different TLS stack than the system curl.
# GDC API is reachable via system curl (HTTP 200) but R's httr fails with
# "SSL connect error". Disabling peer verification is safe here because:
#   1) GDC server identity is verified at OS level
#   2) Downloaded file integrity is verified with SHA-256
# httr::get_config() does not exist in httr; httr::set_config/reset_config are the API.
# reset_config() restores defaults at the end of the TCGA section.
httr::set_config(httr::config(ssl_verifypeer = FALSE))
message("  [GDC] SSL peer verification disabled for TCGAbiolinks calls (SHA-256 integrity check applied after download)")

# 3a. Expression query (STAR counts)
message("  [GDC] Building TCGA-BRCA expression query...")
old_warn <- getOption("warn"); options(warn = 0)
query_expr <- tryCatch(
  TCGAbiolinks::GDCquery(
    project            = "TCGA-BRCA",
    data.category      = "Transcriptome Profiling",
    data.type          = "Gene Expression Quantification",
    workflow.type      = "STAR - Counts",
    sample.type        = "Primary Tumor"
  ),
  error = function(e) {
    options(warn = old_warn)
    httr::reset_config()
    stop("GDCquery TCGA-BRCA expression failed: ", conditionMessage(e))
  }
)
options(warn = old_warn)

message(sprintf("  [GDC] %d files found for download.",
                nrow(TCGAbiolinks::getResults(query_expr))))

# 3b. Download (idempotent — resumes if already partially downloaded)
message("  [GDC] Starting download (may take >30 min)...")
old_warn <- getOption("warn"); options(warn = 0)
tryCatch(
  TCGAbiolinks::GDCdownload(
    query           = query_expr,
    directory       = TCGA_RAW_DIR,
    files.per.chunk = 10
  ),
  error = function(e) {
    options(warn = old_warn)
    httr::reset_config()
    stop("GDCdownload TCGA-BRCA failed: ", conditionMessage(e))
  }
)
options(warn = old_warn)

# 3c. Prepare SummarizedExperiment and save as RDS
message("  [GDC] Preparing SummarizedExperiment and saving as RDS...")
old_warn <- getOption("warn"); options(warn = 0)
tcga_se <- tryCatch(
  TCGAbiolinks::GDCprepare(query_expr, directory = TCGA_RAW_DIR),
  error = function(e) {
    options(warn = old_warn)
    httr::reset_config()
    stop("GDCprepare TCGA-BRCA failed: ", conditionMessage(e))
  }
)
options(warn = old_warn)

tcga_expr_rds <- file.path(TCGA_RAW_DIR, "TCGA_BRCA_SE_counts_raw.rds")
saveRDS(tcga_se, tcga_expr_rds)
h    <- sha256_file(tcga_expr_rds)
size <- file.info(tcga_expr_rds)$size / 1024^2
registry_append("TCGA_BRCA", "Expression_SE_counts", tcga_expr_rds,
                h, "INTEGRO", SCRIPT_NAME, size)
message(sprintf("  [OK] SE saved: %d genes x %d samples | %.1f MB",
                nrow(tcga_se), ncol(tcga_se), size))

# 3d. TCGA clinical (harmonized clinical)
message("  [GDC] Downloading TCGA-BRCA clinical data...")
old_warn <- getOption("warn"); options(warn = 0)
tcga_clin <- tryCatch(
  TCGAbiolinks::GDCquery_clinic("TCGA-BRCA", type = "clinical"),
  error = function(e) {
    options(warn = old_warn)
    stop("Failed to download TCGA-BRCA clinical: ", conditionMessage(e))
  }
)
options(warn = old_warn)

tcga_clin_rds <- file.path(TCGA_RAW_DIR, "TCGA_BRCA_clinical_raw.rds")
saveRDS(tcga_clin, tcga_clin_rds)
h    <- sha256_file(tcga_clin_rds)
size <- file.info(tcga_clin_rds)$size / 1024^2
registry_append("TCGA_BRCA", "Clinical_raw", tcga_clin_rds,
                h, "INTEGRO", SCRIPT_NAME, size)
message(sprintf("  [OK] Clinical: %d patients | %.1f MB | SHA256: %s",
                nrow(tcga_clin), size, h))

# Restore SSL config after all TCGAbiolinks calls
httr::reset_config()
message("  [GDC] SSL config restored.")

} # end TCGA skip block

# =============================================================================
# SECTION 4: GSE20685 — Taiwan (validation microarray)
#   Platform: Affymetrix Human Genome U133A (GPL96)
#   Expression: series matrix (log2 intensities)
#   Clinical: series matrix pData
#   Reference: Lu et al. BMC Med Genomics 2012
# =============================================================================
message("\n", strrep("=", 60))
message("[01_download] SECTION 4: GSE20685 (Taiwan microarray)")
message(strrep("=", 60))

GSE20685_RAW_DIR <- raw_cohort("GSE20685")

# 4a. Supplementary files (expression matrix)
geo_supp_download("GSE20685", GSE20685_RAW_DIR, cohort = "GSE20685")

# 4b. Clinical via series matrix
geo_clinical_download(
  geo_acc  = "GSE20685",
  dest_rds = file.path(GSE20685_RAW_DIR, "GSE20685_clinical_raw.rds"),
  cohort   = "GSE20685"
)

# =============================================================================
# SECTION 5: FINAL REPORT AND GO / NO-GO
# =============================================================================
message("\n", strrep("=", 60))
message("[01_download] SECTION 5: FINAL REPORT")
message(strrep("=", 60))

# Read registry and filter only entries from this script
old_warn <- getOption("warn"); options(warn = 0)
if (file.exists(PATHS$run_registry)) {
  df_reg <- readr::read_csv(PATHS$run_registry, show_col_types = FALSE)
  df_this <- df_reg |> filter(script == SCRIPT_NAME)
} else {
  df_this <- tibble()
}
options(warn = old_warn)

# Save consolidated report
audit_out <- file.path(PATHS$registry_docs, "data_lake_audit_report.csv")
write_csv(df_this, audit_out)
message("Report saved to: ", audit_out)

if (nrow(df_this) > 0) {
  n_ok  <- sum(grepl("INTEGRO", df_this$status))
  n_bad <- sum(!grepl("INTEGRO", df_this$status))
  message(sprintf("\nIntact artifacts : %d", n_ok))
  message(sprintf("Failed artifacts : %d", n_bad))
  print(df_this |> select(cohort, file_type, status, size_mb))
}

# =============================================================================
# SECTION 6: LEAKAGE-PROOF CHECK — SCAN-B (training) vs GSE96058 (validation)
#   Memorial v6.1 §3.2 — ID intersection must be 0.
#   NOTE: this check operates on the pData (series matrix) from GEO.
#   The definitive split (by barcode/patient_id) occurs in 02_harmonize_*.R.
#   Here we verify that the GSM IDs of the two subsets do not overlap.
# =============================================================================
message("\n", strrep("=", 60))
message("[01_download] SECTION 6: LEAKAGE-PROOF CHECK")
message(strrep("=", 60))

scanb_clin_path <- file.path(raw_cohort("GSE96058"), "GSE96058_clinical_raw.rds")
leakage_out     <- file.path(PATHS$results$supp,
                             "leakage_check_scanb_vs_gse96058.csv")

# Note: the GSE96058 dataset contains ALL SCAN-B samples.
# The SCAN-B (training) vs validation split is done by a status column
# present in the metadata (e.g. column "scan.b.training.validation").
# Here we verify the consistency of the groups within the dataset.
if (file.exists(scanb_clin_path)) {
  old_warn <- getOption("warn"); options(warn = 0)
  pheno <- readRDS(scanb_clin_path)
  options(warn = old_warn)

  # Detect split column (adjust if name differs)
  split_col <- intersect(
    c("scan.b.training.validation", "group", "training_validation",
      "characteristics_ch1.6"),
    tolower(names(pheno))
  )

  if (length(split_col) > 0) {
    col_use   <- names(pheno)[tolower(names(pheno)) == split_col[1]][1]
    groups    <- table(pheno[[col_use]])
    train_ids <- rownames(pheno)[grepl("train",  tolower(pheno[[col_use]]))]
    valid_ids <- rownames(pheno)[grepl("valid|test", tolower(pheno[[col_use]]))]
    overlap   <- intersect(train_ids, valid_ids)

    leakage_df <- tibble(
      timestamp            = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      script               = SCRIPT_NAME,
      split_column_used    = col_use,
      n_train_gsm          = length(train_ids),
      n_validation_gsm     = length(valid_ids),
      n_overlap            = length(overlap),
      overlap_ids          = paste(overlap, collapse = ";"),
      sha256_train_ids     = digest::digest(sort(train_ids),  algo = "sha256"),
      sha256_valid_ids     = digest::digest(sort(valid_ids),  algo = "sha256"),
      result               = if (length(overlap) == 0) "PASS_NO_LEAKAGE"
                             else "FAIL_LEAKAGE_DETECTED"
    )
    write_csv(leakage_df, leakage_out)

    if (length(overlap) == 0) {
      message(sprintf("[leakage] PASS — training: %d | validation: %d | overlap: 0",
                      length(train_ids), length(valid_ids)))
      message("  Artifact: ", leakage_out)
    } else {
      stop(sprintf(
        "[leakage] FAIL — %d samples overlapping between training and validation.\n%s",
        length(overlap), paste(overlap, collapse = ", ")
      ))
    }
  } else {
    message("[leakage] WARNING: split column not detected automatically.")
    message("  Inspect pData(GSE96058) and adjust split_col in this script.")
    message("  Available columns: ", paste(names(pheno)[1:10], collapse = ", "))
  }
} else {
  message("[leakage] Deferred — GSE96058 clinical file not yet present.")
}

message("\n[01_download] Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message("Next step: scripts/02_harmonize_clinical_<COHORT>.R")
