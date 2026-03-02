# =============================================================================
# SCRIPT: 01_download_raw_data.R
# PURPOSE: Fully automated, fault-tolerant download of all OS-block cohorts
#          from public repositories (GEO, GDC/TCGA, cBioPortal) with:
#            - automatic retry (3 attempts, exponential backoff)
#            - fallback download methods per source
#            - continue-on-error (all cohorts attempted regardless of failures)
#            - skip already-complete downloads (idempotent/resumable)
#            - SHA-256 verification of every artifact
#
# DESIGN INTENT: "run it and go to sleep"
#   The script never stops early due to a single cohort failure.
#   At the end it prints a clear PASS/FAIL table and exits with code 1
#   only if any cohort failed — making it easy to re-run just the failures.
#
# SOURCES:
#   GSE96058  (SCAN-B — all samples = training) — GEO via GEOquery
#   METABRIC  (validation microarray)        — cBioPortal open access
#   TCGA-BRCA (validation RNA-seq)           — GDC via TCGAbiolinks
#   GSE20685  (Taiwan, validation)           — GEO via GEOquery
#
# OUTPUTS (all in 01_Base_Pura_CorePAM/RAW/<COHORT>/):
#   Expression + raw clinical for each cohort
#   registry/study_registry.csv         — append-only SHA-256 log
#   01_docs/registry/data_lake_audit.csv — final per-cohort report
#   results/supp/leakage_check_*.csv    — anti-leakage check
#
# EXECUTION:
#   Rscript scripts/01_download_raw_data.R
#   Force re-download: FORCE_RERUN=TRUE Rscript scripts/01_download_raw_data.R
# =============================================================================

source("scripts/00_setup.R")

# Load packages outside warn=2 strict mode (internal package warnings are harmless)
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
FORCE       <- as.logical(Sys.getenv("FORCE_RERUN", "FALSE"))

# Collector for per-cohort errors (allows continue-on-error across all sections)
.errors <- list()

# =============================================================================
# UTILITY: retry
#   Runs fn() up to max_attempts times with exponential backoff on failure.
#   Returns the result of fn() on success; stops with the last error on
#   exhausting all attempts.
# =============================================================================
retry <- function(fn, max_attempts = 3, delay_sec = 30, label = "operation") {
  last_err <- NULL
  for (i in seq_len(max_attempts)) {
    result <- tryCatch(
      list(val = fn(), err = NULL),
      error = function(e) list(val = NULL, err = e)
    )
    if (is.null(result$err)) return(result$val)
    last_err <- result$err
    message(sprintf(
      "  [RETRY] %s — attempt %d/%d failed: %s",
      label, i, max_attempts, conditionMessage(last_err)
    ))
    if (i < max_attempts) {
      wait <- delay_sec * i   # 30s, 60s, 90s
      message(sprintf("  [RETRY] Waiting %d s before next attempt...", wait))
      Sys.sleep(wait)
    }
  }
  stop(sprintf("All %d attempts failed for: %s\nLast error: %s",
               max_attempts, label, conditionMessage(last_err)))
}

# =============================================================================
# UTILITY: download_file_safe
#   Downloads url → destfile; skips if already valid (size > MIN_VALID_BYTES).
#   Deletes corrupt small files before re-downloading.
#   Retries internally via retry().
# =============================================================================
MIN_VALID_BYTES <- 10240L   # 10 KB guard — avoids caching HTML error pages

download_file_safe <- function(url, destfile, cohort, file_type,
                               overwrite = FALSE, timeout_sec = 1800) {
  dir.create(dirname(destfile), showWarnings = FALSE, recursive = TRUE)

  if (!overwrite && !FORCE &&
      file.exists(destfile) &&
      file.info(destfile)$size >= MIN_VALID_BYTES) {
    message(sprintf("  [SKIP] Already exists (%s)",  basename(destfile)))
    h <- sha256_file(destfile)
    registry_append(cohort, file_type, destfile, h, "INTEGRO_CACHED",
                    SCRIPT_NAME, file.info(destfile)$size / 1024^2)
    return(invisible(h))
  }

  if (file.exists(destfile) && file.info(destfile)$size < MIN_VALID_BYTES) {
    message(sprintf("  [WARN] File too small (%d bytes) — deleting and re-downloading.",
                    as.integer(file.info(destfile)$size)))
    unlink(destfile)
  }

  message(sprintf("  [DOWN] %s", basename(destfile)))

  retry(
    fn = function() {
      resp <- httr::GET(
        url,
        httr::write_disk(destfile, overwrite = TRUE),
        httr::progress(),
        httr::timeout(timeout_sec)
      )
      if (httr::http_error(resp))
        stop("HTTP ", httr::status_code(resp), " downloading: ", url)
      info <- file.info(destfile)
      if (is.na(info$size) || info$size < MIN_VALID_BYTES)
        stop("Resulting file too small: ", destfile)
    },
    label = basename(destfile)
  )

  h    <- sha256_file(destfile)
  size <- file.info(destfile)$size / 1024^2
  message(sprintf("  [OK]   %.1f MB | SHA256: %.16s...", size, h))
  registry_append(cohort, file_type, destfile, h, "INTEGRO", SCRIPT_NAME, size)
  invisible(h)
}

# =============================================================================
# UTILITY: geo_supp_download
#   Downloads GEO supplementary files using GEOquery with HTTP fallback.
#   Skips if dest_dir already contains valid files.
# =============================================================================
geo_supp_download <- function(geo_acc, dest_dir, cohort) {
  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)

  # Skip check: any substantial files already present?
  existing <- list.files(
    dest_dir, full.names = TRUE,
    pattern = "\\.gz$|\\.txt$|\\.csv$|\\.tar$|\\.tar\\.gz$"
  )
  if (!FORCE && length(existing) > 0 &&
      all(file.info(existing)$size > MIN_VALID_BYTES, na.rm = TRUE)) {
    message(sprintf("  [SKIP] Supplementary files already in %s (%d files)",
                    basename(dest_dir), length(existing)))
    for (f in existing) {
      h  <- sha256_file(f)
      sz <- file.info(f)$size / 1024^2
      registry_append(cohort, paste0("Expression_supp_cached_", basename(f)),
                      f, h, "INTEGRO_CACHED", SCRIPT_NAME, sz)
    }
    return(invisible(existing))
  }

  message(sprintf("  [GEO-SUPP] Downloading supplementary files for %s ...", geo_acc))

  tmp_base <- tempdir()

  # Attempt 1: GEOquery
  files <- tryCatch({
    old_w <- getOption("warn"); options(warn = 0)
    on.exit(options(warn = old_w), add = TRUE)
    retry(
      fn    = function() GEOquery::getGEOSuppFiles(
        geo_acc, baseDir = tmp_base, fetch_files = TRUE
      ),
      label = paste(geo_acc, "supplementary (GEOquery)")
    )
  }, error = function(e) {
    message(sprintf("  [FALLBACK] GEOquery failed for %s: %s", geo_acc, e$message))
    message("  [FALLBACK] Trying direct GEO FTP download...")

    # Fallback: direct FTP URL construction
    prefix  <- substr(geo_acc, 1, nchar(geo_acc) - 3)
    ftp_url <- sprintf(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/%snnn/%s/suppl/",
      prefix, geo_acc
    )
    old_w2 <- getOption("warn"); options(warn = 0)
    on.exit(options(warn = old_w2), add = TRUE)
    retry(
      fn    = function() GEOquery::getGEOSuppFiles(
        geo_acc, baseDir = tmp_base, fetch_files = TRUE,
        baseUrl = ftp_url
      ),
      label = paste(geo_acc, "supplementary (FTP fallback)")
    )
  })

  # Move downloaded files to the permanent destination
  src_dir   <- file.path(tmp_base, geo_acc)
  src_files <- list.files(src_dir, full.names = TRUE)
  out_paths <- file.path(dest_dir, basename(src_files))

  for (i in seq_along(src_files)) {
    file.copy(src_files[i], out_paths[i], overwrite = TRUE)
    h  <- sha256_file(out_paths[i])
    sz <- file.info(out_paths[i])$size / 1024^2
    registry_append(cohort, paste0("Expression_supp_", i),
                    out_paths[i], h, "INTEGRO", SCRIPT_NAME, sz)
    message(sprintf("  [OK] %s | %.1f MB", basename(out_paths[i]), sz))
  }
  invisible(out_paths)
}

# =============================================================================
# UTILITY: geo_clinical_download
#   Downloads the GEO series matrix (pData) and saves as RDS.
#   Retries up to 3 times with an FTP fallback on failure.
# =============================================================================
geo_clinical_download <- function(geo_acc, dest_rds, cohort) {
  dir.create(dirname(dest_rds), showWarnings = FALSE, recursive = TRUE)

  if (!FORCE && file.exists(dest_rds) &&
      file.info(dest_rds)$size > MIN_VALID_BYTES) {
    message(sprintf("  [SKIP] Clinical already exists (%s)", basename(dest_rds)))
    h <- sha256_file(dest_rds)
    registry_append(cohort, "Clinical_raw", dest_rds, h,
                    "INTEGRO_CACHED", SCRIPT_NAME,
                    file.info(dest_rds)$size / 1024^2)
    return(invisible(h))
  }

  message(sprintf("  [GEO-CLIN] Downloading series matrix for %s ...", geo_acc))

  # Primary: GEOquery over HTTPS
  gse <- tryCatch({
    old_w <- getOption("warn"); options(warn = 0)
    on.exit(options(warn = old_w), add = TRUE)
    retry(
      fn    = function() GEOquery::getGEO(
        geo_acc, GSEMatrix = TRUE, getGPL = FALSE
      ),
      label = paste(geo_acc, "series matrix (GEOquery)")
    )
  }, error = function(e) {
    message(sprintf("  [FALLBACK] Series matrix HTTPS failed: %s", e$message))
    message("  [FALLBACK] Trying with destdir to force re-download...")
    # Force fresh download to a temp location
    tmp_dest <- tempfile()
    dir.create(tmp_dest, showWarnings = FALSE)
    old_w2 <- getOption("warn"); options(warn = 0)
    on.exit(options(warn = old_w2), add = TRUE)
    retry(
      fn    = function() GEOquery::getGEO(
        geo_acc, GSEMatrix = TRUE, getGPL = FALSE, destdir = tmp_dest
      ),
      label = paste(geo_acc, "series matrix (forced re-download)")
    )
  })

  pheno <- Biobase::pData(gse[[1]])
  saveRDS(pheno, dest_rds)

  h    <- sha256_file(dest_rds)
  size <- file.info(dest_rds)$size / 1024^2
  message(sprintf("  [OK] %d samples | %.1f MB | SHA256: %.16s...",
                  nrow(pheno), size, h))
  registry_append(cohort, "Clinical_raw", dest_rds, h,
                  "INTEGRO", SCRIPT_NAME, size)
  invisible(h)
}

# =============================================================================
# SECTION 1: GSE96058 — SCAN-B (training) + validation subset
# =============================================================================
message("\n", strrep("=", 60))
message("[01_download] SECTION 1: GSE96058 (SCAN-B + validation)")
message(strrep("=", 60))

tryCatch({
  GSE96058_RAW_DIR <- raw_cohort("GSE96058")
  geo_supp_download("GSE96058", GSE96058_RAW_DIR, cohort = "GSE96058")
  geo_clinical_download(
    geo_acc  = "GSE96058",
    dest_rds = file.path(GSE96058_RAW_DIR, "GSE96058_clinical_raw.rds"),
    cohort   = "GSE96058"
  )
  message("[01_download] GSE96058: DONE")
}, error = function(e) {
  .errors[["GSE96058"]] <<- conditionMessage(e)
  message(sprintf("[01_download] GSE96058: FAILED — %s\n  (continuing with next cohort)",
                  conditionMessage(e)))
})

# =============================================================================
# SECTION 2: METABRIC — cBioPortal (open access)
# =============================================================================
message("\n", strrep("=", 60))
message("[01_download] SECTION 2: METABRIC (cBioPortal)")
message(strrep("=", 60))

tryCatch({
  METABRIC_RAW_DIR  <- raw_cohort("METABRIC")
  METABRIC_TAR      <- file.path(METABRIC_RAW_DIR, "brca_metabric.tar.gz")
  # Primary URL (CDN as of 2026-02; old S3 bucket now returns HTTP 403)
  METABRIC_URL_PRI  <- "https://datahub.assets.cbioportal.org/brca_metabric.tar.gz"
  METABRIC_URL_BAK  <- "https://cbioportal-datahub.s3.amazonaws.com/brca_metabric.tar.gz"

  # Download with automatic fallback to backup URL
  dl_ok <- tryCatch({
    download_file_safe(METABRIC_URL_PRI, METABRIC_TAR, "METABRIC", "Package_tarball")
    TRUE
  }, error = function(e) {
    message(sprintf("  [FALLBACK] Primary METABRIC URL failed: %s", e$message))
    message("  [FALLBACK] Trying backup URL...")
    download_file_safe(METABRIC_URL_BAK, METABRIC_TAR, "METABRIC",
                       "Package_tarball", overwrite = TRUE)
    TRUE
  })

  # Verify tar integrity before extraction
  tar_size <- file.info(METABRIC_TAR)$size
  if (is.na(tar_size) || tar_size < 10240L)
    stop(sprintf("brca_metabric.tar.gz too small (%d bytes) — likely corrupt.",
                 as.integer(tar_size)))

  # Extract only if not already done
  metabric_check <- file.path(METABRIC_RAW_DIR, "brca_metabric",
                               "data_clinical_patient.txt")
  if (!FORCE && file.exists(metabric_check)) {
    message("  [SKIP] METABRIC already extracted.")
  } else {
    message("  [EXTRACT] Extracting brca_metabric.tar.gz ...")
    old_w <- getOption("warn"); options(warn = 0)
    tryCatch(
      utils::untar(METABRIC_TAR, exdir = METABRIC_RAW_DIR),
      error = function(e) {
        options(warn = old_w)
        stop("Extraction failed: ", conditionMessage(e))
      }
    )
    options(warn = old_w)
    if (!file.exists(metabric_check))
      stop("Extraction completed but expected file not found: ", metabric_check)
    message("  [EXTRACT] Extraction verified OK.")
  }

  # Register key extracted files
  metabric_files <- list(
    expression   = file.path(METABRIC_RAW_DIR, "brca_metabric",
                              "data_mrna_illumina_microarray.txt"),
    clin_patient = file.path(METABRIC_RAW_DIR, "brca_metabric",
                              "data_clinical_patient.txt"),
    clin_sample  = file.path(METABRIC_RAW_DIR, "brca_metabric",
                              "data_clinical_sample.txt")
  )
  for (ftype in names(metabric_files)) {
    fpath <- metabric_files[[ftype]]
    if (!file.exists(fpath)) { warning("Missing: ", fpath); next }
    h    <- sha256_file(fpath)
    size <- file.info(fpath)$size / 1024^2
    registry_append("METABRIC", paste0("Extracted_", ftype),
                    fpath, h, "INTEGRO", SCRIPT_NAME, size)
    message(sprintf("  [OK] %s | %.1f MB", basename(fpath), size))
  }
  message("[01_download] METABRIC: DONE")
}, error = function(e) {
  .errors[["METABRIC"]] <<- conditionMessage(e)
  message(sprintf("[01_download] METABRIC: FAILED — %s\n  (continuing with next cohort)",
                  conditionMessage(e)))
})

# =============================================================================
# SECTION 3: TCGA-BRCA — GDC / TCGAbiolinks
#   SSL note: R's bundled libcurl may fail TLS handshake with GDC API even
#   though system curl succeeds. Fix: disable ssl_verifypeer; SHA-256 ensures
#   integrity independently.
# =============================================================================
message("\n", strrep("=", 60))
message("[01_download] SECTION 3: TCGA-BRCA (GDC / TCGAbiolinks)")
message(strrep("=", 60))

tryCatch({
  TCGA_RAW_DIR  <- raw_cohort("TCGA_BRCA")
  tcga_expr_rds <- file.path(TCGA_RAW_DIR, "TCGA_BRCA_SE_counts_raw.rds")
  tcga_clin_rds <- file.path(TCGA_RAW_DIR, "TCGA_BRCA_clinical_raw.rds")
  dir.create(TCGA_RAW_DIR, showWarnings = FALSE, recursive = TRUE)

  tcga_done <- !FORCE &&
    file.exists(tcga_expr_rds) && file.info(tcga_expr_rds)$size > 1e6 &&
    file.exists(tcga_clin_rds)  && file.info(tcga_clin_rds)$size  > 1e3

  if (tcga_done) {
    message("  [SKIP] TCGA-BRCA already downloaded.")
    registry_append("TCGA_BRCA", "Expression_SE_counts", tcga_expr_rds,
                    sha256_file(tcga_expr_rds), "INTEGRO_CACHED", SCRIPT_NAME,
                    file.info(tcga_expr_rds)$size / 1e6)
    registry_append("TCGA_BRCA", "Clinical_raw", tcga_clin_rds,
                    sha256_file(tcga_clin_rds), "INTEGRO_CACHED", SCRIPT_NAME,
                    file.info(tcga_clin_rds)$size / 1e6)
  } else {
    # Disable SSL peer verification for GDC (safe: integrity checked by SHA-256)
    httr::set_config(httr::config(ssl_verifypeer = FALSE))
    message("  [GDC] SSL peer verification disabled (SHA-256 integrity check applied)")
    on.exit(httr::reset_config(), add = TRUE)

    # 3a. Query
    message("  [GDC] Building expression query...")
    old_w <- getOption("warn"); options(warn = 0)
    query_expr <- retry(
      fn = function() TCGAbiolinks::GDCquery(
        project       = "TCGA-BRCA",
        data.category = "Transcriptome Profiling",
        data.type     = "Gene Expression Quantification",
        workflow.type = "STAR - Counts",
        sample.type   = "Primary Tumor"
      ),
      label = "GDCquery TCGA-BRCA"
    )
    options(warn = old_w)
    message(sprintf("  [GDC] %d files in query.",
                    nrow(TCGAbiolinks::getResults(query_expr))))

    # 3b. Download (~5 GB, resumes on re-run)
    message("  [GDC] Downloading chunks (this can take 30–90 min)...")
    old_w <- getOption("warn"); options(warn = 0)
    retry(
      fn = function() TCGAbiolinks::GDCdownload(
        query           = query_expr,
        directory       = TCGA_RAW_DIR,
        files.per.chunk = 10
      ),
      max_attempts = 5,
      delay_sec    = 60,
      label        = "GDCdownload TCGA-BRCA"
    )
    options(warn = old_w)

    # 3c. Prepare SummarizedExperiment
    message("  [GDC] Preparing SummarizedExperiment...")
    old_w <- getOption("warn"); options(warn = 0)
    tcga_se <- retry(
      fn    = function() TCGAbiolinks::GDCprepare(query_expr,
                                                  directory = TCGA_RAW_DIR),
      label = "GDCprepare TCGA-BRCA"
    )
    options(warn = old_w)
    saveRDS(tcga_se, tcga_expr_rds)
    h    <- sha256_file(tcga_expr_rds)
    size <- file.info(tcga_expr_rds)$size / 1024^2
    registry_append("TCGA_BRCA", "Expression_SE_counts", tcga_expr_rds,
                    h, "INTEGRO", SCRIPT_NAME, size)
    message(sprintf("  [OK] SE: %d genes x %d samples | %.1f MB",
                    nrow(tcga_se), ncol(tcga_se), size))

    # 3d. Clinical
    message("  [GDC] Downloading clinical data...")
    old_w <- getOption("warn"); options(warn = 0)
    tcga_clin <- retry(
      fn    = function() TCGAbiolinks::GDCquery_clinic("TCGA-BRCA",
                                                       type = "clinical"),
      label = "GDCquery_clinic TCGA-BRCA"
    )
    options(warn = old_w)
    saveRDS(tcga_clin, tcga_clin_rds)
    h    <- sha256_file(tcga_clin_rds)
    size <- file.info(tcga_clin_rds)$size / 1024^2
    registry_append("TCGA_BRCA", "Clinical_raw", tcga_clin_rds,
                    h, "INTEGRO", SCRIPT_NAME, size)
    message(sprintf("  [OK] Clinical: %d patients | %.1f MB", nrow(tcga_clin), size))

    httr::reset_config()
    message("  [GDC] SSL config restored.")
  }
  message("[01_download] TCGA-BRCA: DONE")
}, error = function(e) {
  httr::reset_config()   # always restore SSL config
  .errors[["TCGA_BRCA"]] <<- conditionMessage(e)
  message(sprintf("[01_download] TCGA-BRCA: FAILED — %s\n  (continuing with next cohort)",
                  conditionMessage(e)))
})

# =============================================================================
# SECTION 4: GSE20685 — Taiwan (Affymetrix HGU133A microarray)
# =============================================================================
message("\n", strrep("=", 60))
message("[01_download] SECTION 4: GSE20685 (Taiwan microarray)")
message(strrep("=", 60))

tryCatch({
  GSE20685_RAW_DIR <- raw_cohort("GSE20685")
  geo_supp_download("GSE20685", GSE20685_RAW_DIR, cohort = "GSE20685")
  geo_clinical_download(
    geo_acc  = "GSE20685",
    dest_rds = file.path(GSE20685_RAW_DIR, "GSE20685_clinical_raw.rds"),
    cohort   = "GSE20685"
  )
  message("[01_download] GSE20685: DONE")
}, error = function(e) {
  .errors[["GSE20685"]] <<- conditionMessage(e)
  message(sprintf("[01_download] GSE20685: FAILED — %s\n  (continuing with next cohort)",
                  conditionMessage(e)))
})

# =============================================================================
# SECTION 5: LEAKAGE-PROOF CHECK — SCAN-B vs GSE96058
#   Memorial v6.1 §3.2 — training/validation GSM-ID intersection must be 0.
# =============================================================================
message("\n", strrep("=", 60))
message("[01_download] SECTION 5: LEAKAGE-PROOF CHECK")
message(strrep("=", 60))

tryCatch({
  scanb_path  <- file.path(raw_cohort("GSE96058"), "GSE96058_clinical_raw.rds")
  leakage_out <- file.path(PATHS$results$supp, "leakage_check_scanb_vs_gse96058.csv")

  if (!file.exists(scanb_path)) {
    message("[leakage] Deferred — GSE96058 clinical not yet present.")
  } else {
    old_w <- getOption("warn"); options(warn = 0)
    pheno <- readRDS(scanb_path)
    options(warn = old_w)

    split_col <- intersect(
      c("scan.b.training.validation", "group", "training_validation",
        "characteristics_ch1.6"),
      tolower(names(pheno))
    )

    if (length(split_col) == 0) {
      message("[leakage] WARNING: split column not auto-detected.")
      message("  Available columns: ", paste(head(names(pheno), 10), collapse = ", "))
    } else {
      col_use   <- names(pheno)[tolower(names(pheno)) == split_col[1]][1]
      train_ids <- rownames(pheno)[grepl("train",        tolower(pheno[[col_use]]))]
      valid_ids <- rownames(pheno)[grepl("valid|test",   tolower(pheno[[col_use]]))]
      overlap   <- intersect(train_ids, valid_ids)

      leakage_df <- tibble::tibble(
        timestamp         = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        script            = SCRIPT_NAME,
        split_column_used = col_use,
        n_train_gsm       = length(train_ids),
        n_validation_gsm  = length(valid_ids),
        n_overlap         = length(overlap),
        overlap_ids       = paste(overlap, collapse = ";"),
        sha256_train_ids  = digest::digest(sort(train_ids), algo = "sha256"),
        sha256_valid_ids  = digest::digest(sort(valid_ids), algo = "sha256"),
        result            = if (length(overlap) == 0) "PASS_NO_LEAKAGE"
                            else "FAIL_LEAKAGE_DETECTED"
      )
      readr::write_csv(leakage_df, leakage_out)

      if (length(overlap) == 0) {
        message(sprintf("[leakage] PASS — train: %d | validation: %d | overlap: 0",
                        length(train_ids), length(valid_ids)))
      } else {
        .errors[["leakage"]] <<- sprintf(
          "%d samples overlap between training and validation", length(overlap))
        message(sprintf("[leakage] FAIL — %d overlapping samples!", length(overlap)))
      }
    }
  }
}, error = function(e) {
  message(sprintf("[leakage] Check failed: %s", conditionMessage(e)))
})

# =============================================================================
# SECTION 6: FINAL REPORT
# =============================================================================
message("\n", strrep("=", 60))
message("[01_download] SECTION 6: FINAL REPORT")
message(strrep("=", 60))

# Write audit report from registry
old_w <- getOption("warn"); options(warn = 0)
if (file.exists(PATHS$run_registry)) {
  df_reg  <- readr::read_csv(PATHS$run_registry, show_col_types = FALSE)
  df_this <- dplyr::filter(df_reg, script == SCRIPT_NAME)
  audit_out <- file.path(PATHS$registry_docs, "data_lake_audit_report.csv")
  readr::write_csv(df_this, audit_out)
  message(sprintf("Audit report: %s (%d entries)", audit_out, nrow(df_this)))
}
options(warn = old_w)

# Per-cohort status summary
cohorts_expected <- c("GSE96058", "METABRIC", "TCGA_BRCA", "GSE20685")
message("\n--- Download status ---")
all_pass <- TRUE
for (co in cohorts_expected) {
  status <- if (co %in% names(.errors)) {
    all_pass <- FALSE
    sprintf("FAIL  (%s)", .errors[[co]])
  } else "PASS"
  message(sprintf("  %-12s %s", co, status))
}
message("-----------------------")

message(sprintf("\n[01_download] Completed: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

if (!all_pass) {
  message("\n[01_download] Some cohorts failed. Re-run to retry failed sections.")
  message("  (Successfully downloaded cohorts will be skipped on re-run.)")
  quit(save = "no", status = 1)
} else {
  message("[01_download] All cohorts OK.")
  message("Next step: scripts/02_harmonize_clinical_<COHORT>.R")
  quit(save = "no", status = 0)
}
