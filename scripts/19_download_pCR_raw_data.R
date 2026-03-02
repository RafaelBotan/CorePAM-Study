# =============================================================================
# SCRIPT: 19_download_pCR_raw_data.R
# PURPOSE: Download raw GEO data for all pCR-block cohorts (NACT response):
#            GSE25066  (Hatzis 2011)  — HGU133Plus2, N=508
#            GSE20194  (Mina)         — HGU133Plus2, N=278
#            GSE32646  (Tabchy 2010)  — HGU133Plus2, N=154
#            GSE22226  (I-SPY1)       — Agilent 44K, N=149
#          Downloads series matrices via GEOquery; skips if already present.
#          SHA-256 + registry for all downloaded files.
#
# PROJETO: Core-PAM (Memorial v6.1 / Freeze Core-PAM — pCR block)
# EXECUTION:
#   Rscript scripts/19_download_pCR_raw_data.R
#   Force re-download: FORCE_RERUN=TRUE Rscript scripts/19_download_pCR_raw_data.R
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "19_download_pCR_raw_data.R"
FORCE       <- as.logical(Sys.getenv("FORCE_RERUN", "FALSE"))

old_warn <- getOption("warn"); options(warn = 0)
suppressPackageStartupMessages(library(GEOquery))
options(warn = old_warn)

message(sprintf("[%s] Starting pCR cohort downloads | FORCE=%s", SCRIPT_NAME, FORCE))

# --------------------------------------------------------------------------
# pCR cohort registry
# --------------------------------------------------------------------------
PCR_COHORTS <- list(
  list(cohort    = "GSE25066",
       accession = "GSE25066",
       platform  = "hgu133plus2",
       note      = "Hatzis 2011; Lancet Oncol; N=508; NACT (EC ± docetaxel)"),
  list(cohort    = "GSE20194",
       accession = "GSE20194",
       platform  = "hgu133plus2",
       note      = "Park 2006 + NCNN; HGU133Plus2; N=278; NACT (paclitaxel/FAC)"),
  list(cohort    = "GSE32646",
       accession = "GSE32646",
       platform  = "hgu133plus2",
       note      = "Tabchy 2010; N=154; NACT (paclitaxel+FAC)"),
  list(cohort    = "ISPY1",
       accession = "GSE22226",
       platform  = "hgu133a",
       note      = "I-SPY1 / Rubin 2011; Agilent 44K; N=149; NACT (AC→T)")
)

.errors <- list()

for (cdef in PCR_COHORTS) {
  cohort    <- cdef$cohort
  accession <- cdef$accession
  platform  <- cdef$platform

  dest_dir <- raw_pcr_cohort(cohort)
  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)

  # Expected series matrix path
  sm_path <- file.path(dest_dir, sprintf("%s_series_matrix.txt.gz", accession))

  if (!FORCE && file.exists(sm_path) && file.info(sm_path)$size > 1e5) {
    message(sprintf("[%s] [%s] Series matrix already present (%s) — skipping.",
                    SCRIPT_NAME, cohort, basename(sm_path)))
    next
  }

  message(sprintf("[%s] [%s] Downloading %s series matrix from GEO...",
                  SCRIPT_NAME, cohort, accession))

  tryCatch({
    old_warn2 <- getOption("warn"); options(warn = 0)
    GEOquery::getGEO(
      GEO     = accession,
      destdir = dest_dir,
      GSEMatrix = TRUE,
      getGPL    = FALSE,
      AnnotGPL  = FALSE
    )
    options(warn = old_warn2)

    # Locate downloaded series matrix (may be named differently)
    sm_files <- list.files(dest_dir, pattern = "series_matrix.txt.gz",
                           full.names = TRUE)
    if (length(sm_files) == 0) {
      stop("Series matrix not found after GEOquery download.")
    }
    sm_actual <- sm_files[1]

    # Normalise filename if GEOquery saved with a different name
    if (sm_actual != sm_path) {
      file.rename(sm_actual, sm_path)
    }

    h    <- sha256_file(sm_path)
    sz   <- file.info(sm_path)$size / 1e6
    registry_append(cohort, "series_matrix_gz", sm_path, h, "ok",
                    SCRIPT_NAME, sz)
    message(sprintf("[%s] [%s] Downloaded: %s (%.1f MB | SHA256: %s)",
                    SCRIPT_NAME, cohort, basename(sm_path), sz, h))

  }, error = function(e) {
    .errors[[cohort]] <<- conditionMessage(e)
    message(sprintf("[%s] [%s] ERROR: %s", SCRIPT_NAME, cohort,
                    conditionMessage(e)))
  })
}

# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------
message(sprintf("\n[%s] === DOWNLOAD SUMMARY ===", SCRIPT_NAME))
for (cdef in PCR_COHORTS) {
  sm_path <- file.path(raw_pcr_cohort(cdef$cohort),
                       sprintf("%s_series_matrix.txt.gz", cdef$accession))
  ok <- file.exists(sm_path) && file.info(sm_path)$size > 1e5
  message(sprintf("  %-10s %-12s  %s",
                  cdef$cohort, cdef$accession, if (ok) "OK" else "FAIL/MISSING"))
}

if (length(.errors) > 0) {
  message("\n[%s] ERRORS encountered:", SCRIPT_NAME)
  for (nm in names(.errors)) message(sprintf("  %s: %s", nm, .errors[[nm]]))
  warning(sprintf("[%s] %d cohort(s) failed. Re-run to retry.", SCRIPT_NAME,
                  length(.errors)), call. = FALSE)
} else {
  message(sprintf("[%s] COMPLETED — all pCR cohorts downloaded.", SCRIPT_NAME))
}
