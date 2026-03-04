# =============================================================================
# SCRIPT: 00_setup.R
# PURPOSE: Reproducible environment — paths, strict I/O, helpers, frozen
#          parameters. Must be sourced at the beginning of all other scripts.
# PROJETO: Core-PAM (Memorial v6.1 / Freeze Core-PAM)
# REGRA:   Never contain gene counts (e.g. PAM29) in paths or variables.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(digest)
  library(arrow)      # parquet read/write
  library(jsonlite)   # training card JSON
})

# --------------------------------------------------------------------------
# 1) PATHS — single source of truth (adjust ROOT_REPO once)
# --------------------------------------------------------------------------
ROOT_REPO <- normalizePath(".", mustWork = TRUE)
DATA_LAKE <- file.path(ROOT_REPO, "01_Base_Pura_CorePAM")

PATHS <- list(
  raw            = file.path(DATA_LAKE, "RAW"),
  processed      = file.path(DATA_LAKE, "PROCESSED"),
  docs           = file.path(ROOT_REPO, "01_docs"),
  registry_docs  = file.path(ROOT_REPO, "01_docs", "registry"),
  scripts        = file.path(ROOT_REPO, "scripts"),
  config         = file.path(ROOT_REPO, "config"),
  results = list(
    corepam    = file.path(ROOT_REPO, "results", "corepam"),
    corepam_os = file.path(ROOT_REPO, "results", "corepam_os"),
    main       = file.path(ROOT_REPO, "results", "main"),
    supp       = file.path(ROOT_REPO, "results", "supp"),
    pcr        = file.path(ROOT_REPO, "results", "pcr")
  ),
  figures = list(
    # Base section dirs (backward-compat; prefer fig_dir() for new code)
    main   = file.path(ROOT_REPO, "figures", "main"),
    supp   = file.path(ROOT_REPO, "figures", "supp"),
    pcr    = file.path(ROOT_REPO, "figures", "pcr"),
    # Structured subdirs: figures/{section}/{lang}/{ext}/
    main_en_pdf = file.path(ROOT_REPO, "figures", "main", "en", "pdf"),
    main_en_png = file.path(ROOT_REPO, "figures", "main", "en", "png"),
    main_pt_pdf = file.path(ROOT_REPO, "figures", "main", "pt", "pdf"),
    main_pt_png = file.path(ROOT_REPO, "figures", "main", "pt", "png"),
    supp_en_pdf = file.path(ROOT_REPO, "figures", "supp", "en", "pdf"),
    supp_en_png = file.path(ROOT_REPO, "figures", "supp", "en", "png"),
    supp_pt_pdf = file.path(ROOT_REPO, "figures", "supp", "pt", "pdf"),
    supp_pt_png = file.path(ROOT_REPO, "figures", "supp", "pt", "png"),
    pcr_en_pdf  = file.path(ROOT_REPO, "figures", "pcr",  "en", "pdf"),
    pcr_en_png  = file.path(ROOT_REPO, "figures", "pcr",  "en", "png"),
    pcr_pt_pdf  = file.path(ROOT_REPO, "figures", "pcr",  "pt", "pdf"),
    pcr_pt_png  = file.path(ROOT_REPO, "figures", "pcr",  "pt", "png"),
    artigo = file.path(ROOT_REPO, "06_plots", "artigo"),
    tese   = file.path(ROOT_REPO, "06_plots", "tese")
  ),
  run_registry = file.path(ROOT_REPO, "registry", "study_registry.csv")
)

# Cohort path helpers — OS block
raw_cohort    <- function(cohort) file.path(PATHS$raw, cohort)
proc_cohort   <- function(cohort) file.path(PATHS$processed, cohort)

# Cohort path helpers — pCR block (separate sub-tree under RAW/pCR/ and PROCESSED/pCR/)
raw_pcr_cohort  <- function(cohort) file.path(PATHS$raw,       "pCR", cohort)
proc_pcr_cohort <- function(cohort) file.path(PATHS$processed, "pCR", cohort)

# Figure path helper — returns the figures/{section}/{lang}/{ext} directory
# section: "main" | "supp" | "pcr"
# lang:    "en"   | "pt"
# ext:     "pdf"  | "png"  | NULL (returns lang dir)
fig_dir <- function(section = "supp", lang = "en", ext = NULL) {
  base <- file.path(ROOT_REPO, "figures", section, lang)
  if (!is.null(ext)) base <- file.path(base, ext)
  dir.create(base, showWarnings = FALSE, recursive = TRUE)
  base
}

# Save figure as PDF + PNG into figures/{section}/{lang}/pdf|png/
# Returns list(pdf=..., png=...)  invisibly
fig_save <- function(plot, name, lang = "en", section = "supp",
                     w = 8, h = 6, dpi = 300) {
  pdf_path <- file.path(fig_dir(section, lang, "pdf"), paste0(name, ".pdf"))
  png_path <- file.path(fig_dir(section, lang, "png"), paste0(name, ".png"))
  old_warn <- getOption("warn"); options(warn = 0)
  pdf(pdf_path, width = w, height = h); print(plot); dev.off()
  png(png_path, width = round(w * dpi), height = round(h * dpi), res = dpi)
  print(plot); dev.off()
  options(warn = old_warn)
  invisible(list(pdf = pdf_path, png = png_path))
}

# --------------------------------------------------------------------------
# 2) STRICT I/O — warning = error (required; Memorial v6.1 §1.4 / §9.1)
# --------------------------------------------------------------------------
options(warn = 2)

# --------------------------------------------------------------------------
# 3) FROZEN PARAMETERS (analysis_freeze.csv — do not edit directly)
# --------------------------------------------------------------------------
freeze_path <- file.path(PATHS$registry_docs, "analysis_freeze.csv")

if (!file.exists(freeze_path)) {
  stop(
    "FREEZE NOT FOUND: ", freeze_path,
    "\nCreate 01_docs/registry/analysis_freeze.csv before any analysis."
  )
}

.freeze_raw <- read_csv(freeze_path, show_col_types = FALSE)
FREEZE <- setNames(
  lapply(.freeze_raw$value, function(x) {
    num <- suppressWarnings(as.numeric(x))
    if (!is.na(num)) num else x
  }),
  .freeze_raw$parameter
)
rm(.freeze_raw)

# Minimum validation
.required <- c("delta_c", "alpha", "k_folds", "seed_folds",
               "min_genes_fraction", "time_unit_divisor")
.missing  <- setdiff(.required, names(FREEZE))
if (length(.missing) > 0) {
  stop("Parameters missing in analysis_freeze.csv: ",
       paste(.missing, collapse = ", "))
}
rm(.required, .missing)

# --------------------------------------------------------------------------
# 4) HELPERS
# --------------------------------------------------------------------------

#' SHA-256 of a file — checks existence and size > 0
sha256_file <- function(path) {
  if (!file.exists(path))       stop("File not found: ", path)
  if (file.info(path)$size == 0) stop("Empty file (0 bytes): ", path)
  digest::digest(path, algo = "sha256", file = TRUE)
}

#' Normalize patient/sample ID
normalize_id <- function(x) trimws(toupper(as.character(x)))

#' Strict CSV read
strict_csv <- function(path, ...) {
  if (!file.exists(path)) stop("CSV not found: ", path)
  readr::read_csv(path, show_col_types = FALSE, ...)
}

#' Strict Parquet read
strict_parquet <- function(path) {
  if (!file.exists(path)) stop("Parquet not found: ", path)
  if (file.info(path)$size == 0) stop("Empty parquet (0 bytes): ", path)
  arrow::read_parquet(path)
}

#' Strict RDS read
strict_rds <- function(path) {
  if (!file.exists(path)) stop("RDS not found: ", path)
  if (file.info(path)$size == 0) stop("Empty RDS (0 bytes): ", path)
  # Temporarily lower warn to avoid converting harmless decompression warnings
  # to hard errors under the global warn=2 setting.
  old_warn <- getOption("warn"); on.exit(options(warn = old_warn))
  options(warn = 0)
  readRDS(path)
}

#' Append to registry (append-only; creates header if file does not yet exist)
registry_append <- function(cohort, file_type, file_path, sha256,
                            status, script, size_mb = NA_real_,
                            extra = list()) {
  path <- PATHS$run_registry
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)

  row <- tibble(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    script    = script,
    cohort    = cohort,
    file_type = file_type,
    path      = file_path,
    sha256    = sha256,
    status    = status,
    size_mb   = round(size_mb, 3)
  )
  if (length(extra) > 0) row <- bind_cols(row, as_tibble(extra))

  write_csv(row, path,
            append   = file.exists(path),
            col_names = !file.exists(path))
  invisible(row)
}

# --------------------------------------------------------------------------
# 5) DIRECTORY STRUCTURE VALIDATION (creates if missing)
# --------------------------------------------------------------------------
.required_dirs <- c(
  PATHS$raw, PATHS$processed,
  PATHS$docs, PATHS$registry_docs,
  PATHS$scripts,
  PATHS$results$corepam, PATHS$results$corepam_os,
  PATHS$results$main,    PATHS$results$supp,    PATHS$results$pcr,
  PATHS$config,
  PATHS$figures$artigo,  PATHS$figures$tese,
  dirname(PATHS$run_registry),
  # Structured figure subdirectories: figures/{section}/{lang}/{ext}/
  PATHS$figures$main_en_pdf, PATHS$figures$main_en_png,
  PATHS$figures$main_pt_pdf, PATHS$figures$main_pt_png,
  PATHS$figures$supp_en_pdf, PATHS$figures$supp_en_png,
  PATHS$figures$supp_pt_pdf, PATHS$figures$supp_pt_png,
  PATHS$figures$pcr_en_pdf,  PATHS$figures$pcr_en_png,
  PATHS$figures$pcr_pt_pdf,  PATHS$figures$pcr_pt_png
)
for (.d in .required_dirs) dir.create(.d, showWarnings = FALSE, recursive = TRUE)
rm(.required_dirs, .d)

# --------------------------------------------------------------------------
message(sprintf(
  "[00_setup] OK | ROOT: %s | delta_c=%.3f | alpha=%.1f | K=%d | seed=%d",
  ROOT_REPO, FREEZE$delta_c, FREEZE$alpha, FREEZE$k_folds, FREEZE$seed_folds
))
