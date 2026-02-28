# =============================================================================
# SCRIPT: 18_make_submission_bundle.R
# PURPOSE: Create a self-contained, re-executable submission bundle:
#          - All scripts (00-18)
#          - All result CSVs and JSONs (no raw data)
#          - All figures (PDF + PNG)
#          - Registry snapshot
#          - SHA-256 manifest of all bundled files
#          - README with exact reproduction steps
# PROJECT: Core-PAM (Memorial v6.1 §Submission)
#
# OUTPUTS:
#   submission_bundle/CorePAM_bundle_YYYYMMDD.zip
#   submission_bundle/CorePAM_manifest_YYYYMMDD.csv
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "18_make_submission_bundle.R"

BUNDLE_DATE <- format(Sys.Date(), "%Y%m%d")
BUNDLE_NAME <- sprintf("CorePAM_bundle_%s", BUNDLE_DATE)

# Skip if bundle already exists (set FORCE_RERUN=TRUE to recreate)
FORCE <- as.logical(Sys.getenv("FORCE_RERUN", "FALSE"))
bundle_dir <- file.path(ROOT_REPO, "submission_bundle")
zip_path   <- file.path(bundle_dir, sprintf("%s.zip", BUNDLE_NAME))

if (!FORCE && file.exists(zip_path)) {
  message(sprintf("[18] Bundle already exists: %s", zip_path))
  message("[18] Skipping. Set FORCE_RERUN=TRUE to recreate.")
  quit(save = "no", status = 0)
}

message(sprintf("[18] Creating submission bundle: %s", BUNDLE_NAME))

dir.create(bundle_dir, showWarnings = FALSE, recursive = TRUE)

# Staging directory inside bundle_dir
stage_dir <- file.path(bundle_dir, BUNDLE_NAME)
if (dir.exists(stage_dir)) unlink(stage_dir, recursive = TRUE)
dir.create(stage_dir, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Helper: copy a file to staging, preserving relative path structure
# ---------------------------------------------------------------------------
stage_file <- function(src_path, base_dir = ROOT_REPO) {
  if (!file.exists(src_path)) return(invisible(NULL))
  rel_path <- sub(paste0("^", normalizePath(base_dir, winslash = "/"), "/?"), "",
                  normalizePath(src_path, winslash = "/"))
  dst_path <- file.path(stage_dir, rel_path)
  dir.create(dirname(dst_path), showWarnings = FALSE, recursive = TRUE)
  file.copy(src_path, dst_path, overwrite = TRUE)
  invisible(dst_path)
}

manifest_rows <- list()

record_file <- function(src_path, category) {
  if (!file.exists(src_path)) return(invisible(NULL))
  h   <- sha256_file(src_path)
  sz  <- file.info(src_path)$size
  rel <- sub(paste0("^", normalizePath(ROOT_REPO, winslash = "/"), "/?"), "",
             normalizePath(src_path, winslash = "/"))
  manifest_rows[[length(manifest_rows) + 1]] <<- list(
    category  = category,
    file      = rel,
    sha256    = h,
    size_bytes = sz,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  stage_file(src_path)
  message(sprintf("[18] Staged [%s]: %s", category, basename(src_path)))
}

# ---------------------------------------------------------------------------
# 1) Scripts (all *.R)
# ---------------------------------------------------------------------------
message("[18] Staging scripts...")
script_files <- list.files(file.path(ROOT_REPO, "scripts"), pattern = "\\.R$",
                            full.names = TRUE)
for (f in script_files) record_file(f, "scripts")

# ---------------------------------------------------------------------------
# 2) Documentation
# ---------------------------------------------------------------------------
message("[18] Staging documentation...")
doc_files <- c(
  file.path(ROOT_REPO, "CLAUDE.md"),
  file.path(ROOT_REPO, "METHODS_DRAFT.md"),
  file.path(ROOT_REPO, "RESULTS_SUMMARY.md"),
  list.files(file.path(ROOT_REPO, "01_docs"), pattern = "\\.(md|csv|json)$",
             recursive = TRUE, full.names = TRUE)
)
for (f in doc_files[file.exists(doc_files)]) record_file(f, "documentation")

# ---------------------------------------------------------------------------
# 3) Analysis freeze
# ---------------------------------------------------------------------------
message("[18] Staging analysis freeze...")
freeze_path <- file.path(ROOT_REPO, "01_docs", "registry", "analysis_freeze.csv")
if (file.exists(freeze_path)) record_file(freeze_path, "freeze")

# ---------------------------------------------------------------------------
# 4) CorePAM artifacts (no .rds — weights and card only)
# ---------------------------------------------------------------------------
message("[18] Staging CorePAM artifacts...")
corepam_files <- list.files(PATHS$results$corepam,
                              pattern = "\\.(csv|json)$",
                              full.names = TRUE)
for (f in corepam_files) record_file(f, "corepam_artifacts")

# ---------------------------------------------------------------------------
# 5) Results (CSV only — no .rds, no .parquet)
# ---------------------------------------------------------------------------
message("[18] Staging results (CSV only)...")
result_dirs <- c(PATHS$results$main, PATHS$results$supp)
for (rdir in result_dirs) {
  if (!dir.exists(rdir)) next
  csv_files <- list.files(rdir, pattern = "\\.csv$", full.names = TRUE)
  for (f in csv_files) record_file(f, "results_csv")
}

# ---------------------------------------------------------------------------
# 6) Figures (PDF + PNG)
# ---------------------------------------------------------------------------
message("[18] Staging figures...")
figure_dirs <- c(PATHS$figures$main, PATHS$figures$supp)
for (fdir in figure_dirs) {
  if (!dir.exists(fdir)) next
  fig_files <- list.files(fdir, pattern = "\\.(pdf|png)$", full.names = TRUE)
  for (f in fig_files) record_file(f, "figures")
}

# ---------------------------------------------------------------------------
# 7) Registry snapshot (copy, not the gitignored original)
# ---------------------------------------------------------------------------
registry_path <- file.path(ROOT_REPO, "registry", "study_registry.csv")
if (file.exists(registry_path)) record_file(registry_path, "registry_snapshot")

# ---------------------------------------------------------------------------
# 8) Manuscript (if present)
# ---------------------------------------------------------------------------
ms_files <- list.files(file.path(ROOT_REPO, "manuscript"),
                        pattern = "\\.(qmd|pdf|html)$",
                        full.names = TRUE, recursive = TRUE)
for (f in ms_files) record_file(f, "manuscript")

# ---------------------------------------------------------------------------
# 9) Create README with reproduction steps
# ---------------------------------------------------------------------------
readme_path <- file.path(stage_dir, "README_REPRODUCTION.md")
readme_text <- sprintf(
'# Core-PAM Submission Bundle — %s

## Reproduction Instructions

1. **Install R** (≥ 4.4) and required packages (see `scripts/00_setup.R`).
2. **Install Quarto** (≥ 1.4) for manuscript rendering.
3. **Set working directory** to the repository root.
4. **Run scripts in order:**

```bash
Rscript scripts/01_download_raw_data.R
Rscript scripts/02_harmonize_clinical_SCANB.R
Rscript scripts/02_harmonize_clinical_GSE96058.R
Rscript scripts/02_harmonize_clinical_TCGA_BRCA.R
Rscript scripts/02_harmonize_clinical_METABRIC.R
Rscript scripts/02_harmonize_clinical_GSE20685.R
Rscript scripts/03_expression_preprocess_SCANB.R
Rscript scripts/03_expression_preprocess_GSE96058.R
Rscript scripts/03_expression_preprocess_TCGA_BRCA.R
Rscript scripts/03_expression_preprocess_METABRIC.R
Rscript scripts/03_expression_preprocess_GSE20685.R
Rscript scripts/04_gene_audit_freeze.R
Rscript scripts/05_reduce_pam50_to_corepam_FINAL.R
Rscript scripts/06_zscore_and_score_SCANB.R
Rscript scripts/06_zscore_and_score_GSE96058.R
Rscript scripts/06_zscore_and_score_TCGA_BRCA.R
Rscript scripts/06_zscore_and_score_METABRIC.R
Rscript scripts/06_zscore_and_score_GSE20685.R
Rscript scripts/07A_preflight_files_strict.R
Rscript scripts/07_survival_analysis_SCANB.R
Rscript scripts/07_survival_analysis_GSE96058.R
Rscript scripts/07_survival_analysis_TCGA_BRCA.R
Rscript scripts/07_survival_analysis_METABRIC.R
Rscript scripts/07_survival_analysis_GSE20685.R
Rscript scripts/08_meta_survival.R
Rscript scripts/11_incremental_value_and_dca.R
Rscript scripts/13_qc_correlations_offdiag.R
Rscript scripts/14_qc_metabric_pca_forensics.R
Rscript scripts/15_qc_schema_range_checks.R
Rscript scripts/16_qc_text_vs_results_assert.R
Rscript scripts/17_render_manuscript_quarto.R
```

## Frozen Parameters

See `01_docs/registry/analysis_freeze.csv` for all frozen model parameters.

## Bundle Contents

- `scripts/`        — All analysis scripts (00-18)
- `results/`        — Result CSVs (no raw data)
- `figures/`        — All figures (PDF + PNG)
- `manuscript/`     — Manuscript QMD and rendered outputs
- `registry/`       — Study registry snapshot (SHA-256 audit trail)
- `01_docs/`        — Freeze parameters, artifact inventory

## Integrity

All files listed in `CorePAM_manifest_%s.csv` with SHA-256 checksums.

## Contact

Rafael Botan — PhD Candidate
GitHub: https://github.com/RafaelBotan/Phd-Genomic
',
  BUNDLE_DATE, BUNDLE_DATE
)
writeLines(readme_text, readme_path)
record_file(readme_path, "readme")

# ---------------------------------------------------------------------------
# 10) Save manifest CSV
# ---------------------------------------------------------------------------
manifest_df <- bind_rows(lapply(manifest_rows, as_tibble))
manifest_path_staged <- file.path(stage_dir,
                                   sprintf("CorePAM_manifest_%s.csv", BUNDLE_DATE))
readr::write_csv(manifest_df, manifest_path_staged)
message(sprintf("[18] Manifest: %d files staged", nrow(manifest_df)))

# ---------------------------------------------------------------------------
# 11) Create ZIP archive
# ---------------------------------------------------------------------------
message("[18] Creating ZIP archive...")

old_wd <- getwd()
setwd(bundle_dir)
old_warn <- getOption("warn"); options(warn = 0)
zip_result <- tryCatch(
  zip(zipfile = basename(zip_path),
      files   = BUNDLE_NAME,
      flags   = "-r"),
  error = function(e) {
    # Fallback: use system zip
    system2("zip", c("-r", basename(zip_path), BUNDLE_NAME),
            stdout = FALSE, stderr = FALSE)
  }
)
options(warn = old_warn)
setwd(old_wd)

if (file.exists(zip_path)) {
  h_zip <- sha256_file(zip_path)
  sz_zip <- file.info(zip_path)$size / 1e6
  message(sprintf("[18] Bundle ZIP: %s (%.1f MB | SHA256: %s)", zip_path, sz_zip, h_zip))
  registry_append("ALL", "submission_bundle", zip_path, h_zip, "ok",
                  SCRIPT_NAME, sz_zip)
} else {
  message("[18] WARNING: ZIP file not created. Staged files are in: ", stage_dir)
}

# ---------------------------------------------------------------------------
# 12) Save manifest to results/supp
# ---------------------------------------------------------------------------
manifest_final <- file.path(PATHS$results$supp,
                              sprintf("CorePAM_manifest_%s.csv", BUNDLE_DATE))
dir.create(PATHS$results$supp, showWarnings = FALSE, recursive = TRUE)
readr::write_csv(manifest_df, manifest_final)
h_mf <- sha256_file(manifest_final)
registry_append("ALL", "submission_manifest", manifest_final, h_mf, "ok",
                SCRIPT_NAME, file.info(manifest_final)$size / 1e6)

message(sprintf("[18] Submission bundle complete: %d files | %s",
                nrow(manifest_df), BUNDLE_DATE))
message("[18] Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
