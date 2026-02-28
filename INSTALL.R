# =============================================================================
# INSTALL.R — Core-PAM Pipeline Package Installer
# =============================================================================
# Run this script ONCE before executing any pipeline script.
# Safe to re-run: already-installed packages are skipped.
# Requires: R >= 4.4.0 and internet access.
# =============================================================================

cat("=============================================================\n")
cat("Core-PAM Pipeline — Package Installer\n")
cat("=============================================================\n\n")

# --- Check R version ---------------------------------------------------------
r_version <- as.numeric(paste0(R.version$major, ".", R.version$minor))
if (r_version < 4.4) {
  stop(sprintf(
    "R version %s detected. R >= 4.4.0 is required. Please upgrade R.",
    R.version.string
  ))
}
cat(sprintf("[OK] R version: %s\n\n", R.version.string))

# --- Install BiocManager if not present --------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("[INSTALL] Installing BiocManager...\n")
  install.packages("BiocManager", repos = "https://cran.r-project.org")
}
cat(sprintf("[OK] BiocManager: %s\n", packageVersion("BiocManager")))

# Set Bioconductor version matching R 4.5
BiocManager::install(version = "3.21", ask = FALSE, update = FALSE)

# --- CRAN packages -----------------------------------------------------------
cran_packages <- c(
  # Core data wrangling
  "tidyverse",    # dplyr, ggplot2, readr, tibble, tidyr, stringr, purrr
  "arrow",        # parquet read/write (analysis_ready.parquet files)
  "jsonlite",     # JSON read/write (training card, results)
  "digest",       # SHA-256 hashing (artifact integrity)
  # Statistical modelling
  "glmnet",       # Elastic-net Cox regression (Core-PAM derivation)
  "survival",     # Survival analysis (Cox models, KM)
  "survminer",    # Kaplan-Meier plots
  "metafor",      # Random-effects meta-analysis
  # Plotting
  "ggplot2",      # Base plotting (included in tidyverse, listed explicitly)
  "ggrepel",      # Non-overlapping labels on plots
  "gridExtra",    # Multi-panel figure layouts
  "reshape2",     # Data reshaping for heatmaps
  # Gene annotation
  "HGNChelper",   # HGNC gene symbol validation and correction
  # Network & web
  "httr",         # HTTP requests (data download)
  "curl",         # Low-level curl bindings
  # Utilities
  "here",         # Robust file path resolution
  "knitr",        # Report generation
  "quarto"        # Quarto document rendering
)

cat("\n[CRAN] Installing/checking CRAN packages...\n")
installed_cran <- rownames(installed.packages())
to_install_cran <- setdiff(cran_packages, installed_cran)

if (length(to_install_cran) > 0) {
  cat(sprintf("  Installing %d packages: %s\n",
              length(to_install_cran),
              paste(to_install_cran, collapse = ", ")))
  install.packages(to_install_cran,
                   repos     = "https://cran.r-project.org",
                   dependencies = TRUE)
} else {
  cat("  All CRAN packages already installed.\n")
}

# --- Bioconductor packages ---------------------------------------------------
bioc_packages <- c(
  "GEOquery",            # Download data from NCBI GEO (GSE96058, GSE20685)
  "TCGAbiolinks",        # Download TCGA-BRCA from GDC
  "SummarizedExperiment",# Container for RNA-seq data (TCGA)
  "Biobase",             # Base Bioconductor infrastructure (pData)
  "edgeR",               # TMM normalisation + logCPM (RNA-seq)
  "biomaRt",             # Ensembl → HGNC gene symbol mapping
  "AnnotationDbi",       # Annotation database interface
  "org.Hs.eg.db",        # Human gene annotation database (Affymetrix mapping)
  "cmprsk"               # Competing risks (Fine-Gray for METABRIC sensitivity)
)

cat("\n[Bioc] Installing/checking Bioconductor packages...\n")
installed_all <- rownames(installed.packages())
to_install_bioc <- setdiff(bioc_packages, installed_all)

if (length(to_install_bioc) > 0) {
  cat(sprintf("  Installing %d packages: %s\n",
              length(to_install_bioc),
              paste(to_install_bioc, collapse = ", ")))
  BiocManager::install(to_install_bioc, ask = FALSE, update = FALSE)
} else {
  cat("  All Bioconductor packages already installed.\n")
}

# --- Verification ------------------------------------------------------------
cat("\n[CHECK] Verifying all packages load correctly...\n")

all_packages <- c(cran_packages, bioc_packages)
# Remove meta-packages that don't load directly
all_packages <- setdiff(all_packages, c("tidyverse", "quarto"))

failed <- character(0)
for (pkg in all_packages) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  if (!ok) {
    failed <- c(failed, pkg)
    cat(sprintf("  [FAIL] %s\n", pkg))
  }
}

cat("\n=============================================================\n")
if (length(failed) == 0) {
  cat(sprintf("[SUCCESS] All %d packages verified.\n", length(all_packages)))
  cat("You are ready to run the Core-PAM pipeline.\n")
  cat("Next step: source('scripts/01_download_raw_data.R')\n")
} else {
  cat(sprintf("[WARNING] %d packages failed to load:\n", length(failed)))
  for (pkg in failed) cat(sprintf("  - %s\n", pkg))
  cat("\nTry reinstalling failed packages manually:\n")
  cat(sprintf(
    "  BiocManager::install(c(%s))\n",
    paste(sprintf('"%s"', failed), collapse = ", ")
  ))
}
cat("=============================================================\n")
