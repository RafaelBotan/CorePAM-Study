# Core-PAM: Minimum PAM50 Panel for Breast Cancer Prognosis

> **PhD thesis project — Rafael Botan**
> Reducing the PAM50 gene panel to the smallest subset (Core-PAM) that maintains non-inferior
> prognostic performance across five independent cohorts.

---

## What is Core-PAM?

PAM50 is a 50-gene expression panel widely used to classify breast cancer subtypes and predict
prognosis. Core-PAM asks: *how few of those 50 genes do we actually need?*

Using elastic-net Cox regression on the SCAN-B training cohort (n ≈ 3,000), this pipeline
derives the **minimum gene subset** whose prognostic C-index is within ΔC = 0.010 of the full
PAM50 model — then validates it across four independent cohorts covering RNA-seq and microarray
platforms.

---

## Requirements

| Component | Version |
|-----------|---------|
| R | ≥ 4.4.0 (developed on 4.5.2) |
| Quarto | ≥ 1.4 (for manuscript rendering) |
| Internet access | Required for data download (step 1) |
| Disk space | ≥ 15 GB for all raw data |
| RAM | ≥ 16 GB (TCGA-BRCA SummarizedExperiment) |

---

## Quick Start

### Step 0 — Install R packages

```r
source("INSTALL.R")
```

This installs all required CRAN and Bioconductor packages. Run once; safe to re-run.

### Step 1 — Run scripts in order

Open R (or RStudio) in the project root and run each script sequentially:

```r
# 0. Environment setup (sourced automatically by all scripts — do not run directly)
#    scripts/00_setup.R

# 1. Download all raw data (~10 GB, requires internet, ~1–4 hours depending on speed)
source("scripts/01_download_raw_data.R")

# 2. Harmonize clinical data (run for each cohort)
source("scripts/02_harmonize_clinical_SCANB.R")
source("scripts/02_harmonize_clinical_GSE96058.R")
source("scripts/02_harmonize_clinical_TCGA_BRCA.R")
source("scripts/02_harmonize_clinical_METABRIC.R")
source("scripts/02_harmonize_clinical_GSE20685.R")

# 3. Preprocess gene expression (run for each cohort)
source("scripts/03_expression_preprocess_SCANB.R")
source("scripts/03_expression_preprocess_GSE96058.R")
source("scripts/03_expression_preprocess_TCGA_BRCA.R")
source("scripts/03_expression_preprocess_METABRIC.R")
source("scripts/03_expression_preprocess_GSE20685.R")

# 4. PAM50 gene coverage audit (all cohorts, one script)
source("scripts/04_gene_audit_freeze.R")

# 5. Derive Core-PAM panel (training only — SCAN-B)
source("scripts/05_reduce_pam50_to_corepam_FINAL.R")

# 6. Calculate Core-PAM score per cohort
source("scripts/06_zscore_and_score_SCANB.R")
source("scripts/06_zscore_and_score_GSE96058.R")
source("scripts/06_zscore_and_score_TCGA_BRCA.R")
source("scripts/06_zscore_and_score_METABRIC.R")
source("scripts/06_zscore_and_score_GSE20685.R")

# 7. Pre-flight check, then survival analysis per cohort
source("scripts/07A_preflight_files_strict.R")
source("scripts/07_survival_analysis_SCANB.R")
source("scripts/07_survival_analysis_GSE96058.R")
source("scripts/07_survival_analysis_TCGA_BRCA.R")
source("scripts/07_survival_analysis_METABRIC.R")
source("scripts/07_survival_analysis_GSE20685.R")

# 8. Random-effects meta-analysis
source("scripts/08_meta_survival.R")

# 11. Incremental value and decision curve analysis
source("scripts/11_incremental_value_and_dca.R")

# 13–16. Quality control checks
source("scripts/13_qc_correlations_offdiag.R")
source("scripts/14_qc_metabric_pca_forensics.R")
source("scripts/15_qc_schema_range_checks.R")
source("scripts/16_qc_text_vs_results_assert.R")

# 17–18. Render manuscript and create submission bundle
source("scripts/17_render_manuscript_quarto.R")
source("scripts/18_make_submission_bundle.R")
```

> **Tip:** Every script checks whether its outputs already exist before running. To force
> re-execution of any script, set the environment variable `FORCE_RERUN=TRUE` before sourcing:
> ```r
> Sys.setenv(FORCE_RERUN = "TRUE")
> source("scripts/05_reduce_pam50_to_corepam_FINAL.R")
> Sys.unsetenv("FORCE_RERUN")
> ```

### Step 2 — Validate the pipeline

At any point, run the master validator to check syntax and output file presence:

```r
source("scripts/00_validate_all.R")
```

---

## Data Sources

All data are publicly available. The download script handles everything automatically.

| Cohort | Role | Platform | Accession | n (approx) |
|--------|------|----------|-----------|-----------|
| SCAN-B | **Training** | RNA-seq (Illumina) | [GSE96058](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058) | ~3 000 |
| GSE96058 | Validation | RNA-seq (Illumina) | [GSE96058](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058) | ~1 000 |
| TCGA-BRCA | Validation | RNA-seq (STAR) | GDC portal | ~1 100 |
| METABRIC | Validation | Microarray (Illumina HT-12) | [cBioPortal](https://www.cbioportal.org/study/summary?id=brca_metabric) | ~2 000 |
| GSE20685 | Validation | Microarray (Affymetrix HGU133A) | [GSE20685](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20685) | ~327 |

> **Note:** SCAN-B and GSE96058 share the same GEO accession (GSE96058). They are split
> deterministically by the `training/validation` column in the series matrix. There is zero
> sample overlap between training and validation (verified in script 01).

---

## Frozen Parameters

All analytical parameters were pre-specified and locked before any analysis.
They are stored in `01_docs/registry/analysis_freeze.csv` and loaded automatically.

| Parameter | Value | Meaning |
|-----------|-------|---------|
| `delta_c` | 0.010 | Non-inferiority margin (C-index OOF) |
| `alpha` | 0.5 | Elastic-net mixing (0 = ridge, 1 = lasso) |
| `k_folds` | 10 | Stratified cross-validation folds |
| `seed_folds` | 42 | Random seed for reproducibility |
| `min_genes_fraction` | 0.80 | Minimum panel coverage required per cohort |
| `bootstrap_n` | 1000 | Bootstrap resamples for 95% CI |

---

## Repository Structure

```
.
├── scripts/                  R analysis pipeline (00–18)
├── 01_docs/
│   └── registry/             Frozen parameters + artifact inventory
├── 01_Base_Pura_CorePAM/     Data lake (git-ignored)
│   ├── RAW/                  Immutable raw downloads (never modified)
│   └── PROCESSED/            Processed outputs per cohort
├── results/                  Analysis outputs (CSV, JSON, RDS)
│   ├── corepam/              Core-PAM model artifacts
│   ├── corepam_os/           Per-cohort survival results
│   ├── main/                 Main manuscript tables
│   └── supp/                 Supplementary tables & QC
├── figures/                  Generated plots (PDF, PNG)
│   ├── main/                 Main figures
│   └── supp/                 Supplementary figures
├── registry/
│   └── study_registry.csv    Append-only SHA-256 artifact log
├── logs/                     Execution logs
├── INSTALL.R                 Package installer (run first)
├── CLAUDE.md                 Full project specification and rules
├── Memorial_v6_1_CorePAM.md  Governance protocol (source of truth)
├── ERRORS_LOG.md             Error history with solutions
└── RESULTS_SUMMARY.md        Running results log
```

---

## Reproducibility Guarantees

- **Frozen parameters:** all analytical choices locked in `analysis_freeze.csv` before analysis
- **SHA-256 hashing:** every raw file and output artifact is hashed and registered in `registry/study_registry.csv`
- **No inter-cohort pooling:** expression matrices are never combined; Z-scoring is always intra-cohort
- **Strict I/O:** `options(warn = 2)` — any warning is treated as a fatal error
- **Skip/force pattern:** scripts skip existing outputs unless `FORCE_RERUN=TRUE`
- **Anti-hard-code check:** script 16 verifies that no number in the manuscript was typed manually

---

## Key Output Artifacts

After a complete run, the main deliverables are:

| File | Location | Description |
|------|----------|-------------|
| `CorePAM_weights.csv` | `results/corepam/` | Final gene panel with Cox weights |
| `CorePAM_model.rds` | `results/corepam/` | Fitted glmnet model object |
| `CorePAM_training_card.json` | `results/corepam/` | Derivation metadata (n, C-index, genes) |
| `survival_results_*.csv` | `results/corepam_os/` | Per-cohort Cox + C-index results |
| `meta_survival_summary.csv` | `results/corepam_os/` | Random-effects meta-analysis |
| `submission_bundle_*.zip` | `results/` | Complete reproducible submission package |

---

## Contact

Rafael Botan — PhD candidate
For issues with this pipeline, open a GitHub issue or check `ERRORS_LOG.md` for known problems.

---

## License

Data are from public repositories (GEO, GDC, cBioPortal) under their respective terms of use.
Analysis code: see repository license.
