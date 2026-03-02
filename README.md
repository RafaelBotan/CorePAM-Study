# CorePAM: A PAM50-Derived Minimum Gene Expression Score for Breast Cancer Prognosis

[![DOI](https://img.shields.io/badge/status-under%20review-yellow)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![R ≥ 4.4](https://img.shields.io/badge/R-%E2%89%A5%204.4-276DC3?logo=r)](https://cran.r-project.org/)

> **PhD thesis project — Rafael de Negreiros Botan**
>
> Reducing the 50-gene PAM50 panel to the smallest prognostic subset (CorePAM, 24 genes)
> whose out-of-fold C-index is non-inferior (ΔC ≤ 0.010) to the full 50-gene elastic-net
> Cox model. Validated across four independent cohorts spanning RNA-seq and microarray platforms.

---

## Overview

PAM50 is the gold-standard 50-gene expression panel for breast cancer molecular classification.
CorePAM asks: **how few of those 50 genes are actually needed for prognostic scoring?**

Using elastic-net Cox regression on the SCAN-B cohort (n = 3,069), we derived a **24-gene
subset** that retains non-inferior prognostic discrimination, then validated it on three
independent cohorts (TCGA-BRCA, METABRIC, GSE20685) covering both RNA-seq and microarray
platforms. A secondary analysis evaluated the CorePAM score as a predictor of pathological
complete response (pCR) to neoadjuvant chemotherapy across five additional cohorts.

**Key results:**
- Random-effects meta-analysis HR = 1.61 (95% CI 1.17–2.22) per 1-SD increase in CorePAM score
- Concordance (C-index) 0.58–0.64 across validation cohorts
- pCR meta-analysis OR = 1.69 (95% CI 1.39–2.05), I² = 0%

---

## Cohorts

### Survival (Overall Survival / Disease-Specific Survival)

| Cohort | Role | Platform | Endpoint | N | Events | Median FU |
|--------|------|----------|----------|---|--------|-----------|
| SCAN-B (GSE96058) | **Training** | RNA-seq | OS | 3,069 | 322 (10.5%) | 54.9 mo |
| TCGA-BRCA | Validation | RNA-seq | OS | 1,072 | 150 (14.0%) | 25.6 mo |
| METABRIC | Validation | Microarray (Illumina) | DSS | 1,980 | 1,144 (57.8%) | 157.9 mo |
| GSE20685 | Validation | Microarray (Affymetrix) | OS | 327 | 83 (25.4%) | 97.2 mo |

### pCR (Pathological Complete Response to NACT)

| Cohort | Role | Platform | N | pCR rate |
|--------|------|----------|---|----------|
| GSE25066 | Primary | HGU133A | 508 | 15.4% |
| GSE20194 | Primary | HGU133Plus2 | 278 | 25.9% |
| GSE32646 | Primary | HGU133Plus2 | 154 | 19.5% |
| I-SPY 1 (GSE22226) | Primary | Agilent 44K | 122 | 27.9% |
| I-SPY 2 (GSE194040) | Exploratory | Agilent 44K | 986 | 32.4% |

---

## Repository Structure

```
.
├── scripts/                     54 R scripts (numbered 00–23b)
│   ├── 00_setup.R               Environment + frozen parameters (sourced by all)
│   ├── 00_colors.R              Colour palette for all figures
│   ├── 01–04_*.R                Data download, clinical harmonisation, expression preprocessing
│   ├── 05_reduce_pam50_to_corepam_FINAL.R   CorePAM derivation (FREEZE)
│   ├── 06_*.R                   Z-score and scoring per cohort
│   ├── 07_*.R                   Survival analysis per cohort + figures
│   ├── 08_meta_survival.R       Random-effects meta-analysis (validation only)
│   ├── 11_*.R                   Incremental value (ΔC-index) + DCA
│   ├── 13–16_*.R                Quality control (correlations, PCA, schema, anti-hard-code)
│   ├── 17–18_*.R                Manuscript rendering + submission bundle
│   └── 19–23b_*.R               pCR analysis block (download → analysis → meta → figures)
│
├── results/
│   ├── corepam/                 Model artefacts (weights, Pareto, training card)
│   ├── main/                    Manuscript tables (Table 3, meta-analysis summary)
│   ├── supp/                    Supplementary tables + QC reports
│   └── pcr/                     pCR logistic results + meta-analysis
│
├── figures/
│   ├── main/{en,pt}/{pdf,png}/  Main figures (EN + PT versions, vector + raster)
│   ├── supp/{en,pt}/{pdf,png}/  Supplementary figures
│   ├── pcr/{en,pt}/{pdf,png}/   pCR figures
│   └── interpretations/         Markdown narrative for each figure
│
├── 01_docs/
│   ├── registry/                Frozen parameters, artefact inventory, audit reports
│   └── endpoint_mapping_templates/  Clinical endpoint mapping per cohort
│
├── manuscript/
│   ├── CorePAM_manuscript.qmd   Quarto source (all numbers read from CSVs)
│   └── CorePAM_manuscript.html  Rendered HTML
│
├── config/                      pCR cohort manifest
├── docs/                        Figures/tables index, reviewer response letter
│
├── CLAUDE.md                    Full project specification and analytical rules
├── Memorial_v6_1_CorePAM.md     Governance protocol (source of truth)
├── COREPAM_REPRO_RUNBOOK_v2.md  Reproducibility runbook
├── CITATION.cff                 Citation metadata
├── INSTALL.R                    One-step package installer
├── MANIFEST.txt                 File manifest
└── README.md                    This file
```

> **Data lake** (`01_Base_Pura_CorePAM/`) is git-ignored. All raw and processed data
> are downloaded automatically by `scripts/01_download_raw_data.R`.

---

## Frozen Parameters

All analytical parameters were pre-specified and locked before any data analysis.
They are stored in [`01_docs/registry/analysis_freeze.csv`](01_docs/registry/analysis_freeze.csv).

| Parameter | Value | Description |
|-----------|-------|-------------|
| `delta_c` | 0.010 | Non-inferiority margin for C-index |
| `alpha` | 0.5 | Elastic-net mixing parameter |
| `k_folds` | 10 | Stratified CV folds |
| `seed_folds` | 42 | Fold assignment seed |
| `min_genes_fraction` | 0.80 | Minimum gene coverage per cohort |
| `time_unit_divisor` | 30.4375 | Days to months conversion |

---

## Reproducibility

### Requirements

| Component | Version |
|-----------|---------|
| R | ≥ 4.4.0 (developed on 4.5.2) |
| Quarto | ≥ 1.4 (for manuscript rendering only) |
| Internet | Required for initial data download |
| Disk space | ≥ 15 GB |
| RAM | ≥ 16 GB |

### Step 0 — Install packages

```r
source("INSTALL.R")
```

### Step 1 — Run the full pipeline

Scripts are numbered and must be executed in order. Each script sources
`00_setup.R` automatically and checks whether its outputs already exist
(skip/force pattern — set `FORCE_RERUN=TRUE` to re-execute).

```r
# Download raw data (requires internet, ~1–4 hours)
source("scripts/01_download_raw_data.R")

# Clinical harmonisation (4 cohorts)
source("scripts/02_harmonize_clinical_SCANB.R")
source("scripts/02_harmonize_clinical_TCGA_BRCA.R")
source("scripts/02_harmonize_clinical_METABRIC.R")
source("scripts/02_harmonize_clinical_GSE20685.R")

# Expression preprocessing (4 cohorts)
source("scripts/03_expression_preprocess_SCANB.R")
source("scripts/03_expression_preprocess_TCGA_BRCA.R")
source("scripts/03_expression_preprocess_METABRIC.R")
source("scripts/03_expression_preprocess_GSE20685.R")

# Gene audit + CorePAM derivation
source("scripts/04_gene_audit_freeze.R")
source("scripts/05_reduce_pam50_to_corepam_FINAL.R")

# Scoring (4 cohorts)
source("scripts/06_zscore_and_score_SCANB.R")
source("scripts/06_zscore_and_score_TCGA_BRCA.R")
source("scripts/06_zscore_and_score_METABRIC.R")
source("scripts/06_zscore_and_score_GSE20685.R")

# Survival analysis (preflight + 4 cohorts + meta)
source("scripts/07A_preflight_files_strict.R")
source("scripts/07_survival_analysis_SCANB.R")
source("scripts/07_survival_analysis_TCGA_BRCA.R")
source("scripts/07_survival_analysis_METABRIC.R")
source("scripts/07_survival_analysis_GSE20685.R")
source("scripts/08_meta_survival.R")

# Figures
source("scripts/07x_extra_figures.R")
source("scripts/07y_pam50full_comparison.R")
source("scripts/07z_table_figures_survival.R")

# Incremental value + DCA
source("scripts/11_incremental_value_and_dca.R")
source("scripts/11b_dca_corepam.R")

# Quality control
source("scripts/13_qc_correlations_offdiag.R")
source("scripts/14_qc_metabric_pca_forensics.R")
source("scripts/15_qc_schema_range_checks.R")
source("scripts/16_qc_text_vs_results_assert.R")

# Manuscript rendering
source("scripts/17_render_manuscript_quarto.R")

# --- pCR block (secondary analysis) ---
source("scripts/19_download_pCR_raw_data.R")
source("scripts/20_prepare_pCR_GSE25066.R")
source("scripts/20_prepare_pCR_GSE20194.R")
source("scripts/20_prepare_pCR_GSE32646.R")
source("scripts/20_prepare_pCR_ISPY1.R")
source("scripts/20_prepare_pCR_ISPY2.R")
source("scripts/21_pCR_logistic_analysis.R")
source("scripts/21c_pcr_er_interaction.R")
source("scripts/21_ispy2_analysis.R")
source("scripts/22_meta_pCR.R")
source("scripts/22b_meta_pCR_with_ispy2.R")
source("scripts/23_pCR_figures.R")
source("scripts/23b_pCR_ispy2_figures.R")
```

### Step 2 — Validate

```r
source("scripts/00_validate_all.R")
```

This checks syntax of all scripts and verifies that expected output files exist.

---

## Reproducibility Safeguards

| Safeguard | Implementation |
|-----------|----------------|
| **Frozen parameters** | All choices locked in `analysis_freeze.csv` before analysis |
| **SHA-256 audit trail** | Every raw file and output artefact hashed in `registry/study_registry.csv` |
| **No inter-cohort pooling** | Expression matrices never combined; Z-scoring always intra-cohort |
| **Strict I/O** | `options(warn = 2)` — any warning is a fatal error |
| **Skip/force pattern** | Scripts check for existing outputs; `FORCE_RERUN=TRUE` to override |
| **Anti-hard-code check** | Script 16 verifies no manuscript number was typed manually |
| **Leakage-proof** | Training (SCAN-B) and validation cohorts are never mixed |
| **No coefficient refitting** | Validation cohorts use training weights only |

---

## Key Output Artefacts

| Artefact | Path | Description |
|----------|------|-------------|
| `CorePAM_weights.csv` | `results/corepam/` | 24 genes with Cox elastic-net weights |
| `CorePAM_model.rds` | `results/corepam/` | Fitted glmnet model (git-ignored; regenerated by script 05) |
| `selected_CorePAM_summary.json` | `results/corepam/` | Derivation metadata (C-index, genes, lambda) |
| `pareto_df_cindex_oof.csv` | `results/corepam/` | Pareto frontier: gene count vs. OOF C-index |
| `Table3_*.csv` | `results/main/` | Survival performance by cohort (HR, C-index, events) |
| `meta_survival_summary.csv` | `results/main/` | Random-effects meta-analysis results |
| `pCR_results_by_cohort.csv` | `results/pcr/` | Logistic regression per pCR cohort |
| `meta_pCR_results.csv` | `results/pcr/` | pCR meta-analysis (OR, I²) |

---

## CorePAM Score Formula

For a given sample in cohort *c*:

```
G_present = CorePAM genes ∩ genes available in cohort c
z_i       = intra-cohort Z-score for gene i
score     = Σ(w_i · z_i) / Σ|w_i|     (over G_present)
score_z   = scale(score)                (for HR per 1-SD)
```

Weights (`w_i`) are fixed from training. Score direction is established in SCAN-B
and never re-inverted per cohort.

---

## Data Sources

All data are publicly available:

- **SCAN-B**: [GSE96058](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058)
- **TCGA-BRCA**: [GDC Data Portal](https://portal.gdc.cancer.gov/)
- **METABRIC**: [cBioPortal](https://www.cbioportal.org/study/summary?id=brca_metabric)
- **GSE20685**: [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20685)
- **GSE25066**: [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25066)
- **GSE20194**: [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20194)
- **GSE32646**: [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32646)
- **I-SPY 1**: [GSE22226](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22226)
- **I-SPY 2**: [GSE194040](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194040)

---

## Citation

If you use this work, please cite:

> Botan RN, Sousa JB. CorePAM: a PAM50-derived gene expression score validated across
> RNA-seq and microarray platforms for breast cancer prognosis. *Breast Cancer Res*. 2026
> (under review).

See also [`CITATION.cff`](CITATION.cff) for machine-readable citation metadata.

---

## License

Analysis code: MIT License.
Data are from public repositories (GEO, GDC, cBioPortal) under their respective terms of use.
