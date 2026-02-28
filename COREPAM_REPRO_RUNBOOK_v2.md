---
title: "Core-PAM — Reproducibility Runbook & Repo Operations (v2)"
date: "2026-02-28"
repo: "RafaelBotan/Phd-Genomic (Private)"
project: "Core-PAM (PAM50 → Core-PAM reduced panel; gene-count is derived)"
---

# 1) GitHub repository: what exists and how to “push” correctly

## 1.1 Repository identity (current)
- **Owner:** `RafaelBotan`
- **Repository:** `Phd-Genomic`
- **Visibility:** Private
- **Default branch:** `main`
- **Working branch:** `dev` (create now)

## 1.2 Remote URLs (choose one)
**SSH**
```bash
git@github.com:RafaelBotan/Phd-Genomic.git
```
**HTTPS**
```bash
https://github.com/RafaelBotan/Phd-Genomic.git
```

## 1.3 First-time local setup (recommended, clean)
### Option A — you do not have a local folder yet
```bash
cd "Y:\GIT Projeto"
git clone git@github.com:RafaelBotan/Phd-Genomic.git
cd Phd-Genomic
git checkout -b dev
git push -u origin dev
```

### Option B — you already have a local project folder (no git yet)
Run inside your local folder (the one containing `scripts/`, `01_docs/`, etc.):
```bash
git init
git branch -M main
git remote add origin git@github.com:RafaelBotan/Phd-Genomic.git
git add .
git commit -m "init: Core-PAM reproducibility scaffold"
git push -u origin main
git checkout -b dev
git push -u origin dev
```

## 1.4 Uploading the generated documentation files (minimum)
Place these files in the repo (suggested locations):
- `01_docs/`
  - `Memorial_v6_1_CorePAM.md`
  - `CHECKLIST_MASTER_v6_1_CorePAM.md`
  - `FIGURES_TABLES_PLAN_v6_1_CorePAM.md`
  - `SCRIPTS_PLAN_v6_1_CorePAM.md`
  - `COREPAM_REPRO_RUNBOOK_v2.md` *(this file)*
  - `rastreabilidade_metodologica_corepam.qmd` *(from template)*

Then commit:
```bash
git add 01_docs/*
git commit -m "docs: v6.1 Core-PAM memorial + checklists + runbook"
git push
```

## 1.5 Core rule: gene-count must never be encoded in filenames
- Prohibited: `PAM29_weights.csv`, `PAM32_model.rds`, etc. *(use neutral names)*
- Allowed: `CorePAM_weights.csv`, `CorePAM_model.rds`, `selected_CorePAM_summary.json`

# 2) Local data lake (raw/processed separation)

## 2.1 Required separation
- **Git repo**: code + manifests + small tables + documentation.
- **Data lake (outside Git)**: raw downloads and heavy processed objects.

## 2.2 Suggested paths (Windows; adjust once and freeze in 00_setup.R)
- Repo (code): `Y:\GIT Projeto\Phd-Genomic\`
- Data lake (raw+processed): `Y:\GIT Projeto\Phd-Genomic\01_Base_Pura_CorePAM\`

## 2.3 .gitignore (mandatory)
At minimum:
- `data/raw/**`
- `data/processed/**`
- `**/*.parquet`
- `**/*.rds`
- `**/*.zip`
- `**/.DS_Store`

# 3) Cohorts: what must exist locally, what is optional

## 3.1 Mandatory (OS/DSS focus)
- SCAN-B (TRAIN; RNA-seq; GEO accession GSE96058, ALL 3069 samples = training; no internal split; confirmed 2026-02-28)
- METABRIC (VALIDATION; microarray; OS + cause-of-death field for DSS mapping)
- TCGA-BRCA (VALIDATION; RNA-seq)
- GSE20685 (VALIDATION; microarray; OS)

## 3.2 Optional later (pCR; separate block)
- GSE25066, GSE20194, GSE32646, I-SPY1, I-SPY2

## 3.3 Accessory registries (must be in repo)
Inside `01_docs/registry/`:
- `cohort_manifest.csv`
- `analysis_freeze.csv`
- `expression_input_spec.csv`
- `pcr_dictionary_registry.csv` *(even if pCR deferred, keep registry ready)*
- `study_registry.csv` *(execution registry; append-only)*
- `artifact_inventory.md` *(inventory of deliverables and meaning)*

# 4) The QMD (living manuscript) — mandatory routine

## 4.1 Required file
- `01_docs/rastreabilidade_metodologica_corepam.qmd`

## 4.2 Non-negotiable rule
- No hard-coded Ns, HRs, C-index, follow-up, or gene-count.
- The QMD must read:
  - summary JSONs (selected Core-PAM)
  - final tables CSV
  - figure paths

## 4.3 When to render
Render after every sealed step (see section 9).

# 5) Core-PAM derivation: minimal gene-count while preserving function

## 5.1 Goal
Reduce PAM50 → Core-PAM with the **smallest feasible gene-count** while maintaining non-inferior performance.

## 5.2 Primary derivation metric
- **OOF Harrell C-index** in the TRAIN cohort using fixed stratified folds.

## 5.3 Non-inferiority selection rule (frozen)
- Compute **Cmax** across candidate dfs.
- Choose smallest df with **C-index ≥ (Cmax − ΔC)**.
- ΔC is fixed in `analysis_freeze.csv`.

## 5.4 Output contract (must exist after derivation)
- `results/pamxx/pareto_df_cindex_oof.csv`
- `results/pamxx/selected_CorePAM_summary.json`
- `results/pamxx/CorePAM_weights.csv`
- `results/pamxx/CorePAM_model.rds`
- `results/pamxx/artifact_hashes.csv`

# 6) OS/DSS validation: what must be produced per cohort

## 6.1 Per-cohort deliverables (minimum)
For each OS/DSS cohort:
- cohort passport (N, events, follow-up median via Reverse KM)
- gene coverage report (present/missing; denominator used for effective score)
- score direction audit (ensure higher score = higher risk)
- Cox univariate: HR per 1 SD
- CORE-A + score (Age + ER) where available
- C-index bootstrap CI
- KM panels (median; sensitivity quartiles)

## 6.2 Meta-analysis deliverables
- random-effects meta-analysis of log(HR) across cohorts
- I² and τ²
- leave-one-out analysis

# 7) Figure and table exports (project + journal)

## 7.1 Project exports (always)
- Every figure: **PDF (vector)** + **PNG (high-res)**
- English outputs: `06_plots/artigo/`
- Portuguese outputs: `06_plots/tese/`

## 7.2 Breast Cancer Research constraints (must comply)
- Main manuscript file formats accepted: DOC/DOCX or RTF; double line spacing; line and page numbering; avoid page breaks.
- Research article abstract: max 350 words; structured (Background/Methods/Results/Conclusions) and no references in abstract; 3–10 keywords; include Declarations section with required subheadings.
- Figures: numbered in order; multi-panel figures must be a single composite file; titles ≤15 words; legends ≤300 words in the manuscript (not inside figure); include figure keys inside the graphic; crop whitespace; each figure file ≤10 MB.
- Accepted figure formats: EPS, PDF, Word, PowerPoint, TIFF, JPEG, PNG, BMP, CDX; prefer PDF for vector.
- PDF sizing: 85 mm half-width, 170 mm full-width; max height 225 mm (figure+legend); ~300 dpi at final size; embed fonts; lines >0.25 pt.
- Tables: numbered and cited in order; avoid color/shading; titles ≤15 words and legends ≤300 words; wide/large tables can be additional files; additional tabular data can be XLS or CSV.
- Additional files: cite sequentially; max 20 MB each; avoid patient consent forms; do not include individual participant details; prefer repository deposition over personal websites.
- Data availability: include 'Availability of data and materials' statement; cite public datasets with persistent identifiers; describe conditions for restricted data access.

# 8) Peer-review watchlist (from predecessor study; track continuously)
Each item below must have an explicit artifact proving compliance.
- Include at least one additional RNA-seq cohort for validation (RNA-seq is current standard).
- Share full source code/workflow (public link or controlled-access repo).
- Methods must specify exactly which expression matrices were downloaded (counts/TPM/intensity) and all transformations (normalization/log/batch/Z-score).
- Define thresholds objectively (e.g., 'weak expression') and define how cross-platform reproducibility was assessed.
- Ensure absolute consistency between manuscript text and figure annotations (HR, follow-up, N).
- Avoid unreadable overplotted multi-cohort curves; use panels/facets.
- Report intrinsic subtype frequencies per cohort and compare full PAM50 vs reduced Core-derived calls where possible.
- Correlation summaries must exclude diagonal autocorrelation (self-correlation=1).
- Investigate unexplained PCA clusters in METABRIC; test association with technical/clinical variables and discuss.
- Keep dataset counts consistent throughout the manuscript.
- Do not pool distinct endpoints (OS vs DRFS) without adjustment; keep endpoints separated.
- If TCGA C-index is low, do sensitivity analyses (e.g., fixed horizon) rather than speculation.
- Include clinical covariates to demonstrate incremental value beyond simple clinicopathological factors.

# 9) The “sealed step” ritual (Definition of Done)
A script/task is only DONE when all items are satisfied:

1. Run script in a clean R session
2. Outputs generated exactly as declared (no silent failures)
3. SHA-256 hashes computed for each output
4. Append rows to `registry/study_registry.csv` (append-only)
5. Update documentation:
   - `01_docs/rastreabilidade_metodologica_corepam.qmd`
   - `01_docs/CHECKLIST_MASTER_*`
   - `01_docs/FIGURES_TABLES_PLAN_*` (if deliverables changed)
   - `01_docs/SCRIPTS_PLAN_*` (if scripts changed)
6. Render QMD (HTML+PDF)
7. Commit + push to `dev`
8. If milestone: merge to `main` and tag freeze (`freeze-corepam-YYYYMMDD`)

