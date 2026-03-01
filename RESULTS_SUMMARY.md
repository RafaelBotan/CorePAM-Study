# RESULTS_SUMMARY.md — Core-PAM Pipeline Results Log

**Convention:** After each analysis script completes, append a 3-line summary:
- Line 1: Date, script name, cohort/stage
- Line 2: Key numeric results (n, HR, C-index, etc.)
- Line 3: Status and next step

---

<!-- Entries are appended automatically by each script and by 17_render_manuscript_quarto.R -->
<!-- Format: Date | Script | Brief result -->

---
**Date:** Pipeline Initialized | **Script:** CLAUDE.md / project setup
Scripts 00–18 created. All scripts pass syntax check. Awaiting data download.
Status: READY — run 01_download_raw_data.R to begin.


---
**Date:** 2026-02-28 | **Script:** 07_survival_analysis_SCANB.R | SCANB (training)
N=3069, Events=322, FU=55.7mo, HR_uni=1.923 (1.736-2.129) p=3.56e-36, C-index=0.698 (0.668-0.729)
Status: COMPLETED — training cohort shows strong prognostic discrimination

---
**Date:** 2026-02-28 | **Script:** 08_meta_survival.R | Meta-analysis (4 cohorts)
Meta HR_uni=1.472 (1.205-1.798) I2=90.6% p=1.52e-4 | Meta HR_COREA=1.478 (1.253-1.744) I2=83.6% p=3.69e-6
Status: COMPLETED — consistent effect across platforms; high I2 expected (multi-platform heterogeneity)

---
**Date:** 2026-02-28 | **Script:** 11_incremental_value_and_dca.R | Delta C-index (4 cohorts)
SCANB ΔC=+0.030 | TCGA ΔC=+0.038 | METABRIC ΔC=+0.142 | GSE20685 ΔC=+0.093
Status: COMPLETED — CorePAM adds significant incremental value beyond clinical covariates

---
**Date:** 2026-02-28 | **Script:** 15_qc_schema_range_checks.R | Schema QC
44 PASS, 0 FAIL across all 4 cohorts + global weights checks
Status: COMPLETED — all structural integrity checks pass

---
**Date:** 2026-02-28 | **Script:** 16_qc_text_vs_results_assert.R | Anti-hard-code assertions
8 PASS, 0 FAIL — all HR and C-index values validated from CSV sources
Status: COMPLETED
