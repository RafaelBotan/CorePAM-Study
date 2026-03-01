# Core-PAM — Figures & Tables Index

> **Auto-generated index.** Update after generating any new figure or table.
> Last updated: 2026-03-01
> Branch: `dev`

---

## MAIN FIGURES (Manuscript)

| Label | Title | File (EN) | File (PT) | Script | Status |
|-------|-------|-----------|-----------|--------|--------|
| **Fig 1** | Study design & freeze diagram | `figures/main/Fig1_StudyDesign_EN.pdf` | `_PT.pdf` | `07z_fig1_study_design.R` | ✅ Done |
| **Fig 2a** | KM — TCGA-BRCA (OS, median split) | `figures/main/Fig2_KM_TCGA_BRCA_OS_EN.pdf` | `_PT.pdf` | `07_survival_analysis_TCGA_BRCA.R` | ✅ Done |
| **Fig 2b** | KM — METABRIC (DSS, median split) | `figures/main/Fig2_KM_METABRIC_DSS_EN.pdf` | `_PT.pdf` | `07_survival_analysis_METABRIC.R` | ✅ Done |
| **Fig 2c** | KM — GSE20685 (OS, median split) | `figures/main/Fig2_KM_GSE20685_OS_EN.pdf` | `_PT.pdf` | `07_survival_analysis_GSE20685.R` | ✅ Done |
| **Fig 3** | KM — SCAN-B (OS, training, median split) | `figures/main/Fig3_KM_SCANB_OS_CorePAM.pdf` | — | `07_survival_analysis_SCANB.R` | ✅ Done |
| **Fig 4** | Meta-analysis forest: HR per 1 SD (validation only + FE+RE+I²) | `figures/main/Fig4_Meta_Forest_HR_per1SD_CorePAM.pdf` | `_PT.pdf` | `07z_table_figures_survival.R` / fixed `07x_extra_figures.R` | ✅ Done |
| **Fig 5a** | Incremental value: ΔC-index (CORE-A vs CORE-A+CorePAM) | `figures/main/Fig5_DeltaCindex_COREA_vs_COREAplus_CorePAM.pdf` | — | `11_incremental_value_and_dca.R` | ✅ Done |
| **Fig 5b** | Calibration: predicted vs observed (60m panels) | `figures/main/Fig5_Calibration_60m_Panels_CorePAM.pdf` | — | `11_incremental_value_and_dca.R` | ✅ Done |

> **Note on Fig numbering:** The final journal submission may renumber. The current structure follows the FIGURES_TABLES_PLAN_v6_1.

---

## SUPPLEMENTARY FIGURES

### Block OS (Prognostic — Primary Analysis)

| Label | Title | File (EN) | Script | Status |
|-------|-------|-----------|--------|--------|
| **FigS1** | Pareto curve: C-index OOF vs df | `figures/supp/FigS_Pareto_CorePAM_EN.pdf` | `05_reduce_pam50_to_corepam_FINAL.R` | ✅ Done |
| **FigS2** | CorePAM gene weights (lollipop) | `figures/supp/FigS_Weights_CorePAM_Lollipop.pdf` | `07x_extra_figures.R` | ✅ Done |
| **FigS3** | Off-diagonal score correlations (QC) | `figures/supp/FigS3_Correlation_OffDiagonal.pdf` | `13_qc_correlations_offdiag.R` | ✅ Done |
| **FigS4** | METABRIC PCA forensics | `figures/supp/FigS4_METABRIC_PCA_forensics.pdf` | `14_qc_metabric_pca_forensics.R` | ✅ Done |
| **FigS5** | METABRIC sensitivity: OS + Fine–Gray | `figures/supp/FigS5_METABRIC_Sensitivity_EN.pdf` | `07_survival_analysis_METABRIC.R` | ✅ Done |
| **FigS6** | TCGA-BRCA sensitivity: 24m horizon | `figures/supp/FigS6_TCGA_24m_Sensitivity_EN.pdf` | `07_survival_analysis_TCGA_BRCA.R` | ✅ Done |
| **FigS7** | KM by quartiles — SCAN-B (OS) | `figures/supp/FigS_KM_SCANB_OS_Quartis_CorePAM.pdf` | `07_survival_analysis_SCANB.R` | ✅ Done |
| **FigS8** | KM by quartiles — TCGA-BRCA (OS) | `figures/supp/FigS_KM_TCGA_BRCA_OS_Quartis_CorePAM_EN.pdf` | `07_survival_analysis_TCGA_BRCA.R` | ✅ Done |
| **FigS9** | KM by quartiles — METABRIC (DSS) | `figures/supp/FigS_KM_METABRIC_DSS_Quartis_CorePAM_EN.pdf` | `07_survival_analysis_METABRIC.R` | ✅ Done |
| **FigS10** | KM by quartiles — GSE20685 (OS) | `figures/supp/FigS_KM_GSE20685_OS_Quartis_CorePAM_EN.pdf` | `07_survival_analysis_GSE20685.R` | ✅ Done |
| **FigS11** | C-index by cohort (bootstrap) | `figures/supp/FigS_Cindex_ByCohort.pdf` | `07x_extra_figures.R` | ✅ Done |
| **FigS12** | Forest HR — validation cohorts only | `figures/supp/FigS_Forest_HR_ValidationCohorts.pdf` | `07x_extra_figures.R` | ✅ Done |
| **FigS13** | Gene z-score heatmap — TCGA-BRCA | `figures/supp/FigS_Heatmap_TCGA_BRCA.pdf` | `07x_extra_figures.R` | ✅ Done |
| **FigS14** | Gene z-score heatmap — METABRIC | `figures/supp/FigS_Heatmap_METABRIC.pdf` | `07x_extra_figures.R` | ✅ Done |
| **FigS15** | Gene z-score heatmap — GSE20685 | `figures/supp/FigS_Heatmap_GSE20685.pdf` | `07x_extra_figures.R` | ✅ Done |
| **FigS16** | CorePAM score distribution by cohort (intra-cohort z) | `figures/supp/FigS_CorePAM_RawScore_Distribution_EN.pdf` | `07x_pcr_os_extra_figures.R` | ✅ Done |
| **FigS17** | Score by ER status (SCAN-B + METABRIC, Wilcoxon) | `figures/supp/FigS_ScoreByER_ByCohort.pdf` | `07x_extra_figures.R` (fixed) | ✅ Done |
| **FigS18** | CorePAM vs PAM50-full: scatter (Spearman ρ per cohort) | `figures/supp/FigS_CorePAM_vs_PAM50full_EN.pdf` | `07y_pam50full_comparison.R` | ✅ Done |
| **FigS19** | CorePAM vs PAM50-full: HR forest comparison | `figures/supp/FigS_CorePAM_vs_PAM50full_HR_EN.pdf` | `07y_pam50full_comparison.R` | ✅ Done |
| **FigS20** | Decision Curve Analysis — OS/DSS (4 cohorts, 60m/24m) | `figures/supp/FigS_DCA_OS_EN.pdf` | `11b_dca_corepam.R` | ✅ Done |

### Block pCR (Secondary Analysis) — Directory: `figures/pcr/`

| Label | Title | File (EN) | Script | Status |
|-------|-------|-----------|--------|--------|
| **Fig-pCR1** | pCR forest: OR per 1 SD per cohort + meta pooled (4 primary cohorts) | `figures/pcr/Fig_pCR1_Forest_OR_EN.pdf` | `23_pCR_figures.R` | ✅ Done |
| **FigS-pCR1b** | Extended pCR forest: 4 primary + I-SPY2 (Exploratory) + 2 pooled rows | `figures/pcr/Fig_pCR1_Extended_OR_EN.pdf` | `23b_pCR_ispy2_figures.R` | ✅ Done |
| **FigS-pCR2** | pCR ROC curves by cohort (4 panels) | `figures/pcr/Fig_pCR2_ROC_EN.pdf` | `23_pCR_figures.R` | ✅ Done |
| **FigS-pCR3** | pCR rate by score quartile per cohort | `figures/pcr/Fig_pCR3_QuartileRate_EN.pdf` | `23_pCR_figures.R` | ✅ Done |
| **FigS-pCR4** | Score distribution by pCR status (violin) | `figures/pcr/Fig_pCR4_ScoreDist_EN.pdf` | `23_pCR_figures.R` | ✅ Done |
| **FigS-pCR5** | Precision-recall curves | `figures/pcr/FigS_pCR_PR_EN.pdf` | `23_pCR_figures.R` | ✅ Done |
| **FigS-pCR6** | ER × CorePAM interaction (GSE32646 + ISPY1) | `figures/supp/FigS_pCR_ER_Interaction_EN.pdf` | `21c_pcr_er_interaction.R` | ✅ Done |
| **FigS-pCR7** | Decision Curve Analysis — pCR (4 cohorts) | `figures/supp/FigS_DCA_pCR_EN.pdf` | `11b_dca_corepam.R` | ✅ Done |

---

## MAIN TABLES (Manuscript)

| Label | Title | File | Script | Status |
|-------|-------|------|--------|--------|
| **Table 1** | Cohort characteristics | `results/supp/cohort_summary_table.csv` | `02_harmonize_clinical_*.R` | ✅ Data ready (formatting pending) |
| **Table 2** | CorePAM gene panel (24 genes + weights) | `results/corepam/CorePAM_weights.csv` | `05_reduce_pam50_to_corepam_FINAL.R` | ✅ Done |
| **Table 3** | Survival performance by cohort (HR, C-index, p-value) | `results/main/Table3_Survival_Performance_ByCohort.csv` / `.xlsx` | `07z_table_figures_survival.R` | ✅ Done |
| **Table 4** | Incremental value: ΔC-index per cohort | `results/main/incremental_value_by_cohort.csv` | `11_incremental_value_and_dca.R` | ✅ Done |
| **Table 5** | Meta-analysis summary (validation only: FE+RE+I²) | `results/main/meta_survival_summary_validation_only.csv` | `08_meta_survival.R` / fixed | ✅ Done |

---

## SUPPLEMENTARY TABLES

| Label | Title | File | Script | Status |
|-------|-------|------|--------|--------|
| **TableS1** | PAM50 gene list audit (50 genes vs SCAN-B coverage) | `results/supp/pam50_gene_list_audit.csv` | session diagnostic | ✅ Done |
| **TableS2** | Gene audit by cohort (PAM50 coverage per platform) | `results/supp/gene_audit_by_cohort.csv` | `04_gene_audit_freeze.R` | ✅ Done |
| **TableS3** | Survival results per cohort (detailed) | `results/supp/survival_results_*.csv` | `07_survival_analysis_*.R` | ✅ Done |
| **TableS4** | Calibration data (60m panels) | `results/supp/calibration_*.csv` | `11_incremental_value_and_dca.R` | ✅ Done |
| **TableS5** | METABRIC PCA forensics | `results/supp/metabric_pca_forensics.csv` | `14_qc_metabric_pca_forensics.R` | ✅ Done |
| **TableS6** | QC: off-diagonal score correlations | `results/supp/qc_score_correlations_offdiag.csv` | `13_qc_correlations_offdiag.R` | ✅ Done |
| **TableS7** | PAM50-full comparison (HR + C-index per cohort) | `results/supp/pam50full_comparison.csv` | `07y_pam50full_comparison.R` | ✅ Done |
| **TableS8** | Training cohort inclusion log | `results/supp/train_inclusion_log_SCANB.csv` | session diagnostic | ✅ Done |
| **TableS-pCR1** | pCR results by cohort (OR, AUC, CI) | `results/pcr/pcr_validation_table.csv` | `21_pCR_logistic_analysis.R` | ✅ Done |
| **TableS-pCR2** | pCR meta-analysis (RE/FE, I²) | `results/pcr/pcr_meta_OR.csv` | `22_meta_pCR.R` | ✅ Done |
| **TableS-pCR3** | pCR ER-stratified ORs | `results/pcr/pcr_er_stratified.csv` | `21_pCR_logistic_analysis.R` | ✅ Done |
| **TableS-pCR4** | ER × CorePAM interaction LRT | `results/pcr/pcr_er_interaction.csv` | `21c_pcr_er_interaction.R` | ✅ Done |
| **TableS-pCR5** | pCR gene coverage audit (22/24 per cohort) | `results/pcr/pcr_gene_coverage_audit.csv` | `logs/check_pcr_gene_coverage.R` | ✅ Done |
| **TableS-pCR6** | Sensitivity: 21-gene intersection OR | `results/pcr/pcr_sensitivity_intersection.csv` | `21b_sensitivity_intersection.R` | ✅ Done |
| **TableS-pCR7** | I-SPY2 analysis: univariate + adjusted + control-only | `results/pcr/ispy2_results.csv` | `21_ispy2_analysis.R` | ✅ Done |
| **TableS-pCR8** | Meta-analysis (k=4 without I-SPY2) | `results/pcr/meta_pCR_without_ispy2.csv` | `22b_meta_pCR_with_ispy2.R` | ✅ Done |
| **TableS-pCR9** | Meta-analysis (k=5 with I-SPY2) | `results/pcr/meta_pCR_with_ispy2.csv` | `22b_meta_pCR_with_ispy2.R` | ✅ Done |

---

## KEY RESULTS SUMMARY

| Metric | Value | Source |
|--------|-------|--------|
| CorePAM genes | 24 / 50 PAM50 | `results/corepam/CorePAM_weights.csv` |
| C-index OOF SCAN-B | 0.6699 | `results/corepam/selected_CorePAM_summary.json` |
| HR per 1 SD (SCAN-B) | 1.923 (train) | `results/supp/survival_results_SCANB.csv` |
| HR per 1 SD (TCGA-BRCA OS) | 1.204 | `results/supp/survival_results_TCGA_BRCA.csv` |
| HR per 1 SD (METABRIC DSS) | 1.412 | `results/supp/survival_results_METABRIC.csv` |
| HR per 1 SD (GSE20685 OS) | 1.401 | `results/supp/survival_results_GSE20685.csv` |
| Meta HR RE (val. only) | 1.346 (1.212–1.494) | `results/main/meta_survival_summary_validation_only.csv` |
| Meta I² (val. only) | 42.3% | `results/main/meta_survival_summary_validation_only.csv` |
| pCR meta OR (RE, k=4) | 1.686 (1.385–2.052) I²=0% | `results/pcr/meta_pCR_without_ispy2.csv` |
| pCR meta OR (RE, k=5+ISPY2) | 1.685 (1.498–1.896) I²=0% | `results/pcr/meta_pCR_with_ispy2.csv` |
| I-SPY2 OR univariate | 1.685 (1.454–1.953) p<0.001 AUC=0.648 | `results/pcr/ispy2_results.csv` |
| I-SPY2 OR adjusted (HR+HER2+arm) | 1.485 (1.244–1.771) p<0.001 | `results/pcr/ispy2_results.csv` |
| I-SPY2 OR control-only (Paclitaxel n=179) | 1.443 (0.960–2.168) p=0.078 | `results/pcr/ispy2_results.csv` |
| pCR ER interaction ISPY1 | p=0.0005 | `results/pcr/pcr_er_interaction.csv` |
| CorePAM vs PAM50full (METABRIC C) | 0.638 vs 0.616 | `results/supp/pam50full_comparison.csv` |
| Spearman ρ CorePAM vs PAM50full | 0.77–0.90 | `results/supp/pam50full_comparison.csv` |

---

## PENDING / FLAGGED

| Item | Issue | Priority |
|------|-------|----------|
| FigS-pCR1..5 | All generated in `figures/pcr/` ✅ confirmed | — |
| Table 1 formatting | cohort_summary_table.csv exists; needs journal-ready formatting | Medium |
| Fig 5 PT version | Only EN generated; PT needed for thesis | Low |
| I-SPY2 integration | ✅ Complete: 986 samples, 24/24 genes, OR=1.685 p<0.001, meta with/without ISPY2 | Done |
| pCR calibration plot | Not yet implemented (reliability diagram) | Medium |

---

*Index automatically generated and should be updated whenever a script produces new outputs.*
