# FIGURES & TABLES PLAN — Core-PAM (Memorial v6.1)

Regra: **nenhum arquivo deve conter o número de genes no nome**. O tamanho final é reportado dentro das tabelas/resultados.

---

## Figuras principais (manuscrito)

### Fig 1 — Study design + freeze diagram (Core-PAM)
- `figures/main/Fig1_StudyDesign_Freeze_CorePAM.pdf/png`
- Fonte: `cohort_manifest.csv`, `analysis_freeze.csv`.

### Fig 2 — Pareto (df vs Cadj OOF) e escolha do Core-PAM
- `figures/main/Fig2_Pareto_CorePAM_Cadj_OOF.pdf/png`
- Tabela fonte: `results/corepam/pareto_df_cindex_oof.csv`
- Script: `05_reduce_pam50_to_corepam_FINAL.R`

### Fig 3 — Kaplan–Meier por coorte (painéis, sem sobreposição)
- `figures/main/Fig3_KM_<COHORT>_<ENDPOINT>_CorePAM.pdf/png`
- Script: `07_survival_analysis_<COHORT>.R`

### Fig 4 — Meta-análise random-effects (HR por 1 SD)
- `figures/main/Fig4_Meta_Forest_HR_per1SD_CorePAM.pdf/png`
- Script: `08_meta_survival.R`

### Fig 5 — Valor incremental (CORE-A vs CORE-A+Core-PAM) + calibração 60m
- `figures/main/Fig5_DeltaCindex_COREA_vs_COREAplus_CorePAM.pdf/png`
- `figures/main/Fig5_Calibration_60m_Panels_CorePAM.pdf/png`
- Script: `11_incremental_value_and_dca.R`

---

## Suplementos (auditabilidade)

- Fig S1 Fluxo por coorte: `figures/supp/FigS1_Flow_PerCohort_CorePAM.pdf`
- Fig S3 Correlação off-diagonal: `figures/supp/FigS3_Correlation_OffDiagonal.pdf`
- Fig S4 METABRIC PCA forensics: `figures/supp/FigS4_METABRIC_PCA_Forensics.pdf`
- Fig S5 METABRIC OS + Fine–Gray (sensibilidade): `figures/supp/FigS5_METABRIC_Sensitivity.pdf`
- Fig S6 TCGA 24m (sensibilidade): `figures/supp/FigS6_TCGA_24m_Sensitivity.pdf`

---

## Tabelas principais

- Table 1 Cohort characteristics: `results/main/Table1_CohortCharacteristics.csv`
- Table 2 Core-PAM genes + pesos: `results/main/Table2_CorePAM_Coefficients.csv` (fonte `results/corepam/CorePAM_weights.csv`)
- Table 3 Performance por coorte: `results/main/Table3_Survival_Performance_ByCohort.csv`
- Table 4 Incremental value: `results/main/Table4_IncrementalValue_DeltaCindex.csv`

---

## Tabelas suplementares (exemplos)

- Table S2 `expression_input_spec.csv`
- Table S3 `gene_audit_by_cohort.csv`
- Table S6 Leakage check
- Table S8 METABRIC competing risks
