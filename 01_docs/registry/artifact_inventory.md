# Inventário de Artefatos — Core-PAM (Memorial v6.1)

Cada artefato listado aqui deve ter SHA-256 registrado em `registry/study_registry.csv`.

## Scripts e seus artefatos esperados

| Script | Artefatos gerados |
|---|---|
| `00_setup.R` | (sourced; sem outputs diretos) |
| `01_download_raw_data.R` | `01_Base_Pura_CorePAM/RAW/<COHORT>/` (expressão + clínica por coorte), `01_docs/registry/data_lake_audit_report.csv`, `results/supp/leakage_check_scanb_vs_gse96058.csv`, append em `registry/study_registry.csv` |
| `02_harmonize_clinical_<COHORT>.R` | `01_Base_Pura_CorePAM/PROCESSED/<COHORT>/clinical_FINAL.parquet`, `01_docs/endpoint_mapping_templates/endpoint_mapping_<COHORT>.csv` |
| `03_expression_preprocess_<COHORT>.R` | `01_Base_Pura_CorePAM/PROCESSED/<COHORT>/expression_genelevel_preZ.parquet`, `results/supp/gene_mapping_audit_<COHORT>.csv` |
| `04_gene_audit_freeze.R` | `results/supp/gene_audit_by_cohort.csv` |
| `05_reduce_pam50_to_corepam_FINAL.R` | `results/corepam/pareto_df_cindex_oof.csv`, `results/corepam/CorePAM_weights.csv`, `results/corepam/CorePAM_model.rds`, `results/corepam/CorePAM_training_card.json`, `results/corepam/artifact_hashes.csv` |
| `06_zscore_and_score_<COHORT>.R` | `01_Base_Pura_CorePAM/PROCESSED/<COHORT>/analysis_ready.parquet`, `results/supp/risk_score_summary_<COHORT>.csv` |
| `07_survival_analysis_<COHORT>.R` | `results/supp/survival_results_<COHORT>.csv`, `figures/main/Fig3_KM_<COHORT>_<ENDPOINT>_CorePAM.pdf/png` |
| `08_meta_survival.R` | `results/main/meta_survival_summary.csv`, `figures/main/Fig4_Meta_Forest_HR_per1SD_CorePAM.pdf/png` |
| `11_incremental_value_and_dca.R` | `results/main/incremental_value_by_cohort.csv`, `figures/main/Fig5_*.pdf/png` |
| `16_qc_text_vs_results_assert.R` | hard-fail se divergência (sem output) |
| `17_render_manuscript_quarto.R` | PDF/HTML final; hash registrado |

## Regras de nomenclatura (Memorial v6.1 §1.5 / Runbook §1.5)

- PROIBIDO: qualquer número de genes no nome de arquivo (`PAM29`, `PAM32`, etc.)
- CORRETO: `CorePAM_weights.csv`, `CorePAM_model.rds`

## Coortes e papéis (congelados)

| Coorte | Papel | Plataforma | Endpoint primário |
|---|---|---|---|
| SCAN-B | TRAIN | RNA-seq | OS |
| GSE96058 | VALIDATION | RNA-seq | OS |
| TCGA-BRCA | VALIDATION | RNA-seq | OS |
| METABRIC | VALIDATION | Microarray | DSS (OS = sensibilidade) |
| GSE20685 | VALIDATION | Microarray | OS |
