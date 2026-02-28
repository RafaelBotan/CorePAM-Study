# Inventário de Artefatos — Core-PAM (Memorial v6.1 | Plano A selado 2026-02-28)

Cada artefato listado aqui deve ter SHA-256 registrado em `registry/study_registry.csv`.

## Scripts e seus artefatos esperados

| Script | Artefatos gerados | Status |
|---|---|---|
| `00_setup.R` | (sourced; sem outputs diretos) | ✅ |
| `00_colors.R` | (sourced; sem outputs diretos — paleta padrão) | ✅ |
| `01_download_raw_data.R` | `01_Base_Pura_CorePAM/RAW/<COHORT>/`, `01_docs/registry/data_lake_audit_report.csv`, `results/supp/leakage_check_scanb_vs_gse96058.csv` | ✅ |
| `02_harmonize_clinical_<COHORT>.R` | `01_Base_Pura_CorePAM/PROCESSED/<COHORT>/clinical_FINAL.parquet`, `01_docs/endpoint_mapping_templates/endpoint_mapping_<COHORT>.csv` | ✅ |
| `03_expression_preprocess_<COHORT>.R` | `01_Base_Pura_CorePAM/PROCESSED/<COHORT>/expression_genelevel_preZ.parquet`, `results/supp/gene_mapping_audit_<COHORT>.csv` | ✅ |
| `04_gene_audit_freeze.R` | `results/supp/gene_audit_by_cohort.csv`, `results/supp/pam50_coverage_summary.csv` | ✅ |
| `05_reduce_pam50_to_corepam_FINAL.R` | `results/corepam/pareto_df_cindex_oof.csv`, `results/corepam/CorePAM_weights.csv`, `results/corepam/CorePAM_model.rds`, `results/corepam/selected_CorePAM_summary.json`¹, `results/corepam/artifact_hashes.csv` | ✅ |
| `logs/generate_audit_files.R` | `results/supp/pam50_gene_list_audit.csv`, `results/supp/pam50_gene_list_audit_wide.csv`, `results/supp/train_inclusion_log.csv`, `results/supp/dedup_techrep_scanb.csv` | ✅ |
| `logs/plot_pareto_curve.R` | `results/figures/supp/FigS2_Pareto_CorePAM_EN.pdf/png`, `results/figures/supp/FigS2_Pareto_CorePAM_PT.pdf/png` | ✅ Validada (2026-02-28) |
| `06_zscore_and_score_<COHORT>.R` | `01_Base_Pura_CorePAM/PROCESSED/<COHORT>/analysis_ready.parquet`, `results/supp/risk_score_summary_<COHORT>.csv` | ⬜ |
| `07_survival_analysis_<COHORT>.R` | `results/corepam_os/survival_results_<COHORT>.csv`, `results/figures/main/Fig2_KM_<COHORT>_<ENDPOINT>_EN.pdf/png`, `..._PT.pdf/png` | ⬜ |
| `08_meta_survival.R` | `results/corepam_os/meta_summary.csv`, `results/figures/main/Fig3_Meta_Forest_HR1SD_EN.pdf/png`, `..._PT.pdf/png` | ⬜ |
| `11_incremental_value_and_dca.R` | `results/corepam_os/incremental_value.csv`, `results/figures/main/Fig4_DeltaCindex_COREA_EN.pdf/png`, `..._PT.pdf/png`, `results/figures/main/Fig4_Calibration_60m_EN.pdf/png`, `..._PT.pdf/png` | ⬜ |
| `16_qc_text_vs_results_assert.R` | hard-fail se divergência (sem output) | ⬜ |
| `17_render_manuscript_quarto.R` | PDF/HTML final; hash registrado | ⬜ |

¹ `selected_CorePAM_summary.json` = equivalente ao `CorePAM_training_card.json` mencionado no Memorial §5.4. Nome real do arquivo gerado pelo script.

---

## Regras de nomenclatura (Memorial v6.1 §1.5 / Runbook §1.5)

- **PROIBIDO:** qualquer número de genes no nome de arquivo (`PAM24`, `PAM29`, `PAM32`, etc.)
- **CORRETO:** `CorePAM_weights.csv`, `CorePAM_model.rds`, `FigS2_Pareto_CorePAM_EN.pdf`
- **BILÍNGUE:** toda figura tem sufixo `_EN` (journal) e `_PT` (ABNT) — 4 arquivos por figura (2 PDF + 2 PNG)
- **RESOLUÇÃO:** PNG 600 DPI + PDF cairo_pdf (vetorial editável)
- **CORES:** obrigatório `source("scripts/00_colors.R")` em todos os scripts de figuras

---

## Coortes e papéis (congelados — Opção A aprovada 2026-02-28)

| Coorte | Papel | Plataforma | Endpoint primário | N | Eventos |
|---|---|---|---|---|---|
| SCAN-B (GSE96058) | **TRAIN** | RNA-seq | OS | 3.069 | 322 |
| TCGA-BRCA | VALIDATION | RNA-seq | OS | 1.072 | 150 |
| METABRIC | VALIDATION | Microarray | DSS (OS = sensibilidade) | 1.980 | 646 DSS |
| GSE20685 | VALIDATION | Microarray | OS | 327 | 83 |

> **GSE96058 não existe como coorte de validação separada.** Todas as 3.069 amostras = SCAN-B treino.

---

## PAM50 canônico (Parker et al. 2009) — congelado

- **50 genes exatos** (Parker JS et al., J Clin Oncol 2009;27:1160–1167, Table 1)
- KRT8 **excluído** — não está na lista canônica
- Auditoria documentada em `results/supp/pam50_gene_list_audit.csv`

---

## Legenda de status
- ✅ Gerado e registrado
- 🔄 Gerado — aguardando renomeação/validação
- ⬜ Pendente
