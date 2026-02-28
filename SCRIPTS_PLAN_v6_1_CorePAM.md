# SCRIPTS_PLAN_v6_1 — Core-PAM (Memorial v6.1)

| Script | Descrição sumária (o que faz / entrega) |
|---|---|
| `scripts/00_setup.R` | Ambiente reprodutível, paths, seeds, helpers (hash/ID/strict read), carrega `analysis_freeze.csv`. |
| `scripts/01_download_raw_data.R` | Download de todas as coortes OS/DSS (SCAN-B via GSE96058 [ALL 3069 amostras = treino], METABRIC, TCGA-BRCA, GSE20685): skip/force, SHA-256, registry. Script único monolítico. |
| `scripts/01_download_raw_data_pcr.R` | Download das coortes pCR/NACT (GSE25066, GSE20194, GSE32646, I-SPY1, I-SPY2). **Executar somente após bloco OS selado.** *(a criar)* |
| `scripts/02_harmonize_clinical_<COHORT>.R` | Harmoniza clínica + endpoints OS/DSS; 4 coortes: SCANB, TCGA_BRCA, METABRIC, GSE20685; METABRIC: DSS correto e parquet final primário. |
| `scripts/03_expression_preprocess_<COHORT>.R` | RNA-seq: counts→TMM/logCPM; microarray: log2; mapeamento HGNC; probes→gene; outputs preZ + audit. 4 coortes: SCANB, TCGA_BRCA, METABRIC, GSE20685. Scripts GSE96058 removidos (Option A, 2026-02-28). |
| `scripts/04_gene_audit_freeze.R` | Auditoria transversal de genes (PAM50 no treino; Core-PAM nas validações). |
| `scripts/05_reduce_pam50_to_corepam_FINAL.R` | Reduz PAM50 para **Core-PAM** (tamanho a determinar): Pareto df×Cadj OOF, não-inferioridade (ΔC), congela `CorePAM_weights.csv`, `CorePAM_model.rds`, `CorePAM_training_card.json`. |
| `scripts/06_zscore_and_score_<COHORT>.R` | Aplica Core-PAM: Z-score intra-coorte, score re-escalonado, direção padronizada; gera `analysis_ready.parquet`. 4 coortes: SCANB, TCGA_BRCA, METABRIC, GSE20685. |
| `scripts/07_survival_analysis_<COHORT>.R` | Cox univ + CORE-A, C-index bootstrap, KM (painéis), calibração 60m; exporta tabelas/figuras. 4 coortes: SCANB, TCGA_BRCA, METABRIC, GSE20685. |
| `scripts/08_meta_survival.R` | Meta-análise random-effects (logHR+SE), I²/τ², forest plot. |
| `scripts/11_incremental_value_and_dca.R` | ΔC-index bootstrap CORE-A vs CORE-A+score; DCA sensibilidade. |
| `scripts/14_qc_metabric_pca_forensics.R` | METABRIC PCA forensics + sensibilidade estratificada. |
| `scripts/15_qc_schema_range_checks.R` | QC estrutural (schema/ranges/IDs/time<=0/missingness). |
| `scripts/16_qc_text_vs_results_assert.R` | Anti hard-code: falha se texto divergir dos CSVs finais. |
| `scripts/17_render_manuscript_quarto.R` | Render final e registro de hash do PDF. |
| `scripts/18_make_submission_bundle.R` | Bundle reexecutável + manifest de hashes. |

## Micro-scripts obrigatórios
| Script | Uso |
|---|---|
| `scripts/07A_preflight_files_strict.R` | I/O strict: warning=erro + SHA-256. |
| `scripts/07B_find_metabric_clinical_fallback.R` | Localiza fallback íntegro para clínica METABRIC. |
| `scripts/07D_validate_one_cohort_corepam.R` | Teste mínimo por coorte (join, cobertura, score, HR). |

## Scripts removidos (Option A, aprovado 2026-02-28)
Os scripts abaixo foram removidos porque GSE96058 não existe mais como coorte de validação separada. SCAN-B = TODAS as 3069 amostras do GEO GSE96058 (treinamento).
- `scripts/03_expression_preprocess_GSE96058.R` — removido
- `scripts/06_zscore_and_score_GSE96058.R` — removido
- `scripts/07_survival_analysis_GSE96058.R` — removido
