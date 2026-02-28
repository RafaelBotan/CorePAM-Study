# SCRIPTS_PLAN_v6_1 — Core-PAM (Memorial v6.1)

| Script | Descrição sumária (o que faz / entrega) |
|---|---|
| `scripts/00_setup.R` | Ambiente reprodutível, paths, seeds, helpers (hash/ID/strict read), carrega `analysis_freeze.csv`. |
| `scripts/01_ingest_raw_<COHORT>.R` | Ingestão determinística + RAW imutável + hashes no registry. |
| `scripts/02_harmonize_clinical_<COHORT>.R` | Harmoniza clínica + endpoints OS/DSS; METABRIC: DSS correto e parquet final primário. |
| `scripts/03_expression_preprocess_<COHORT>.R` | RNA-seq: counts→TMM/logCPM; microarray: log2; mapeamento HGNC; probes→gene; outputs preZ + audit. |
| `scripts/04_gene_audit_freeze.R` | Auditoria transversal de genes (PAM50 no treino; Core-PAM nas validações). |
| `scripts/05_reduce_pam50_to_corepam_FINAL.R` | Reduz PAM50 para **Core-PAM** (tamanho a determinar): Pareto df×Cadj OOF, não-inferioridade (ΔC), congela `CorePAM_weights.csv`, `CorePAM_model.rds`, `CorePAM_training_card.json`. |
| `scripts/06_zscore_and_score_<COHORT>.R` | Aplica Core-PAM: Z-score intra-coorte, score re-escalonado, direção padronizada; gera `analysis_ready.parquet`. |
| `scripts/07_survival_analysis_<COHORT>.R` | Cox univ + CORE-A, C-index bootstrap, KM (painéis), calibração 60m; exporta tabelas/figuras. |
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
