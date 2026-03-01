# CLAUDE.md — Core-PAM Project

## Visão geral

Projeto de doutorado de Rafael Botan. Redução do painel PAM50 para um núcleo mínimo (**Core-PAM**, tamanho derivado — não pré-fixado) com validação prognóstica multicoorte em câncer de mama. Estudo retrospectivo 100% público, multi-coorte, multi-plataforma (RNA-seq + microarray).

## Regras absolutas (nunca violar)

- **NUNCA** usar número de genes em nomes de arquivo ou variável (`PAM29`, `PAM32` etc.) — usar sempre `CorePAM`
- **NUNCA** fazer pooling de expressão entre coortes
- **NUNCA** fazer batch correction inter-coorte
- **NUNCA** refitar coeficientes nas coortes de validação
- **Z-score sempre intra-coorte** (ou intra-subcoorte em multi-plataforma)
- **Warning de leitura = erro** (`options(warn=2)` ativo via `00_setup.R`)
- **Tempo sempre em meses** (dias / 30.4375)
- **time ≤ 0 → DROP** (reportar contagem)
- **METABRIC endpoint primário = DSS** (não OS)
- **Leakage-proof obrigatório** entre SCAN-B e coortes de validação externas (TCGA-BRCA, METABRIC, GSE20685)
- Todo artefato recebe **SHA-256** e é registrado em `registry/study_registry.csv`
- Manuscrito **sem números hard-coded** (lê de CSVs/JSONs)

## Parâmetros congelados (analysis_freeze.csv)

| Parâmetro | Valor |
|---|---|
| delta_c | 0.010 |
| alpha (elastic-net) | 0.5 |
| k_folds | 10 |
| seed_folds | 42 |
| min_genes_fraction | 0.80 |
| time_unit_divisor | 30.4375 |

## Estrutura de diretórios

```
Y:/Phd-Genomic-claude/          ← ROOT_REPO (working directory)
├── scripts/                    ← todos os scripts R numerados
├── 01_docs/
│   ├── registry/               ← analysis_freeze.csv, manifests, audit reports
│   └── endpoint_mapping_templates/
├── 01_Base_Pura_CorePAM/       ← DATA_LAKE (gitignored)
│   ├── RAW/<COHORT>/           ← dados brutos imutáveis
│   └── PROCESSED/<COHORT>/     ← expression_genelevel_preZ.parquet, clinical_FINAL.parquet, analysis_ready.parquet
├── results/
│   ├── corepam/                ← CorePAM_weights.csv, CorePAM_model.rds, pareto_df_cindex_oof.csv, CorePAM_training_card.json
│   ├── corepam_os/             ← resultados de sobrevida por coorte
│   ├── main/                   ← tabelas principais do manuscrito
│   └── supp/                   ← tabelas suplementares, auditorias, QC
├── figures/
│   ├── main/                   ← Fig1-5 (PDF + PNG)
│   └── supp/                   ← FigS1-S8
├── 06_plots/
│   ├── artigo/                 ← figuras em inglês
│   └── tese/                   ← figuras em português
├── registry/
│   └── study_registry.csv      ← append-only; todos os hashes (gitignored)
└── Memorial_v6_1_CorePAM.md    ← fonte da verdade metodológica
```

## Coortes e papéis (congelados)

| Coorte | Papel | Plataforma | Endpoint |
|---|---|---|---|
| SCAN-B | TRAIN | RNA-seq | OS |
| TCGA-BRCA | VALIDATION | RNA-seq | OS (sensib. 24m) |
| METABRIC | VALIDATION | Microarray Illumina | DSS (OS = sensib.) |
| GSE20685 | VALIDATION | Microarray Affymetrix | OS |

**SCAN-B usa o acesso GEO GSE96058 (TODAS as 3069 amostras = treinamento). Não há split interno.** O pData do GEO não possui coluna de separação training/validation. Decisão aprovada em 2026-02-28 (Option A): usar todas as amostras como treino.

## Pipeline de execução (ordem)

```
00_setup.R                          ← sourced por todos; nunca executar diretamente
01_download_raw_data.R              ← GEO + GDC + cBioPortal + SHA-256
02_harmonize_clinical_<COHORT>.R    ← clínica + endpoints + endpoint_mapping (4 coortes: SCANB, TCGA_BRCA, METABRIC, GSE20685; SCANB cobre dataset GSE96058 completo)
03_expression_preprocess_<COHORT>.R ← TMM/logCPM (RNA-seq) ou log2 as-is (array) + HGNC
04_gene_audit_freeze.R              ← auditoria transversal PAM50 (GO/NO-GO)
05_reduce_pam50_to_corepam_FINAL.R  ← derivação Core-PAM (FREEZE)
06_zscore_and_score_<COHORT>.R      ← Z-score + score re-escalonado + direção
07A_preflight_files_strict.R        ← verificação estrita de inputs (rodar antes do 07)
07D_validate_one_cohort_corepam.R   ← smoke test por coorte
07_survival_analysis_<COHORT>.R     ← Cox, C-index, KM, CORE-A
08_meta_survival.R                  ← meta-análise random-effects
11_incremental_value_and_dca.R      ← ΔC-index + calibração
13_qc_correlations_offdiag.R        ← QC correlações off-diagonal
14_qc_metabric_pca_forensics.R      ← PCA forense METABRIC
15_qc_schema_range_checks.R         ← QC estrutural hard-fail
16_qc_text_vs_results_assert.R      ← anti hard-code hard-fail
17_render_manuscript_quarto.R       ← render QMD final
18_make_submission_bundle.R         ← bundle reexecutável

# pCR block (NACT — análise secundária, separada do bloco OS)
19_download_pCR_raw_data.R          ← download GEO: GSE25066, GSE20194, GSE32646, GSE22226
20_prepare_pCR_<COHORT>.R           ← harmonize + preprocess + Z-score + CorePAM score (4 scripts)
21_pCR_logistic_analysis.R          ← glm(pcr ~ score_z) + AUC + bootstrap OR CI por coorte
22_meta_pCR.R                       ← meta-análise RE/FE de log(OR)
23_pCR_figures.R                    ← Fig6 forest + FigS9 ROC + FigS10 distribuição
```

## Helpers disponíveis (carregados via 00_setup.R)

```r
sha256_file(path)                   # SHA-256 de um arquivo
normalize_id(x)                     # trimws + toupper
strict_csv(path)                    # read_csv com verificação
strict_parquet(path)                # read_parquet com verificação
strict_rds(path)                    # readRDS com verificação
registry_append(cohort, file_type,  # append ao study_registry.csv
                file_path, sha256,
                status, script, size_mb)
raw_cohort(cohort)                  # path para RAW/<cohort>
proc_cohort(cohort)                 # path para PROCESSED/<cohort>
FREEZE$delta_c / $alpha / $k_folds / $seed_folds / ...
PATHS$results$corepam / $main / $supp / ...
```

## Score (fórmula única para todas as coortes)

```
G_present = genes Core-PAM ∩ genes da coorte
z_i       = Z-score por gene intra-coorte
score     = Σ(wᵢ · zᵢ) / Σ|wᵢ|   (apenas genes presentes)
score_z   = scale(score)            (para HR por 1 SD)
```

Se HR(score) < 1 → inverter sinal; registrar `score_direction = -1`.

## Artefatos principais do freeze (results/corepam/)

- `CorePAM_weights.csv` — genes + pesos do painel derivado
- `CorePAM_model.rds` — glmnet.fit no λ escolhido (full-data path)
- `CorePAM_training_card.json` — parâmetros, λ, df, Cadj, hashes
- `pareto_df_cindex_oof.csv` — fronteira Pareto df × Cadj OOF

## Coding conventions (mandatory for all new scripts)

- **Language:** All code comments and documentation in **English**
- **Variables:** English names following tidyverse snake_case conventions
- **Skip/force pattern:** Every script that writes output files must check:
  ```r
  FORCE <- as.logical(Sys.getenv("FORCE_RERUN", "FALSE"))
  if (!FORCE && file.exists(out_path)) {
    message("[SCRIPT] Output already exists, skipping. Set FORCE_RERUN=TRUE to rerun.")
    quit(save = "no", status = 0)
  }
  ```
- **Error logging:** Before fixing any script error, append entry to `ERRORS_LOG.md`
  with: date, script, error text, root cause, solution.
  (Note: ERRORS_LOG.md, SESSION_NOTES.md, METHODS_DRAFT.md, RESULTS_SUMMARY.md are
   local-only files — they are in .gitignore and should NOT be committed to git.)
- **Results summary:** After each analysis, append 3-line summary to `RESULTS_SUMMARY.md`.
- **Methods draft:** After finishing each script, append academic Portuguese paragraph
  to `METHODS_DRAFT.md` (kept local only).
- **Session notes:** At end of session update `SESSION_NOTES.md` (local only; not committed).

## Support files

| File | Purpose |
|---|---|
| `ERRORS_LOG.md` | Log of all pipeline errors + solutions (append before fixing) |
| `RESULTS_SUMMARY.md` | 3-line summaries of each analysis result |
| `METHODS_DRAFT.md` | Academic Portuguese methods paragraphs for thesis |
| `SESSION_NOTES.md` | End-of-session progress notes (auto-committed) |
| `scripts/00_validate_all.R` | Syntax check + output presence check for all scripts |

## Script status

| Script | Status |
|---|---|
| 00_setup.R | Complete and tested |
| 00_validate_all.R | Complete (syntax + output presence check) |
| 01_download_raw_data.R | Complete (requires data) |
| 02_harmonize_clinical_*.R | Complete (4 cohorts: SCANB, TCGA_BRCA, METABRIC, GSE20685) |
| 03_expression_preprocess_*.R | Complete (4 cohorts: SCANB, TCGA_BRCA, METABRIC, GSE20685) |
| 04_gene_audit_freeze.R | Complete |
| 05_reduce_pam50_to_corepam_FINAL.R | Complete |
| 06_zscore_and_score_*.R | Complete (4 cohorts: SCANB, TCGA_BRCA, METABRIC, GSE20685) |
| 07A_preflight_files_strict.R | Complete |
| 07D_validate_one_cohort_corepam.R | Complete |
| 07_survival_analysis_*.R | Complete (4 cohorts: SCANB, TCGA_BRCA, METABRIC, GSE20685) |
| 08_meta_survival.R | Complete |
| 11_incremental_value_and_dca.R | Complete |
| 13_qc_correlations_offdiag.R | Complete |
| 14_qc_metabric_pca_forensics.R | Complete |
| 15_qc_schema_range_checks.R | Complete |
| 16_qc_text_vs_results_assert.R | Complete |
| 17_render_manuscript_quarto.R | Complete |
| 18_make_submission_bundle.R | Complete |
| 19_download_pCR_raw_data.R | Complete (pCR block — requires internet) |
| 20_prepare_pCR_*.R | Complete (4 cohorts: GSE25066, GSE20194, GSE32646, ISPY1) |
| 20_utils_pcr_extract.R | Complete (shared pCR helper) |
| 21_pCR_logistic_analysis.R | Complete |
| 22_meta_pCR.R | Complete |
| 23_pCR_figures.R | Complete (Fig6, FigS9, FigS10) |

## Git

- **Repo:** https://github.com/RafaelBotan/Phd-Genomic.git (privado)
- **Branch de trabalho:** `dev`
- **Milestones:** merge para `main` + tag `freeze-corepam-YYYYMMDD`
- **R:** `C:/Program Files/R/R-4.5.2/bin/Rscript.exe`
- **Testar sintaxe:** `Rscript --no-save -e "parse(file='scripts/X.R'); cat('OK')"`

## Referências dos documentos de governança

- `Memorial_v6_1_CorePAM.md` — protocolo completo (fonte da verdade)
- `CHECKLIST_MASTER_v6_1_CorePAM.md` — checklist operacional
- `FIGURES_TABLES_PLAN_v6_1_CorePAM.md` — plano de figuras e tabelas
- `SCRIPTS_PLAN_v6_1_CorePAM.md` — scripts e artefatos
- `COREPAM_REPRO_RUNBOOK_v2.md` — runbook Git, data lake, sealed step ritual
- `01_docs/registry/analysis_freeze.csv` — parâmetros congelados
- `01_docs/registry/artifact_inventory.md` — inventário de artefatos

## Autonomy rules

- Do not ask for confirmation between scripts in the pipeline
- When multiple approaches exist, choose the most robust one and proceed
- Always continue to the next step automatically after completing a task
- Only stop if a critical error cannot be resolved without user input
- Commit and push to `dev` automatically after each completed pipeline step
- Never pause to ask "shall I continue?" — just continue

## Sealed Step Ritual (Definition of Done para cada script)

1. Rodar em sessão R limpa
2. Outputs gerados conforme declarado
3. SHA-256 de cada output
4. Append ao `registry/study_registry.csv`
5. Atualizar QMD + checklists
6. Render QMD
7. Commit + push para `dev`
8. Se milestone: merge `main` + tag freeze
