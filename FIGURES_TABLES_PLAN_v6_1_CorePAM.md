# FIGURES & TABLES PLAN — Core-PAM (Memorial v6.1 | Plano A selado 2026-02-28)

## Regras universais (não-negociáveis)

1. **Sem número de genes no nome** — nunca `PAM24`, `PAM32`, etc. Sempre `CorePAM`.
2. **Bilíngue obrigatório** — cada figura gera dois arquivos: `_EN` (inglês, estilo journal) e `_PT` (português, ABNT).
3. **Resolução mínima** — PNG 600 DPI + PDF vetorial editável (`cairo_pdf`).
4. **Paleta padronizada** — todas as figuras usam `source("scripts/00_colors.R")` antes de qualquer `ggplot`.
5. **Sem hard-code numérico** — todos os números nas figuras/tabelas devem ser lidos de CSV/JSON, nunca embutidos no script.
6. **Validação humana** — após gerar cada figura, enviar caminho ao usuário e aguardar aprovação antes de prosseguir.
7. **Versões anteriores não aprovadas** ficam em disco (não apagar); versões aprovadas substituem no mesmo nome de arquivo.
8. **Diretórios canônicos:**
   - Figuras principais: `results/figures/main/`
   - Figuras suplementares: `results/figures/supp/`
   - Tabelas principais: `results/main/`
   - Tabelas suplementares: `results/supp/`

---

## Figuras principais (manuscrito — 4 figuras)

### Fig 1 — Desenho do estudo + diagrama de congelamento
- **Arquivos:** `results/figures/main/Fig1_StudyDesign_EN.pdf/png` | `Fig1_StudyDesign_PT.pdf/png`
- **Conteúdo:** esquema das 4 coortes (1 treino SCAN-B + 3 validação), setas de fluxo, parâmetros congelados (ΔC=0.010, α=0.5, K=10, fold SHA-256)
- **Fonte de dados:** `01_docs/registry/analysis_freeze.csv`, `01_docs/registry/cohort_manifest.csv`
- **Script gerador:** `logs/plot_study_design.R` (a criar)
- **Status:** ⬜ Pendente

### Fig 2 — Kaplan–Meier por coorte (painéis, sem sobreposição)
- **Arquivos:** `results/figures/main/Fig2_KM_<COHORT>_<ENDPOINT>_EN.pdf/png` | `..._PT.pdf/png`
- **Coortes e endpoints:**
  - TCGA-BRCA: OS
  - METABRIC: DSS (primário); OS (sensibilidade → FigS5)
  - GSE20685: OS
- **Cutoff:** mediana intra-coorte (primário); quartis (sensibilidade → FigS6)
- **Script gerador:** `scripts/07_survival_analysis_<COHORT>.R`
- **Status:** ⬜ Pendente

### Fig 3 — Meta-análise random-effects (HR por 1 SD)
- **Arquivos:** `results/figures/main/Fig3_Meta_Forest_HR1SD_EN.pdf/png` | `..._PT.pdf/png`
- **Conteúdo:** forest plot com HR e IC95% por coorte + diamante resumo; I² e τ²
- **Fonte de dados:** `results/corepam_os/meta_summary.csv`
- **Script gerador:** `scripts/08_meta_survival.R`
- **Status:** ⬜ Pendente

### Fig 4 — Valor incremental (CORE-A vs CORE-A + Core-PAM) + calibração 60m
- **Arquivos:** `results/figures/main/Fig4_DeltaCindex_COREA_EN.pdf/png` | `..._PT.pdf/png`
- **Arquivos (calibração):** `results/figures/main/Fig4_Calibration_60m_EN.pdf/png` | `..._PT.pdf/png`
- **Fonte de dados:** `results/corepam_os/incremental_value.csv`
- **Script gerador:** `scripts/11_incremental_value_and_dca.R`
- **Status:** ⬜ Pendente

---

## Figuras suplementares (auditabilidade e sensibilidade)

### Fig S1 — Fluxograma de inclusão por coorte
- **Arquivos:** `results/figures/supp/FigS1_FlowDiagram_EN.pdf/png` | `..._PT.pdf/png`
- **Conteúdo:** N_raw → N_dedup → N_endpoint → N_final para cada coorte; 4 painéis
- **Fonte de dados:** `results/supp/train_inclusion_log.csv`, logs de cada script 02/03
- **Script gerador:** `logs/plot_flow_diagram.R` (a criar)
- **Status:** ⬜ Pendente

### Fig S2 — Curva Pareto (evidência do platô) ← *movida do main*
- **Arquivos:** `results/figures/supp/FigS2_Pareto_CorePAM_EN.pdf/png` | `..._PT.pdf/png`
- **Conteúdo:** df × C_adj OOF; triângulo = Core-PAM selecionado (df=24); círculo = C_max (df=36); linha tracejada = limiar ΔC
- **Fonte de dados:** `results/corepam/pareto_df_cindex_oof.csv`, `results/corepam/selected_CorePAM_summary.json`
- **Script gerador:** `logs/plot_pareto_curve.R` ✅ Implementado
- **Status:** ✅ **Validada (2026-02-28)**

### Fig S3 — Correlação off-diagonal entre scores (coortes)
- **Arquivos:** `results/figures/supp/FigS3_Correlation_OffDiagonal_EN.pdf/png` | `..._PT.pdf/png`
- **Script gerador:** `scripts/13_qc_correlations.R`
- **Status:** ⬜ Pendente

### Fig S4 — METABRIC PCA forense
- **Arquivos:** `results/figures/supp/FigS4_METABRIC_PCA_Forensics_EN.pdf/png` | `..._PT.pdf/png`
- **Script gerador:** `scripts/14_qc_pca_metabric.R`
- **Status:** ⬜ Pendente

### Fig S5 — METABRIC OS + Fine–Gray (sensibilidade)
- **Arquivos:** `results/figures/supp/FigS5_METABRIC_Sensitivity_EN.pdf/png` | `..._PT.pdf/png`
- **Script gerador:** `scripts/07_survival_analysis_METABRIC.R`
- **Status:** ⬜ Pendente

### Fig S6 — TCGA-BRCA horizonte 24m (sensibilidade)
- **Arquivos:** `results/figures/supp/FigS6_TCGA_24m_Sensitivity_EN.pdf/png` | `..._PT.pdf/png`
- **Script gerador:** `scripts/07_survival_analysis_TCGA.R`
- **Status:** ⬜ Pendente

---

## Tabelas principais

### Table 1 — Características das coortes
- **Arquivo:** `results/main/Table1_CohortCharacteristics.csv`
- **Fonte:** `results/supp/cohort_summary_table.csv`
- **Coortes:** SCAN-B (treino), TCGA-BRCA, METABRIC, GSE20685 (validação)
- **Status:** ⬜ Pendente (dados brutos disponíveis em cohort_summary_table.csv)

### Table 2 — Genes e pesos do Core-PAM
- **Arquivo:** `results/main/Table2_CorePAM_Coefficients.csv`
- **Fonte:** `results/corepam/CorePAM_weights.csv`
- **Conteúdo:** gene, peso (coeficiente Cox elastic-net), direção do risco
- **Status:** ✅ Dados disponíveis — formatação final pendente

### Table 3 — Performance por coorte
- **Arquivo:** `results/main/Table3_Survival_Performance_ByCohort.csv`
- **Conteúdo:** N, eventos, HR (IC95%), C-index (IC95% bootstrap), p-valor
- **Status:** ⬜ Pendente

### Table 4 — Valor incremental
- **Arquivo:** `results/main/Table4_IncrementalValue_DeltaCindex.csv`
- **Fonte:** `results/corepam_os/incremental_value.csv`
- **Status:** ⬜ Pendente

---

## Tabelas suplementares

| ID | Arquivo | Fonte | Status |
|----|---------|-------|--------|
| Table S1 | `results/supp/train_inclusion_log.csv` | Script 02/03 SCAN-B | ✅ Gerado |
| Table S2 | `results/supp/pam50_gene_list_audit_wide.csv` | `logs/generate_audit_files.R` | ✅ Gerado |
| Table S3 | `results/supp/gene_audit_by_cohort.csv` | Scripts 03/04 | ✅ Gerado |
| Table S4 | `results/supp/cohort_summary_table.csv` | Scripts 02 | ✅ Gerado |
| Table S5 | `results/supp/dedup_techrep_scanb.csv` | `logs/generate_audit_files.R` | ✅ Gerado (vazio = sem dups) |
| Table S6 | `results/supp/leakage_check_scanb_vs_gse96058.csv` | Script 01 | ✅ Gerado |
| Table S7 | `results/supp/pam50_coverage_summary.csv` | Script 04 | ✅ Gerado |
| Table S8 | METABRIC competing risks (Fine–Gray) | Script 11 | ⬜ Pendente |

---

## Reconciliação de nomes de arquivo (memorial vs. pipeline)

| Memorial v6.1 §5.4 | Arquivo real gerado | Ação |
|--------------------|--------------------|------|
| `CorePAM_training_card.json` | `selected_CorePAM_summary.json` | Usar nome real; atualizar referência no memorial |
| `results/corepam_os/` | `results/corepam/` (derivação) | `results/corepam_os/` reservado para resultados de sobrevida (scripts 06–11) |
| `figures/main/Fig2_Pareto_*.pdf` | `results/figures/supp/FigS2_Pareto_*_EN/PT.pdf` | Pareto movido para suplemento (Plano A) |

---

## Legenda de status
- ✅ Gerado e validado
- 🔄 Gerado — aguardando validação
- ⬜ Pendente
