# SESSION_NOTES.md — Core-PAM Pipeline Session Log

---

## Session: 2026-02-28 (Continuação — Pipeline 07-16)

### What was done
- Retomada do pipeline após sessão anterior (scripts 07x, 07z já existentes)
- **07_survival_analysis_SCANB.R** — corrigido e executado:
  - Bugs: `time_col = "os_time"` → `"os_time_months"`, função `bootstrap_cindex` sem Cadj
  - Resultado: SCANB HR=1.923, C-index=0.698 (0.668-0.729), N=3069, FU=55.7mo
- **08_meta_survival.R** — corrigido e executado:
  - Bug: `geom_errorbarh()` deprecated → `geom_errorbar(orientation="y")`
  - Resultado: Meta HR=1.472 (1.205-1.798), I2=90.6%, 4 coortes
- **11_incremental_value_and_dca.R** — corrigido e executado:
  - Bugs: ENDPOINT_MAP sem `_months`, fallback sem `_months`, `concordance(~matrix)` returns vector, CORE-A NA filter ausente
  - Resultado: ΔC positivo em todas coortes (SCANB=+0.030, TCGA=+0.038, METABRIC=+0.142, GSE20685=+0.093)
- **13_qc_correlations_offdiag.R** — executado (correlações off-diagonal)
- **14_qc_metabric_pca_forensics.R** — executado (1 outlier: MB-5135)
- **15_qc_schema_range_checks.R** — corrigido e executado: 44 PASS, 0 FAIL
  - Bugs: ENDPOINT_MAP, `sample_id` vs `patient_id`, `dss_time` vs `dss_time_months`
- **16_qc_text_vs_results_assert.R** — executado após reinstalar lubridate: 8 PASS, 0 FAIL
- FigS5 regenerada (METABRIC Sensitivity) — foi relatada como "saiu pela metade" em versão anterior
- Agente em background escrevendo interpretações das figuras em `figures/interpretations/`

### New requirement (user request)
- Ao final de cada figura validada: criar arquivo `.md` com o mesmo nome da figura
  em `figures/interpretations/`, ultra didático, explicando "com amor" todos os achados

### What works
- Pipeline completo executado: scripts 07 (SCANB), 07A, 07D, 07 (3 validações), 07x, 07z,
  08, 11, 13, 14, 15, 16
- Novos artefatos: Fig3 (SCANB KM), Fig4 (Meta-análise), Fig5 (ΔC + Calibração)
- Meta-análise: 4 coortes combinadas, efeito consistente
- Valor incremental: CorePAM supera CORE-A clínico em todas coortes
- QC: todos os checks passam
- Commit c2f59e8 com 61 arquivos, push para dev

### Números chave do pipeline
| Coorte    | Endpoint | N    | HR (95%CI)            | C-index (95%CI)       |
|-----------|----------|------|-----------------------|-----------------------|
| SCANB     | OS       | 3069 | 1.923 (1.736-2.129)   | 0.698 (0.668-0.729)   |
| TCGA_BRCA | OS       | 1072 | 1.204 (1.036-1.400)   | 0.624 (0.572-0.678)   |
| METABRIC  | DSS      | 1978 | 1.412 (1.310-1.523)   | 0.638 (0.615-0.659)   |
| GSE20685  | OS       | 327  | 1.401 (1.118-1.755)   | 0.623 (0.562-0.683)   |
| **Meta**  | —        | —    | **1.472 (1.205-1.798)** | I2=90.6%            |

### Bugs recorrentes corrigidos nesta sessão
- `os_time` → `os_time_months` em todos os scripts que faltavam (07_SCANB, 11, 15)
- `geom_errorbarh()` deprecated → `geom_errorbar(orientation="y")` (08, 11)
- `concordance(Surv ~ matrix)$concordance` retorna vetor → usar `coxph()` + `concordance(model)$concordance`
- Colunas all-NA (er_status TCGA/GSE20685) → filtrar com `mean(!is.na()) >= 0.8`
- lubridate corrompido → reinstalar com `install.packages("lubridate", type="binary")`

### What's pending
- **17_render_manuscript_quarto.R** — render do QMD (requer Quarto instalado)
- **18_make_submission_bundle.R** — bundle final
- **figures/interpretations/** — agente em background completando interpretações
- Potencial: revisar I2=90.6% da meta-análise no texto do manuscrito (heterogeneidade esperada multi-plataforma)

