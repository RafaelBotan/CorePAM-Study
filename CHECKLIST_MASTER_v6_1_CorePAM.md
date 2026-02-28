# CHECKLIST_MASTER_v6_1 — Core-PAM (Memorial v6.1)

**Objetivo:** checklist operacional para execução reprodutível, auditável e “A1-proof”.  
**Regra:** nenhum passo é concluído sem (i) artefato, (ii) hash, (iii) commit quando aplicável.  
**Escopo:** Bloco Prognóstico (OS/DSS).

---

## 0) Regras globais (Freeze v6.1)

- [ ] Sem pooling de expressão entre coortes.
- [ ] Z-score: intra-coorte (ou intra-subcoorte em multi-plataforma).
- [ ] RNA-seq: counts → edgeR TMM → logCPM → Z-score.
- [ ] Microarray: log2-intensity → Z-score.
- [ ] Tempo: meses; dias → meses por 30.4375.
- [ ] Follow-up: Reverse Kaplan–Meier.
- [ ] Time ≤ 0: drop (reportar).
- [ ] Score (único):
  - [ ] score = (Σ w·z) / (Σ |w|) em genes presentes.
  - [ ] score_z = scale(score) para HR por 1 SD.
  - [ ] Direção: se HR(score)<1 → inverter sinal e registrar `score_direction`.
- [ ] Cobertura mínima do painel: ≥80% do painel final (limiar congelado em `analysis_freeze.csv`).
- [ ] C-index: reportar Craw e Cadj=max(Craw,1−Craw) (Cadj é diagnóstico/seleção).
- [ ] STRICT I/O: warning de leitura = erro (não prosseguir).

---

## 1) Repositório e estrutura

- [ ] Criar estrutura:
  - [ ] `01_docs/registry/`
  - [ ] `registry/`
  - [ ] `scripts/`
  - [ ] `data/raw/`, `data/processed/`
  - [ ] `results/corepam/`, `results/corepam_os/`
  - [ ] `figures/main/`, `figures/supp/`

---

## 2) Governança (antes de análise)

Salvar em `01_docs/registry/`:
- [ ] `cohort_manifest.csv`
- [ ] `analysis_freeze.csv` (inclui ΔC, alpha, K, seed folds, min_genes_fraction)
- [ ] `expression_input_spec.csv`
- [ ] `study_registry_template.csv` → copiar para `registry/study_registry.csv`

Commit: “Freeze v6.1 governance manifests (Core-PAM)”.

---

## 3) Preflight obrigatório (antes de pipeline grande)

- [ ] `scripts/07A_preflight_files_strict.R`
  - [ ] confirmar que todos os inputs leem sem warning/erro.
- [ ] Se METABRIC clínico falhar:
  - [ ] regenerar `data/processed/METABRIC_clinical_final.parquet` a partir de fonte íntegra.
  - [ ] atualizar `cohort_manifest.csv` com novo caminho e SHA.
- [ ] `scripts/07D_validate_one_cohort_corepam.R` (SCANB e METABRIC)
  - [ ] join > 0
  - [ ] cobertura ≥80%
  - [ ] score_direction registrado
  - [ ] HR/1SD plausível

---

## 4) Ingestão RAW e hashes

- [ ] Executar `01_download_raw_data.R` (coortes OS: SCAN-B via GSE96058 [ALL 3069 amostras = treino], METABRIC, TCGA-BRCA, GSE20685)
- [ ] Registrar SHA-256 dos RAW no `registry/study_registry.csv`.
- [ ] **Bloco pCR (após OS selado):** executar `01_download_raw_data_pcr.R` (GSE25066, GSE20194, GSE32646, I-SPY1, I-SPY2)

---

## 5) Harmonização clínica

- [ ] `02_harmonize_clinical_<COHORT>.R` por coorte (4 coortes: SCANB, TCGA_BRCA, METABRIC, GSE20685).
- [ ] Exportar `clinical_FINAL.parquet` por coorte.
- [ ] METABRIC: DSS (principal) + OS (sensibilidade).

---

## 6) Expressão pré-Z

- [ ] `03_expression_preprocess_<COHORT>.R` (4 coortes: SCANB, TCGA_BRCA, METABRIC, GSE20685)
- [ ] Exportar `expression_genelevel_preZ.parquet`.
- [ ] Exportar `gene_mapping_audit_<COHORT>.csv`.

---

## 7) Auditoria de genes

- [ ] `04_gene_audit_freeze.R`
- [ ] `results/supp/gene_audit_by_cohort.csv`
- [ ] Confirmar: PAM50 completo no SCANB (para derivação).

---

## 8) Derivação do Core-PAM (tamanho a determinar)

- [ ] `05_reduce_pam50_to_corepam_FINAL.R`
- [ ] Outputs em `results/corepam/`:
  - [ ] `pareto_df_cindex_oof.csv`
  - [ ] `CorePAM_weights.csv`
  - [ ] `CorePAM_model.rds`
  - [ ] `CorePAM_training_card.json`
  - [ ] `artifact_hashes.csv`
- [ ] Registrar hashes no `registry/study_registry.csv`.

---

## 9) Score por coorte e sobrevida

- [ ] `06_zscore_and_score_<COHORT>.R` (4 coortes: SCANB, TCGA_BRCA, METABRIC, GSE20685; usa `results/corepam/CorePAM_weights.csv`)
- [ ] `07_survival_analysis_<COHORT>.R` (4 coortes: SCANB, TCGA_BRCA, METABRIC, GSE20685)
- [ ] Exportar `survival_results_<COHORT>.csv` + figuras (KM em painéis).

---

## 10) Meta-análise e TRIPOD

- [ ] `08_meta_survival.R` (forest HR/1SD + I²/τ²)
- [ ] `11_incremental_value_and_dca.R` (ΔC-index + calibração 60m; DCA apenas sensibilidade)

---

## 11) QC final

- [ ] `14_qc_metabric_pca_forensics.R`
- [ ] `16_qc_text_vs_results_assert.R` (hard-fail)
- [ ] `17_render_manuscript_quarto.R` (hash do PDF no registry)
