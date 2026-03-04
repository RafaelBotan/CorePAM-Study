# Memorial de Rastreabilidade e Governança Científica  
**Projeto:** Core-PAM — Redução do PAM50 para um Núcleo Mínimo Funcional (tamanho a determinar) com Validação Prognóstica Multicoorte em Câncer de Mama  
**Versão:** 6.1 (TRIPOD-Ready / A1-Proof / Reexecutável / Auditável / Freeze “Core-PAM”)  
**Natureza:** estudo retrospectivo 100% público, multi-coorte, multi-plataforma (RNA-seq + microarray), com pipeline fechado e governança por hashes.

---

## 0) Sumário executivo (o que está congelado)

### 0.1 Entregável central do estudo (assinatura)
- **Assinatura final:** **Core-PAM** (subconjunto mínimo do PAM50; **número de genes não pré-fixado**).  
- **Critério de minimalismo:** seleção do **menor painel não-inferior** ao melhor modelo dentro do PAM50 na coorte de treino (SCAN-B), usando **C-index Out-of-Fold (OOF)** como métrica primária.  
- **Margem de não-inferioridade:** **ΔC = 0.010** (pré-especificada; pode ser ajustada apenas via emenda formal do freeze).  
- **Algoritmo de derivação:** Cox Elastic-Net (**α = 0.5**) com λ selecionado por envelope de não-inferioridade no path do glmnet e folds estratificados fixos.  
- **Comparabilidade cross-plataforma:** Z-score por gene **intra-coorte** e reporte do efeito como **HR por 1 SD do score**.  
- **Anti-gene-ausente:** score **re-escalonado** por soma de pesos absolutos (**Effective Gene Count**) em todas as coortes e plataformas.

> **Nota de nomenclatura:** o estudo, scripts, outputs e figuras **não** devem conter “PAM29/PAM32” ou qualquer número embutido no nome. O painel é referido como **Core-PAM** e o tamanho final é reportado apenas nos resultados/tabelas, nunca no identificador do método.

### 0.2 Dois blocos (independentes)
- **Bloco Prognóstico (OS/DSS):** objetivo primário (paper principal).  
- **Bloco Preditivo (NACT/pCR):** objetivo secundário (paper 2 opcional; execução sequencial sem contaminar o bloco OS).

---

## 1) Escopo, princípios e premissas (imutáveis)

### 1.1 Separação rígida de desfechos (proibida mistura)
- **Bloco Prognóstico:** endpoints de longo prazo (OS e/ou DSS/BCSS conforme disponibilidade).  
- **Bloco pCR:** endpoint binário de resposta patológica completa.  
- **Regra:** sem pooling de endpoints (OS ≠ DRFS ≠ DSS; pCR analisado separadamente).

### 1.2 Princípio de Não-Agregação (“Zero-Pooling”)
- Proibido combinar matrizes de expressão de coortes diferentes (pooling).  
- Proibido batch correction inter-coorte.  
- Permitido/obrigatório: normalização e padronização apenas **intra-coorte** (ou intra-subcoorte quando houver múltiplas plataformas no mesmo estudo).

### 1.3 Congelamento (“Freeze”) e anti-p-hacking
- Assinatura (genes + pesos + regra de score) é definida **uma única vez** na coorte de treino.  
- Coortes externas são exclusivamente validação (sem refit de coeficientes; sem tuning por coorte).  
- Sensibilidades (λ alternativo, cutoff alternativo, modelos clínicos alternativos) são pré-especificadas e rotuladas como sensibilidade.

### 1.4 Rastreabilidade forense (auditoria reexecutável)
- Todo artefato (raw → processado → resultados → figuras/tabelas) recebe hash SHA-256.  
- `registry/study_registry.csv` registra: caminho, SHA-256, script gerador, commit Git, timestamp, coorte, endpoint, parâmetros.  
- O manuscrito não contém números hard-coded: números são importados de CSVs finais.

---

## 2) Objetivos e claims (hierarquia formal)

### 2.1 Objetivo primário (Prognóstico)
Validar o **Core-PAM Risk Score** em coortes independentes, demonstrando:
- **Efeito:** HR por 1 SD do score (univariado e ajustado em CORE-A).  
- **Discriminação:** C-index por coorte (IC95% via bootstrap) + síntese por meta-análise.  
- **Calibração:** avaliação em horizonte fixo (5 anos) quando aplicável, com recalibração explícita.  
- **Valor incremental:** melhora versus baseline clínico mínimo (CORE-A).

### 2.2 Objetivos secundários
- Transportabilidade: consistência RNA-seq → microarray sem batch correction inter-coorte.  
- pCR (NACT): AUC e OR por coorte, síntese por meta-análise; I-SPY2 com regras especiais.  
- Subtipagem: frequências e concordância (PAM50 completo quando reconstruível vs painel reduzido).

---

## 3) Arquitetura de coortes e “Freeze de papéis”

### 3.1 Bloco Prognóstico (OS/DSS) — papéis congelados
- **Treino / derivação do painel (Discovery):** SCAN-B (RNA-seq).  
- **Validação externa RNA-seq:** TCGA-BRCA.
- **Validação microarray:** METABRIC, GSE20685 (Taiwan).

### 3.2 Regra “Leakage-Proof”
- SCAN-B = acesso GEO completo GSE96058 (n=3069); o pData do GEO NÃO possui coluna de separação training/validation (confirmado em 2026-02-28). Decisão (Option A, aprovada por Rafael Botan em 2026-02-28): usar TODAS as 3069 amostras como treinamento.
- GSE96058 deixa de existir como coorte de validação separada. Os scripts GSE96058 foram removidos.
- Sobreposição entre treino e validação: proibida. A regra leakage-proof aplica-se agora entre SCAN-B e as coortes de validação externas (TCGA-BRCA, METABRIC, GSE20685).
- Deduplicação determinística: cruzamento por Patient_ID + Sample_Barcode (ou chave equivalente) entre SCAN-B e cada coorte de validação.

### 3.3 Bloco NACT/pCR (secundário — implementado em scripts 19–23)
- **Coortes:** GSE25066 (N≈508, HGU133Plus2), GSE20194 (N≈278, HGU133Plus2),
  GSE32646 (N≈154, HGU133Plus2), I-SPY1/GSE22226 (N≈149, HGU133A).
- **Scripts:**
  - `19_download_pCR_raw_data.R` — download GEO series matrices (idempotente)
  - `20_prepare_pCR_<COHORT>.R` — harmonize clínica + preprocess expressão + Z-score + score por coorte
  - `21_pCR_logistic_analysis.R` — glm(pcr ~ score_z) + AUC DeLong + bootstrap OR
  - `22_meta_pCR.R` — meta-análise RE (DerSimonian-Laird) + FE (IVW) de log(OR)
  - `23_pCR_figures.R` — Fig6 (forest OR), FigS9 (ROC), FigS10 (distribuição por status)
- **Regras:**
  - Pesos CorePAM CONGELADOS — sem re-treino nas coortes pCR
  - Z-score intra-coorte (mesmo protocolo do bloco OS)
  - Sem pooling de expressão entre coortes
  - Direção: OR reportado conforme calculado (sem inversão para pCR)
  - I-SPY1 usa HGU133A (coverage ≥80% obrigatório)
- pCR não interfere no bloco prognóstico: scripts, outputs e manuscrito separados.

---

## 4) Inputs de expressão e normalização (fechado)

### 4.1 Unidade oficial de tempo
- Unidade oficial: meses.  
- Conversão dias → meses: 30.4375.  
- Follow-up mediana: Reverse Kaplan–Meier.  
- Tempos ≤ 0: drop com contagem reportada por coorte (QC).

### 4.2 RNA-seq (SCAN-B, TCGA-BRCA)
- Input congelado: raw counts.  
- Normalização pré-Z: edgeR TMM → logCPM (intra-coorte).  
- Padronização final: Z-score por gene intra-coorte.

### 4.3 Microarray (METABRIC, GSE20685, NACT GSEs)
- Input: intensidades log2 conforme disponibilizado.  
- Padronização final: Z-score por gene intra-coorte.

### 4.4 Mapeamento gênico (regra de ouro)
- IDs (Symbol/Entrez/RefSeq/Ensembl) → HGNC official.  
- Múltiplas probes → selecionar probe de maior variância intra-coorte.  
- Auditoria por coorte: `gene_mapping_audit_<coorte>.csv`.

---

## 5) Derivação do painel mínimo (treino) — protocolo “Core-PAM”

### 5.1 Filosofia de derivação (aprendizado incorporado)
- Seleção por estabilidade de gene pode falhar por divisão de votos em módulos colineares (colinearidade PAM50).  
- A derivação robusta do painel mínimo é feita por performance OOF e seleção do menor modelo não-inferior ao melhor dentro do PAM50.

### 5.2 Métrica primária de redução (congelada)
- C-index OOF de Harrell no SCAN-B com folds estratificados fixos.  
- C-index “orientation-free” para seleção: `Cadj = max(Craw, 1 − Craw)`.

### 5.3 Seleção do “XX mínimo” (sem pré-fixar tamanho)
- Calcular `Cadj(df)` para cada tamanho de modelo possível no λ-path.  
- Definir `Cmax = max(Cadj)` e selecionar o menor df com `Cadj >= Cmax − ΔC`.  
- Desempate: maior λ (mais parcimonioso) dentro do mesmo df.

### 5.4 Artefatos congelados da derivação
- `results/corepam/pareto_df_cindex_oof.csv`
- `results/corepam/CorePAM_weights.csv`
- `results/corepam/CorePAM_model.rds` (usar `cv_full$glmnet.fit` no λ escolhido)
- `results/corepam/CorePAM_training_card.json` (inclui df final, λ, ΔC, folds, hashes)
- `results/corepam/artifact_hashes.csv`

### 5.5 Regra de direção do score (congelada)
- Se `coxph(Surv ~ score)` em uma coorte produzir HR < 1, inverter: `score := -score`.  
- Registrar `score_direction` nos outputs.

---

## 6) Cálculo do score (aplicação cega) — regra única para todas as coortes

### 6.1 Z-score intra-coorte + re-escalonamento (Effective Gene Count)
Para cada coorte:
1) `G_present = genes painel ∩ genes coorte`  
2) `z_i` = Z-score por gene intra-coorte  
3) score:

\[
score=rac{\sum_{i \in G_{present}} w_i \cdot z_i}{\sum_{i \in G_{present}} |w_i|}
\]

4) `score_z = scale(score)` para reportar HR por 1 SD.

### 6.2 Critério mínimo e QC
- Cobertura mínima do painel por coorte: **≥80%** (ou limiar absoluto congelado no `analysis_freeze.csv`).  
- Reportar por coorte: genes_present, genes_missing, denom_sum_absw, score_direction, exclusões por `time<=0`.

---

## 7) Endpoints (regras fechadas)

### 7.1 OS (Bloco Prognóstico)
- event=1 óbito por qualquer causa; 0 censura.
- tempo em meses.

### 7.2 DSS/BCSS — **METABRIC endpoint PRIMÁRIO** (OS = sensibilidade)
- **METABRIC:** `dss_event=1` somente para óbitos por câncer de mama; demais = censura.
- `dss_time = OS_time` (mesmo tempo, só o event code difere).
- OS como análise de sensibilidade (FigS5).
- Sensibilidade adicional: competing risks (Fine–Gray) quando aplicável.

### 7.3 pCR — Bloco NACT (secundário)
- Endpoint binário: pCR = 1 (resposta patológica completa); RD = 0 (doença residual).
- Obtido de GEO pData por detecção automática de coluna (padrões: "pCR", "pathologic
  complete response", "RD", "residual disease").
- Análise: glm(pcr ~ score_z, family=binomial); OR + IC 95% bootstrap B=1000; AUC DeLong.

---

## 8) Plano estatístico — Bloco Prognóstico (primário)

### 8.1 Por coorte (sem pooling)
- Cox univariado: `Surv ~ score_z`.  
- Cox multivariado principal: CORE-A + score_z (quando disponível).  
- Discriminação: C-index (IC95% via bootstrap).  
- KM: mediana intra-coorte (principal) + quartis (sensibilidade).

### 8.2 Meta-análise
- Random-effects meta do log(HR) e SE entre coortes; heterogeneidade (I², τ²).  
- Leave-one-out (sensibilidade).

---

## 9) QC forense e contingências (obrigatórios)

### 9.1 Contingência — artefato corrompido (ex.: METABRIC RDS)
- Warning de leitura = erro (strict).
- Regenerar clínica METABRIC como parquet final (`data/processed/METABRIC_clinical_final.parquet`) a partir de fonte íntegra.
- Registrar novo hash no registry e atualizar `cohort_manifest.csv`.

### 9.2 C-index ajustado (C_adj = max(C, 1−C))
- Na seleção: usar `Cadj`.
- Na validação: aplicar `score_direction` (HR>1) e reportar C_adj consistentemente.

### 9.3 Congelamento de coeficientes
- Coeficientes extraídos do `glmnet.fit` no λ escolhido (full-data path).
- Evitar refit single-λ que altere df.

### 9.4 TCGA follow-up curto
- Sensibilidade por horizonte curto (24 meses) separada do resultado principal.

### 9.5 Micro-scripts (pré-execução)
- Preflight strict (I/O)  
- Join QC (IDs)  
- Teste mínimo por coorte (score + HR)

---

## 10) Deliverables finais (prova auditável)

### Bloco Prognóstico (OS/DSS)
- `registry/study_registry.csv`
- `01_docs/registry/analysis_freeze.csv`
- `results/corepam/` (pareto + pesos + modelo + training card + hashes)
- `results/main/` (Table3, meta, incremental value)
- `results/supp/` (survival per cohort, calibration, QC)
- `figures/main/` (Fig1–Fig5, bilingual EN+PT)
- `figures/supp/` (FigS1–FigS8, bilingual)

### Bloco pCR (NACT — secundário)
- `results/pcr/audit_gene_coverage_<COHORT>.csv` (cobertura CorePAM por coorte)
- `results/pcr/audit_pcr_definition.csv` (mapeamento pCR por coorte — ver §11.4)
- `results/pcr/pcr_validation_table.csv` (OR, AUC, n, events, cobertura — tabela consolidada)
- `results/pcr/meta_pCR_results.csv` / `pcr_meta_OR.csv` (pooled OR + heterogeneidade)
- `results/pcr/artifact_hash_manifest.csv` (SHA-256 de todos os outputs pCR)
- `figures/main/Fig6_Forest_pCR_OR_[EN|PT].pdf`
- `figures/supp/FigS9_ROC_pCR_[EN|PT].pdf`
- `figures/supp/FigS10_ScoreDist_pCR_[EN|PT].pdf`
- `01_Base_Pura_CorePAM/PROCESSED/pCR/<COHORT>/analysis_ready.parquet` (local, gitignored)

---

## 11) Bloco pCR — Protocolo Detalhado (Módulo Reprodutível Independente)

> **Nota:** as análises pCR são cientificamente independentes da validação prognóstica (OS) e
> introduzem pontos sensíveis para revisores (heterogeneidade de endpoint, AUC/OR CI,
> ajuste por covariáveis). Para não atrasar o paper OS, o bloco pCR é executado como módulo
> separado com manifesto, QA, scripts, tabelas e figuras próprios.

### 11.1 Coortes pCR (conjunto de decisão atual — congelado)

#### 11.1.1 Coortes primárias (corpo do artigo — congeladas)

| Coorte | GEO | Plataforma | N esperado | NACT |
|--------|-----|-----------|------------|------|
| GSE25066 | GSE25066 | HGU133Plus2 (GPL570) | ≈508 | EC ± docetaxel ± capecitabina |
| GSE32646 | GSE32646 | HGU133Plus2 (GPL570) | ≈154 | paclitaxel + FAC |
| GSE20194 | GSE20194 | HGU133Plus2 (GPL570) | ≈278 | paclitaxel/FAC |
| I-SPY1  | GSE22226 | HGU133A (GPL96)      | ≈149 | AC → paclitaxel |

**Regra editorial (imutável):** O módulo pCR primário é composto **exclusivamente** por essas
4 coortes. Quaisquer coortes adicionais entram **apenas** como análise exploratória ou de
sensibilidade — **nunca migram para o módulo primário**. Essa regra existe para:
(a) evitar comparações múltiplas não planejadas;
(b) manter robustez contra acusação de "picking" por revisores;
(c) preservar a capacidade de concluir o trabalho sem aguardar harmonização adicional.

#### 11.1.2 Coortes exploratórias / sensibilidade pCR (fora do corpo do artigo)

| Coorte | GEO | Plataforma | Status | Razão |
|--------|-----|-----------|--------|-------|
| I-SPY2 | TBD (GSE194040 ou similar) | Agilent 44K/32K | pending_harmonization | plataforma distinta; harmonização não concluída |

**Critérios mínimos para aceitar uma coorte como exploratória:**
1. endpoint pCR claramente mapeável no dicionário de auditoria;
2. cobertura gênica ≥ 80% (FREEZE$min_genes_fraction);
3. N suficiente para estimativas estáveis (≥ 50 amostras, ≥ 10 eventos pCR);
4. variáveis mínimas disponíveis (ER/HER2 para CORE-NACT incremental value).

**I-SPY2:** Candidata exploratória. Plataforma Agilent com ~990 amostras pré-tratamento.
Deve ser incorporada **somente** quando:
- GEO accession e mapeamento de colunas clínicas confirmados;
- cobertura gênica ≥ 80% verificada;
- script `20_prepare_pCR_ISPY2.R` testado em sessão limpa.
Resultado aparecerá apenas em seção de Análise de Sensibilidade (não no Abstract).

### 11.2 Preditor canônico (congelado)

- **Preditor:** CorePAM score (soma ponderada de genes; pesos congelados do treino OS).
- **Direção:** fixada no treino OS (maior score = maior risco prognóstico).
- **Regra crítica:** nunca inverter o sinal do score por coorte pCR para "forçar OR > 1".
  Se OR < 1, reportar e interpretar.
- **Fórmula:** Z-score intra-coorte por gene → score = Σ(wᵢ·zᵢ)/Σ|wᵢ| → score_z = scale(score).

### 11.3 QA Gates (regras rígidas por coorte)

**11.3.1 Integridade de arquivo**
- Arquivo de expressão legível (parquet/series_matrix.txt.gz), coluna ID identificada.
- Arquivo clínico contém: pCR binário mapeado + identificador alinhado à expressão.

**11.3.2 Cobertura de genes**
- Tabela obrigatória: Gene × present_in_matrix → `results/pcr/audit_gene_coverage_<COHORT>.csv`.
- Cobertura mínima: **≥ 80%** dos genes CorePAM presentes.
- Se < 80%: coorte excluída da meta-análise primária; aparece apenas como nota documentada
  de "cobertura insuficiente" no suplemento.

**11.3.3 Completude do desfecho**
- pCR mapeável a 0/1 com definição documentada (§11.4).
- Valores NA quantificados; se extremo (>30%), coorte excluída.

### 11.4 Mapeamento do Endpoint pCR — Auditoria (crítico para revisores)

**Deliverable obrigatório:** `results/pcr/audit_pcr_definition.csv`

Campos por coorte:
- `cohort`, `geo_accession`, `raw_variable_name`, `raw_label_levels`,
  `mapped_binary_rule`, `textual_definition` (da publicação/metadata GEO),
  `insitu_included` (sim/não/desconhecido).

**Heterogeneidade conhecida (deve ser documentada explicitamente):**
- GSE20194: pCR definido como ypT0 ypN0 (sem in-situ explícito).
- GSE25066, GSE32646, I-SPY1: frequentemente ypT0/is ypN0 (DCIS residual permitida).

**Análise de sensibilidade (pré-especificada):**
- Primário: usar definição canônica publicada por coorte.
- Sensibilidade: definição harmonizada "sem doença invasiva" quando reconstruível.

### 11.5 Endpoints estatísticos primários

**Por coorte:**
- Logística univariada: `glm(pcr ~ score_z, family = binomial)`
- OR por 1 DP (score_z) + IC 95% Wald + IC 95% bootstrap (B=1000, seed=42)
- AUC (ROC) + IC 95% DeLong

**Calibração** (diagnóstica, não filtro de seleção):
- Calibração logística: intercepto/slope; curva de calibração.
- Brier score (opcional).

### 11.6 Endpoints secundários — Valor incremental (CORE-NACT)

**Baseline mínimo (CORE-NACT) — apenas variáveis consistentemente disponíveis:**
- ER status (binário) + HER2 status (binário) + [grade/estágio clínico se presente em ≥2 coortes]

**Análise por coorte:**
- Modelo A: `pCR ~ ER + HER2 (+ ...)`
- Modelo B: `pCR ~ ER + HER2 (+ ...) + score_z`
- ΔAUC (B − A) com IC bootstrap
- Likelihood ratio test (modelos aninhados)

### 11.7 Pooling (pré-especificado)

- **Primário:** meta-análise de log(OR) por efeitos aleatórios (DerSimonian-Laird).
- **Sensibilidade:** IVW (efeitos fixos).
- AUC pooling: opcional; se feito, método declarado explicitamente e tratado como
  evidência de suporte.

### 11.8 Checklist anti-armadilha de revisores

1. AUC + IC **e** OR + IC reportados para **cada** coorte pCR (não apenas algumas).
2. Distinguir claramente OR univariado vs. ajustado (ambos podem ser apresentados, mas rotulados).
3. Documentar explicitamente diferenças na definição pCR (ypT0ypN0 vs. ypT0/isypN0).
4. Sem "direction hacking": sem inversão de sinal; sem escolha de cutoff post-hoc.
5. Coortes pCR são **validação only** — sem re-ajuste dos pesos usando desfecho pCR.

### 11.9 Outputs obrigatórios (results/pcr/)

| Arquivo | Gerado por | Conteúdo |
|---------|-----------|----------|
| `audit_gene_coverage_<COHORT>.csv` | 20_prepare_pCR_*.R | Gene × presente/ausente |
| `audit_pcr_definition.csv` | 21_pCR_logistic_analysis.R | Definição pCR por coorte |
| `pcr_validation_table.csv` | 21_pCR_logistic_analysis.R | OR, AUC, n, events, cobertura |
| `pCR_results_by_cohort.csv` | 21_pCR_logistic_analysis.R | Tabela completa por coorte |
| `pcr_meta_OR.csv` | 22_meta_pCR.R | OR poolado + heterogeneidade |
| `meta_pCR_results.csv` | 22_meta_pCR.R | RE + FE consolidado |
| `artifact_hash_manifest.csv` | 23_pCR_figures.R | SHA-256 de todos outputs pCR |

### 11.10 Regras Invioláveis (Hard Rules — sem graus de liberdade do pesquisador)

1. **Nunca re-treinar** pesos CorePAM em desfechos pCR.
2. **Nunca inverter** a direção do score por coorte pCR. Direção congelada no treino OS.
3. **Nunca escolher** cutpoints otimizando p-valores. Se dicotomização mostrada: mediana
   (pré-especificada) ou quartis (sensibilidade).
4. **Para cada coorte:** reportar OR+IC **e** AUC+IC (ou documentar explicitamente por que impossível).

### 11.11 Ordem de Execução (referência de implementação)

- **Step 0** — Manifesto + mapeamento
  `config/pcr_cohort_manifest.csv` (cohort, expr_path, clin_path, endpoint_var, mapping_rule)
  `results/pcr/audit_pcr_definition.csv` (criado ou atualizado)

- **Step 1** — QA por coorte (`20_prepare_pCR_<COHORT>.R`)
  Leitura expressão + clínica → padronização ID → mapeamento pCR → Z-score → score
  Salva: `audit_gene_coverage_<COHORT>.csv`, linha em `audit_pcr_definition.csv`
  Gate: coverage ≥80% + pCR mapeável a 0/1

- **Step 2** — Métricas primárias por coorte (`21_pCR_logistic_analysis.R`)
  Z-score intra-coorte → score_z → logística univariada → OR + IC + AUC DeLong
  Append à `pcr_validation_table.csv`

- **Step 3** — Valor incremental (`21_pCR_logistic_analysis.R`)
  Modelo A: `pCR ~ ER + HER2 (+...)` → Modelo B: `+ score_z`
  ΔAUC bootstrap + LRT → `pcr_incremental_table.csv`

- **Step 4** — Meta-análise (`22_meta_pCR.R`)
  Pool log(OR) RE+FE → `pcr_meta_OR.csv`

- **Step 5** — Figuras (`23_pCR_figures.R`)
  Salvar em `figures/pcr/` (resolução journal-compliant ≥300 dpi) + `artifact_hash_manifest.csv`

### 11.12 Campos mínimos de reporte por coorte (linha de tabela)

| Campo | Descrição |
|-------|-----------|
| cohort | identificador da coorte |
| endpoint | "pCR" |
| n_total | N total analisado |
| n_pcr1 | N pCR=1 |
| n_rd | N pCR=0 |
| genes_total | N genes CorePAM panel |
| genes_present | N genes presentes na coorte |
| genes_missing | genes ausentes |
| coverage_pct | % cobertura |
| or_per_1SD | OR por 1 DP score_z |
| or_ci_lo95 | IC 95% inferior |
| or_ci_hi95 | IC 95% superior |
| p_value | p-valor logístico |
| auc | AUC ROC |
| auc_lo95 | IC 95% AUC inferior |
| auc_hi95 | IC 95% AUC superior |
| baseline_covariates | variáveis usadas no modelo ajustado |
| pcr_definition_note | nota sobre definição pCR |

### 11.13 Condições de parada (stop conditions)

- **pCR não mapeável:** parar coorte, não adivinhar. Aguardar resolução em `audit_pcr_definition.csv`.
- **Cobertura < 80%:** marcar como "fail coverage", pular inferência.
- **Logística não-convergente:** logar; usar fallback penalizado (ridge) **somente como sensibilidade**.

### 11.14 Lista de Figuras pCR (v1 — mínimo obrigatório)

**Figuras principais** (bilingual EN+PT; salvar em `figures/pcr/`):

| Fig | Arquivo | Conteúdo |
|-----|---------|----------|
| Fig_pCR1 | `Fig_pCR1_Forest_OR_[EN\|PT].pdf/png` | Forest plot OR por 1 DP + pooled RE+FE |
| Fig_pCR2 | `Fig_pCR2_ROC_AUC_Panel_[EN\|PT].pdf/png` | AUC (com IC) por coorte — small multiples |
| Fig_pCR3 | `Fig_pCR3_pCR_Rate_by_Quartile_[EN\|PT].pdf/png` | Taxa pCR por quartil de score_z por coorte |
| Fig_pCR4 | `Fig_pCR4_Calibration_[EN\|PT].pdf/png` | Observado vs previsto pCR — maior coorte |

**Figuras suplementares:**

| Fig | Arquivo | Conteúdo |
|-----|---------|----------|
| FigS_pCR1 | `FigS_pCR1_GeneCoverage.pdf/png` | Heatmap/tabela cobertura CorePAM por coorte |
| FigS_pCR2 | `FigS_pCR2_pCRDefinition.pdf/png` | Tabela definição pCR por coorte (ypT0ypN0 vs ypT0/isypN0) |
| FigS_pCR3 | `FigS_pCR3_DeltaAUC.pdf/png` | ΔAUC baseline vs baseline+score (bootstrap IC) |
| FigS_pCR4 | `FigS_pCR4_DCA.pdf/png` | Decision curve (opcional; apenas se pré-especificado) |
| FigS_pCR5 | `FigS_pCR5_Sensitivity_NoStdDef.pdf/png` | Sensibilidade: excluir coortes c/ definição não-padrão |

**Tabelas:**

| Tabela | Arquivo | Conteúdo |
|--------|---------|----------|
| T_pCR_1 | `T_pCR_1_Cohort_Characteristics.csv` | Características das coortes + definição pCR + N pCR/RD |
| T_pCR_2 | `T_pCR_2_OR_AUC.csv` | OR univariado + AUC por coorte |
| T_pCR_3 | `T_pCR_3_Incremental_Value.csv` | ΔAUC + IC bootstrap (baseline vs baseline+score) |

### 11.15 Config/Manifesto (pré-requisito da execução)

Arquivo: `config/pcr_cohort_manifest.csv`

Campos: cohort, geo_accession, platform, expr_path, clin_path, pcr_raw_var, pcr_mapping_rule,
pcr_definition_text, insitu_included, timepoint_notes, n_expected.

**Heterogeneidade pré-documentada:**
- GSE20194: pCR = ypT0 ypN0 (sem in-situ explícito)
- GSE25066, GSE32646, I-SPY1: pCR = ypT0/is ypN0 (DCIS residual permitida)

### 11.16 Modos de falha + Plano-B (contingências pré-autorizadas)

| Situação | Ação |
|----------|------|
| Arquivo ilegível/corrompido | Congelar SHA-256 + size + path; marcar coorte "unavailable" |
| Cobertura < 80% | Excluir da meta-análise primária; documentar no suplemento |
| Modelo logístico não-convergente | Fallback para logística penalizada (ridge, λ fixo) nessa coorte; primário mantido onde convergente |
| Endpoint ambíguo | Coorte pausada até mapeamento pCR resolvido e documentado em `audit_pcr_definition.csv` |
