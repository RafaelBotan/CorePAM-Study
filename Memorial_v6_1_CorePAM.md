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
- **Validação externa RNA-seq:** TCGA-BRCA, GSE96058.  
- **Validação microarray:** METABRIC, GSE20685 (Taiwan).

### 3.2 Regra “Leakage-Proof” (SCAN-B vs GSE96058)
- Sobreposição entre treino e validação: proibida.  
- Deduplicação determinística: cruzamento por Patient_ID + Sample_Barcode (ou chave equivalente).  
- Artefato obrigatório: `results/leakage_check_scanb_vs_gse96058.csv` com hashes das listas de IDs e interseção=0.

### 3.3 Bloco NACT/pCR (secundário)
- GSE25066, GSE20194, GSE32646, I-SPY1 e I-SPY2 (pendente).  
- pCR não interfere no bloco prognóstico: scripts, outputs e manuscrito podem ser separados.

---

## 4) Inputs de expressão e normalização (fechado)

### 4.1 Unidade oficial de tempo
- Unidade oficial: meses.  
- Conversão dias → meses: 30.4375.  
- Follow-up mediana: Reverse Kaplan–Meier.  
- Tempos ≤ 0: drop com contagem reportada por coorte (QC).

### 4.2 RNA-seq (SCAN-B, TCGA, GSE96058)
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
- `results/corepam/CorePAM_training_card.json` (inclui df final, λ, ΔC, seed, folds, hashes)

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

### 7.1 OS
- event=1 óbito por qualquer causa; 0 censura.  
- tempo em meses.

### 7.2 DSS/BCSS (referência METABRIC)
- `dss_time = OS_time`.  
- `dss_event=1` somente causa câncer; demais = censura no tempo do óbito.  
- Sensibilidade: competing risks (Fine–Gray) quando aplicável.

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

## 9) QC forense e Planos B (obrigatórios)

### 9.1 Plano B — artefato corrompido (ex.: METABRIC RDS)
- Warning de leitura = erro (strict).  
- Regenerar clínica METABRIC como parquet final (`data/processed/METABRIC_clinical_final.parquet`) a partir de fonte íntegra.  
- Registrar novo hash no registry e atualizar `cohort_manifest.csv`.

### 9.2 Plano B — inversão de direção / C-index < 0.5
- Na seleção: usar `Cadj`.  
- Na validação: aplicar `score_direction` (HR>1) e reportar.

### 9.3 Plano B — discrepância de df no refit
- Congelar coeficientes a partir do `glmnet.fit` no λ escolhido (full-data path).  
- Evitar refit single-λ que altere df.

### 9.4 Plano B — TCGA follow-up curto
- Sensibilidade por horizonte curto (24 meses) separada do resultado principal.

### 9.5 Micro-scripts (pré-execução)
- Preflight strict (I/O)  
- Join QC (IDs)  
- Teste mínimo por coorte (score + HR)

---

## 10) Deliverables finais (prova auditável)
- `registry/study_registry.csv`  
- `cohort_manifest.csv`  
- `analysis_freeze.csv`  
- `results/corepam/` (pareto + pesos + modelo + training card + hashes)  
- `results/corepam_os/` (resultados por coorte + meta + figuras principais/suplementares)
