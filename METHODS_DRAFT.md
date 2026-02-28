# METHODS_DRAFT.md — Rascunho da Seção de Métodos (Portuguese)

**Convenção:** Após a conclusão de cada script, um parágrafo em português acadêmico
é adicionado aqui para uso direto na tese/manuscrito. Os parágrafos são gerados
automaticamente ou manualmente conforme o pipeline avança.

---

## Script 00 — Configuração e Parâmetros Congelados

Os parâmetros analíticos foram pré-especificados e congelados antes do início
de qualquer análise, seguindo o princípio de inferência confirmatória. O congelamento
incluiu: critério de não-inferioridade (ΔC = 0,010), parâmetro de regularização
elástica (α = 0,5), número de dobras de validação cruzada (K = 10), semente aleatória
(42) e fração mínima de genes do painel exigida nas coortes de validação (80%).
Os parâmetros foram armazenados no arquivo `01_docs/registry/analysis_freeze.csv`
e carregados via o script `00_setup.R` por todos os demais scripts do pipeline.
Verificações de integridade (SHA-256) foram computadas para todos os artefatos gerados.

---

## Script 01 — Coleta e Verificação de Integridade dos Dados Brutos

Os dados foram obtidos integralmente de fontes públicas. A coorte de treinamento
(SCAN-B, acesso GEO: GSE96058) foi baixada via o pacote GEOquery. Os dados de
expressão gênica do TCGA-BRCA foram recuperados via TCGAbiolinks como objeto
SummarizedExperiment contendo contagens brutas do alinhador STAR. Os dados do
METABRIC foram obtidos do cBioPortal por URL direta. O conjunto GSE20685 foi
recuperado do GEO. Após cada download, foi calculado o resumo criptográfico SHA-256,
registrado em `registry/study_registry.csv`, garantindo rastreabilidade e
reprodutibilidade dos dados de origem.

---

## Script 02 — Harmonização Clínica e Definição de Desfechos

Os dados clínicos de cada coorte foram harmonizados de forma independente, sem
cruzamento entre coortes. A coorte de treinamento SCAN-B corresponde ao acesso
GEO GSE96058 em sua totalidade (n=3.069 amostras). O pData do GEO foi inspecionado
e não contém coluna de separação training/validation; portanto, todas as 3.069
amostras foram alocadas ao conjunto de treinamento (Opção A, aprovada em 2026-02-28).
Para as coortes RNA-seq (SCAN-B, TCGA-BRCA) e microarray de validação (GSE20685),
o desfecho primário foi a sobrevida global (OS). Para o METABRIC, o desfecho
primário foi a sobrevida específica por doença (DSS); óbitos por outras causas
foram censurados no momento do óbito. Os tempos foram expressos em meses
(divisor: 30,4375 dias/mês; anos × 12 para GSE20685).

**Resultados (scripts 02–04, execução 2026-02-28):**

| Coorte | N | Eventos OS | % | Mediana FU all | Mediana FU censurados | Max FU | PAM50 |
|--------|---|------------|---|----------------|----------------------|--------|-------|
| SCAN-B | 3.069 | 322 | 10,5% | 53,0 mo | 54,9 mo | 81,3 mo | 50/50 |
| TCGA-BRCA | 1.072 | 150 | 14,0% | 28,3 mo | 25,6 mo | 282,7 mo | 49/50 |
| METABRIC | 1.980 | 1.144 OS / 646 DSS | 57,8% | 116,5 mo | **157,9 mo** | 355,2 mo | 50/50 |
| GSE20685 | 327 | 83 | 25,4% | 97,2 mo | 110,4 mo | 169,2 mo | 44/50† |

† ANLN, CXXC5, GPR160, NUF2, TMEM45B, UBE2T ausentes: sem probes no HGU133A (GPL96);
  44/50 (88%) é o teto absoluto da plataforma — acima do critério mínimo de 80%.

A arquitetura final compreende 4 coortes: SCAN-B (treino) + 3 validação.

---

## Script 03 — Pré-processamento de Expressão Gênica

O pré-processamento foi realizado de forma estritamente independente para cada
coorte, sem pooling ou correção de batch inter-coorte. Para as coortes RNA-seq
(SCAN-B, TCGA-BRCA), as contagens brutas foram normalizadas pelo método
TMM (Trimmed Mean of M-values) utilizando o pacote edgeR, seguido de transformação
logCPM com pseudo-contagem de 1. Para as coortes de microarray (METABRIC: Illumina;
GSE20685: Affymetrix HGU133A), os valores de intensidade foram utilizados na escala
log2 como obtidos, com aplicação de log2(x+1) quando o percentil 95 excedia 20,
indicando dados não transformados. O mapeamento dos identificadores Ensembl para
símbolos HGNC oficiais foi realizado via `org.Hs.eg.db` + `AnnotationDbi`
(anotação offline, sem dependência de rede), com cache local em RDS para
reutilização. O pacote `biomaRt` não foi utilizado por apresentar conflito de
dependências nesta instalação. Os identificadores Affymetrix foram mapeados via
pacote de anotação da plataforma correspondente (`hgu133a.db` para GSE20685). Em casos de múltiplas sondas mapeando para o mesmo símbolo HGNC,
foi retida a sonda com maior variância intra-coorte. Genes com fração de dados
ausentes superior a 20% foram removidos.

---

## Script 04 — Auditoria de Cobertura do Painel PAM50

Antes da derivação do Core-PAM, foi realizada uma auditoria transversal da
cobertura do painel PAM50 (50 genes) em todas as 4 coortes. Para a coorte
de treinamento (SCAN-B), o critério go/no-go exigiu a presença de todos os 50 genes
(resultado: 50/50; GO). Para as coortes de validação, o critério mínimo foi de 80%:
TCGA-BRCA 49/50 (98%; MIA ausente — gene não mapeado pelo org.Hs.eg.db nesta coorte),
METABRIC 50/50 (100%), GSE20685 44/50 (88%; ANLN, CXXC5, GPR160, NUF2, TMEM45B, UBE2T
ausentes — limitação de plataforma: esses 6 genes não possuem probes no array HGU133A,
confirmado via `hgu133a.db`; 44/50 é o teto absoluto desta plataforma).
Os resultados foram exportados em tabela gene × coorte (formato largo) e em
tabela de resumo de cobertura por coorte, ambas registradas com hashes SHA-256.

---

## Script 05 — Derivação do Painel Core-PAM

O painel Core-PAM foi derivado exclusivamente na coorte de treinamento SCAN-B
(n=3.069; Opção A), sem qualquer acesso às coortes de validação. Utilizou-se
regressão de Cox com regularização elástica (α = 0,5, pacote glmnet) com K=10 dobras.

A atribuição de amostras às dobras foi **determinística** — baseada nos primeiros
8 dígitos hexadecimais do SHA-256 do identificador de paciente, com estratificação
por status de evento. Essa abordagem garante que a mesma amostra sempre receba a
mesma dobra independentemente da ordem dos dados, eliminando a dependência de uma
semente aleatória como motor de seleção de genes. Os parâmetros de K e α foram
pré-especificados e congelados antes de qualquer análise.

O C-index de Harrell fora da dobra (OOF) foi calculado para cada nível único de
graus de liberdade (df = número de genes não-zero) na trajetória de regularização.
Para cada df, foi selecionado o λ mais parcimonioso (maior λ com aquele df).
O C-index ajustado pela orientação — Cadj = max(C_bruto, 1 − C_bruto) — torna
a métrica invariante à inversão de escala do score. O painel final (Core-PAM) foi
selecionado como o **menor df** com Cadj ≥ Cadj_máximo − ΔC (ΔC = 0,010),
aplicando o critério de não-inferioridade pré-especificado. O número de genes
não foi pré-especificado: é derivado dos dados e da regra de não-inferioridade.

Os pesos finais foram extraídos da trajetória completa do glmnet no λ escolhido
(não do cv.glmnet), evitando variação de df entre treino completo e OOF.
O máximo de genes não-zero avaliado na trajetória foi df=49 (de 50 candidatos):
com regularização elastic-net (α=0,5), a penalidade L1 manteve um gene com
coeficiente zero mesmo no λ mínimo — comportamento esperado do elastic-net que
não implica limitação metodológica.
O painel resultante, seus pesos e os metadados de derivação foram registrados nos
arquivos `CorePAM_weights.csv`, `CorePAM_model.rds`, `selected_CorePAM_summary.json`
e `pareto_df_cindex_oof.csv` (curva Pareto df × C-index OOF).

---

## Script 06 — Cálculo do Score Core-PAM nas Coortes

Para cada coorte, o score Core-PAM foi calculado de forma independente em três etapas:
(i) Z-score intra-coorte por gene (média e desvio-padrão calculados dentro da própria
coorte, sem referência externa); (ii) cálculo do score ponderado pela contagem efetiva
de genes: score = Σ(wᵢ × zᵢ) / Σ|wᵢ|, considerando apenas os genes do painel
presentes na coorte; (iii) padronização do score para HR por 1 DP (score_z = scale(score)).
A cobertura mínima de 80% dos genes do painel foi verificada em cada coorte. Genes
ausentes foram contabilizados e reportados. A direção do score foi verificada por modelo
de Cox univariado; quando HR < 1, o sinal foi invertido e o parâmetro `score_direction`
foi registrado como "inverted". Os resultados foram exportados como `analysis_ready.parquet`
por coorte.

---

## Scripts 07–08 — Análise de Sobrevida e Meta-análise

Em cada coorte de validação, o score Core-PAM foi avaliado como preditor contínuo
de sobrevida por modelo de Cox univariado e multivariado (ajustado para o modelo
CORE-A: idade e status de receptor de estrogênio). O C-index de Harrell foi estimado
com intervalo de confiança de 95% por bootstrap (n=1.000 reamostras, semente=42).
Curvas de Kaplan-Meier foram construídas dicotomizando as amostras no ponto de corte
da mediana intra-coorte (primário) e nos quartis (sensibilidade). Uma meta-análise
de efeitos aleatórios (método REML, pacote metafor) foi conduzida combinando os
log(HR) univariados e seus erros-padrão das três coortes de validação (TCGA-BRCA,
METABRIC, GSE20685). A heterogeneidade foi quantificada pelo I² e τ². A análise
leave-one-out foi conduzida como análise de sensibilidade.

---

## Script 11 — Valor Incremental e DCA

O valor incremental do Core-PAM além do modelo clínico CORE-A foi avaliado pelo
delta C-index (ΔC = C(CORE-A + score) − C(CORE-A)), com bootstrap para IC95%.
A análise de decisão clínica (DCA) foi conduzida como análise de sensibilidade
para avaliar o benefício líquido em limiares de probabilidade clinicamente relevantes.

---

## Scripts 13–16 — Controle de Qualidade

A qualidade dos dados e análises foi avaliada por quatro componentes complementares:
(i) correlações off-diagonal entre scores das coortes, verificando coerência sem
pooling; (ii) análise de componentes principais forense no METABRIC, para detecção
de efeitos de lote ou amostras outliers; (iii) verificação estrutural de esquema e
intervalos de valores em todos os arquivos `analysis_ready.parquet`; e (iv)
verificação de que nenhum número relevante foi hard-coded no manuscrito, assegurando
que todos os resultados numéricos são lidos diretamente dos arquivos CSV/JSON gerados
pelo pipeline.

---

<!-- Novos parágrafos são adicionados abaixo conforme os scripts são concluídos -->
