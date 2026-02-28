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
cruzamento entre coortes. Para todas as coortes RNA-seq (SCAN-B, GSE96058,
TCGA-BRCA, GSE20685), o desfecho primário foi a sobrevida global (OS), definida
como o intervalo entre o diagnóstico e o óbito por qualquer causa ou a última
observação (censura). Para o METABRIC, o desfecho primário foi a sobrevida
específica por doença (DSS), codificada como evento=1 exclusivamente para óbito
por doença; óbitos por outras causas foram censurados no momento do óbito.
Os tempos de seguimento foram expressos em meses (divisor: 30,4375 dias/mês).
Observações com tempo ≤ 0 foram excluídas e contabilizadas. A separação entre
as amostras de treinamento (SCAN-B) e de validação (GSE96058) foi realizada
utilizando a coluna de split do pData do GEO, com verificação formal de ausência
de sobreposição (intersecção = 0).

---

## Script 03 — Pré-processamento de Expressão Gênica

O pré-processamento foi realizado de forma estritamente independente para cada
coorte, sem pooling ou correção de batch inter-coorte. Para as coortes RNA-seq
(SCAN-B, GSE96058, TCGA-BRCA), as contagens brutas foram normalizadas pelo método
TMM (Trimmed Mean of M-values) utilizando o pacote edgeR, seguido de transformação
logCPM com pseudo-contagem de 1. Para as coortes de microarray (METABRIC: Illumina;
GSE20685: Affymetrix HGU133A), os valores de intensidade foram utilizados na escala
log2 como obtidos, com aplicação de log2(x+1) quando o percentil 95 excedia 20,
indicando dados não transformados. O mapeamento dos identificadores Ensembl para
símbolos HGNC oficiais foi realizado via biomaRt com cache local em RDS. Os
identificadores Affymetrix foram mapeados via pacote de anotação da plataforma
correspondente. Em casos de múltiplas sondas mapeando para o mesmo símbolo HGNC,
foi retida a sonda com maior variância intra-coorte. Genes com fração de dados
ausentes superior a 20% foram removidos.

---

## Script 04 — Auditoria de Cobertura do Painel PAM50

Antes da derivação do Core-PAM, foi realizada uma auditoria transversal da
cobertura dos 50 genes do painel PAM50 em todas as cinco coortes. Para a coorte
de treinamento (SCAN-B), o critério go/no-go exigiu a presença de todos os 50 genes.
Para as coortes de validação, o critério mínimo foi de 80% do painel PAM50 presente,
sendo a cobertura do Core-PAM verificada novamente após a derivação do painel.
Os resultados foram exportados em uma tabela gene × coorte (formato largo) e em
uma tabela de resumo de cobertura por coorte, ambas registradas com seus hashes
SHA-256.

---

## Script 05 — Derivação do Painel Core-PAM

O painel Core-PAM foi derivado exclusivamente na coorte de treinamento SCAN-B,
sem qualquer acesso às coortes de validação. Utilizou-se regressão de Cox com
regularização elástica (α = 0,5, pacote glmnet), com validação cruzada K=10 dobras
estratificadas por evento (semente=42). A estratificação garantiu distribuição
proporcional de eventos e censuras em cada dobra. O C-index fora da dobra (OOF)
foi calculado para cada valor de λ na trajetória do modelo, utilizando predições
lineares mantidas pelo argumento `keep=TRUE` do `cv.glmnet`. Para cada λ, calculou-se
o C-index ajustado pela orientação: Cadj = max(C_bruto, 1 − C_bruto), tornando a
métrica invariante à inversão de escala do score. O painel final foi selecionado
como o modelo com menor número de genes não-zero cujo Cadj ficou dentro de ΔC=0,010
do Cadj máximo observado (critério de não-inferioridade pré-especificado).
Se HR(score) < 1, o sinal dos pesos foi invertido para que scores elevados
correspondam a pior prognóstico. O painel resultante, seus pesos e os parâmetros
de derivação foram registrados nos arquivos `CorePAM_weights.csv`,
`CorePAM_model.rds` e `CorePAM_training_card.json`.

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
log(HR) univariados e seus erros-padrão de todas as cinco coortes. A heterogeneidade
foi quantificada pelo I² e τ². A análise leave-one-out foi conduzida como análise de
sensibilidade.

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
