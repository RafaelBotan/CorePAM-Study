# Interpretação: FigS — Distribuição do Escore CorePAM por Coorte

## O que é este gráfico?

Este gráfico suplementar mostra a **distribuição (histograma ou densidade) do escore CorePAM** em cada uma das quatro coortes do estudo — SCAN-B (treinamento), TCGA-BRCA, METABRIC e GSE20685 (validação). Tipicamente apresentado como painel múltiplo (facets), com uma curva de densidade ou histograma por coorte.

O escore CorePAM é calculado como uma média ponderada dos Z-scores intra-coorte: `score = Σ(wᵢ·zᵢ) / Σ|wᵢ|`. Por ser baseado em Z-scores (que por definição têm média ≈ 0 e DP ≈ 1 por gene), o escore resultante deve ter distribuição aproximadamente centrada em zero em cada coorte.

## Como ler este gráfico

- **Cada painel:** Uma coorte (SCAN-B, TCGA-BRCA, METABRIC, GSE20685).
- **Eixo X:** Valor do escore CorePAM (tipicamente varia de ~-2 a +2, em escala de DP).
- **Eixo Y:** Densidade (ou contagem) — altura da curva em cada ponto do eixo X.
- **Forma da distribuição:** Curva de sino (normal/gaussiana) indica distribuição simétrica. Bimodalidade (dois picos) indicaria dois subgrupos distintos — biologicamente correspondente à divisão luminal/basal.
- **Mediana (linha vertical):** Ponto de corte para a análise KM principal — divide exatamente 50% de cada grupo.
- **Sobreposição de cores:** Se colorido por ER status ou subtipo PAM50, mostra onde cada grupo se posiciona na distribuição.

## Achados principais

Características esperadas da distribuição em cada coorte:

- **SCAN-B:** Distribuição larga com possível bimodalidade sutil, refletindo a diversidade de subtipos numa coorte de RNA-seq com boa representação de todos os subtipos.

- **TCGA-BRCA:** Distribuição similar a SCAN-B — ambas são RNA-seq de população ocidental. A comparação visual entre as duas confirma que os Z-scores são adequadamente calibrados.

- **METABRIC:** Distribuição potencialmente mais estreita ou com pico diferente — reflete a população europeia mais envelhecida com predominância de subtipos luminais.

- **GSE20685:** Distribuição potencialmente mais estreita por N menor (327). A forma geral deve ser similar às outras coortes se a normalização foi adequada.

**Importante:** Por design do Z-score intra-coorte, a média de cada coorte deve ser ≈ 0 e o DP ≈ variável (dependendo dos pesos do modelo). A comparabilidade das distribuições entre coortes **não** é garantida e **não é necessária** — o que importa é que dentro de cada coorte a distribuição seja razoável e que haja variação suficiente para discriminação prognóstica.

## O que isso significa?

Esta figura serve para verificar que o processo de cálculo do escore foi executado corretamente em cada coorte. Especificamente:

1. **Ausência de outliers extremos:** Escores muito distantes (>4 DP) poderiam indicar amostras problemáticas ou erros de processamento.

2. **Distribuição contínua:** Uma distribuição contínua (sem buracos ou gaps) confirma que o escore funciona como um preditor contínuo.

3. **Variação suficiente:** Uma distribuição colapsada (todos os valores muito próximos) indicaria que o escore não está discriminando bem.

4. **Consistência entre coortes:** Embora as distribuições não precisem ser idênticas (Z-scores são intra-coorte), padrões muito diferentes poderiam indicar problemas de processamento em uma coorte específica.

5. **Bimodalidade (se presente):** Uma distribuição bimodal é biologicamente esperada em câncer de mama — um pico luminal (escores baixos, alta expressão de ESR1/PGR/BCL2) e um pico basal (escores altos, alta expressão de EXO1/CENPF/KRT5). A bimodalidade confirmaria que o escore está capturando a divisão luminal/basal.

## Contexto no estudo

Esta figura suplementar é essencial para **transparência metodológica**. Antes de qualquer análise de sobrevida, é boa prática verificar a distribuição das variáveis preditoras. Distribuições problemáticas (muito assimétricas, multimodais de forma inesperada, ou com outliers extremos) poderiam invalidar suposições dos modelos Cox.

A figura também é relevante para interpretar as análises de quartis (figuras KM por quartis): a distribuição mostra se os quartis têm tamanhos aproximadamente iguais e onde estão os pontos de corte dos quartis em relação à distribuição.

## Pontos de atenção

- **Comparação visual entre coortes:** As distribuições não são diretamente comparáveis em escala absoluta (por causa do Z-score intra-coorte). Comparar formas (bimodalidade, assimetria) é válido; comparar posições absolutas não é.
- Uma distribuição **unimodal e simétrica** seria adequada para análises Cox que assumem efeito linear do preditor. Uma distribuição **bimodal** pode violar implicitamente essa suposição se os dois grupos tiverem efeitos muito diferentes.
- O escore CorePAM é a **soma ponderada normalizada** — portanto, a distribuição é mais comportada que a de escores brutos. Isso é uma vantagem da normalização por Σ|wᵢ|.
- A **escala do eixo X** pode variar entre coortes dependendo da variância real dos Z-scores dos genes CorePAM em cada população. Isso é esperado e não é um problema.
