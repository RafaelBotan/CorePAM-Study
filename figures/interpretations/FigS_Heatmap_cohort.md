# Interpretação: FigS — Heatmaps de Expressão dos Genes CorePAM por Coorte (TCGA-BRCA, METABRIC, GSE20685)

## O que é este gráfico?

Este arquivo cobre os **três heatmaps de expressão** dos genes CorePAM, um para cada coorte de validação (TCGA-BRCA, METABRIC, GSE20685). Cada heatmap mostra os **Z-scores de expressão dos 24 genes CorePAM** em todas as amostras da coorte correspondente, com amostras ordenadas pelo escore CorePAM calculado.

Heatmaps de expressão são os "radiografias moleculares" de um conjunto de amostras — eles permitem visualizar padrões de co-expressão entre genes e padrões de estratificação entre amostras de uma única vez.

## Como ler cada heatmap

### Estrutura geral:
- **Linhas:** Os 24 genes do CorePAM.
- **Colunas:** Pacientes da coorte (uma coluna por paciente).
- **Cor de cada célula:** Z-score de expressão do gene i na amostra j.
  - **Vermelho/quente:** Expressão alta (Z > 0, acima da média da coorte).
  - **Azul/frio:** Expressão baixa (Z < 0, abaixo da média da coorte).
  - **Branco/neutro:** Expressão média (Z ≈ 0).

### Ordenação das colunas (amostras):
- As amostras são ordenadas **do menor para o maior escore CorePAM** (esquerda = baixo risco, direita = alto risco).
- Essa ordenação cria o padrão visual característico: genes de peso negativo aparecem vermelhos à esquerda (alta expressão em baixo escore) e azuis à direita; genes de peso positivo aparecem azuis à esquerda e vermelhos à direita.

### Anotações no topo (row/column annotations):
- **Escore CorePAM:** Barra de gradiente colorida acima do heatmap, confirmando a ordenação.
- **Status de ER:** Barra de cores mostrando ER+ (tipicamente azul) e ER− (tipicamente vermelho) por paciente.
- **Subtipo PAM50** (se disponível): Barra de cores com Luminal A/B, HER2-enriched, Basal-like.
- **Evento de sobrevida** (opcional): Óbito vs. censura.

### Agrupamento de genes (linhas):
- Genes podem ser organizados por **hierarquia de clustering** (dendrograma ao lado) ou por **pesos positivos vs. negativos**.
- A separação visual entre o bloco de genes com peso positivo (marcadores basais) e negativo (marcadores luminais) é o padrão esperado.

## Achados principais por coorte

### Heatmap TCGA-BRCA
- **N = 1.072 amostras** — o heatmap mais largo dos três.
- Padrão esperado: transição gradual da esquerda (Q1 — mais luminal: ESR1 vermelho, KRT5 azul) para a direita (Q4 — mais basal-like: ESR1 azul, KRT5 vermelho).
- Subtipos Basal-like concentrados à direita; Luminal A à esquerda.
- O gradiente é contínuo — não há um "ponto de corte" visível, mas sim uma transição suave, confirmando a natureza contínua do escore.

### Heatmap METABRIC
- **N = 1.978 amostras** — o maior heatmap do conjunto.
- A proporção maior de casos ER+ (Luminal) pode tornar o gradiente menos dramático visualmente — mais amostras comprimidas no lado esquerdo.
- A presença de microarray Illumina (em vez de RNA-seq) pode criar padrões de intensidade ligeiramente diferentes, mas a direção e o padrão geral devem ser similares.
- Casos de alto escore (à direita) devem ter pior prognóstico DSS, visível se a anotação de evento estiver incluída.

### Heatmap GSE20685
- **N = 327 amostras** — o heatmap mais estreito.
- As colunas são mais largas visualmente por causa do N menor.
- Padrão similar às outras coortes, mas com menos resolução visual por causa do tamanho.
- A plataforma Affymetrix HGU133A pode criar artefatos visuais se a normalização não foi realizada corretamente — a ausência de artefatos é uma confirmação de qualidade.

## O que isso significa?

Os três heatmaps juntos servem como **validação visual** de que:

1. **O Z-score intra-coorte funcionou corretamente:** Cada gene tem média ≈ 0 e variância padronizada em cada coorte — visível como a ausência de uma coorte toda vermelha ou toda azul.

2. **O padrão biológico é consistente entre plataformas:** O gradiente de expressão dos genes CorePAM (basal vs. luminal) é visualmente similar em RNA-seq (TCGA-BRCA), microarray Illumina (METABRIC), e microarray Affymetrix (GSE20685). Isso é notável e confirma que os genes selecionados capturam um sinal biológico robusto, não tecnológico.

3. **A ordenação pelo escore funciona:** Ao ordenar as amostras pelo escore CorePAM, os genes com peso positivo aparecem em gradiente vermelho→azul (da direita para esquerda) e os genes com peso negativo em gradiente azul→vermelho. Isso confirma que a fórmula do escore está correta e que os pesos têm a direção biológica esperada.

4. **Subtipos moleculares se alinham com o escore:** Anotações de ER status e subtipo PAM50 devem se alinhar com a ordenação: casos ER− e Basal-like à direita (alto escore), casos ER+ e Luminal A à esquerda (baixo escore).

5. **Heterogeneidade dentro dos grupos:** Mesmo dentro do grupo "ER+" ou "Luminal B", há variação no escore — isso é visível no heatmap como amostras ER+ espalhadas por toda a faixa de escore baixo-a-médio, com algumas até no lado direito. Essa heterogeneidade intra-grupo é o que dá ao CorePAM seu valor clínico adicional além do ER.

## Contexto no estudo

Os heatmaps são figuras de "narrativa biológica" — eles transformam estatísticas em padrões visuais que qualquer oncologista ou pesquisador pode reconhecer. A separação luminal/basal no espaço dos genes CorePAM é literalmente visível a olho nu nessas figuras.

Para a defesa de tese, os heatmaps são poderosos porque demonstram visualmente que o CorePAM está capturando o mesmo sinal biológico em três coortes independentes, em duas plataformas diferentes, em populações de diferentes países. Essa concordância visual é a demonstração mais intuitiva de robustez cross-plataforma.

Do ponto de vista técnico, os heatmaps também confirmam que a normalização por Z-score intra-coorte foi executada corretamente — se houvesse erros de normalização, os padrões seriam caóticos ou dominados por efeitos de batch.

## Pontos de atenção

- **Heatmaps com muitas amostras** (TCGA: 1.072, METABRIC: 1.978) podem ser difíceis de visualizar em formato impresso. A versão digital/interativa permite zoom. Para publicação em papel, versões "condensadas" (por decil ou amostragem representativa) podem ser necessárias.

- **A escala de cores** é um parâmetro crítico: se os limites do gradiente de cores forem muito estreitos, padrões sutis ficam invisíveis; se muito largos, diferenças reais ficam saturadas. Tipicamente, os limites são definidos por percentis (ex: −2 a +2 DP) para evitar saturação por outliers.

- **Clustering hierárquico de genes** (dendrograma): Se os genes forem agrupados por similaridade de expressão (em vez de ordenados por peso), o padrão visual muda. O clustering pode revelar grupos de co-expressão biologicamente interessantes mas dificulta a interpretação direta dos pesos.

- **O heatmap ordenado por escore** (em vez de clustering hierárquico de amostras) é a representação mais direta para o objetivo deste estudo — mostrar que o escore captura um gradiente contínuo. Clustering hierárquico de amostras seria mais adequado para análises exploratórias de subtipos.

- Os heatmaps **não têm poder estatístico formal** — são visualizações. A significância estatística dos resultados vem das análises Cox e do C-index apresentados nas outras figuras. Os heatmaps complementam com narrativa visual.

- A **heterogeneidade aparente** dentro de grupos (ex: algumas amostras ER+ com padrão "basal-like" no heatmap) pode refletir erros de classificação clínica, subtipos mistos (ER+/HER2-enriched), ou genuína heterogeneidade molecular — essas amostras são biologicamente interessantes para investigação futura.
