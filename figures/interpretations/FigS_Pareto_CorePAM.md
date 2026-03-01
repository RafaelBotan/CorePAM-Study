# Interpretação: FigS — Fronteira de Pareto: C-index OOF vs. Número de Genes (df)

## O que é este gráfico?

Esta figura suplementar apresenta a **fronteira de Pareto** gerada durante a derivação do CorePAM — o coração metodológico de toda a redução do PAM50. É uma das figuras mais importantes do ponto de vista metodológico, pois mostra **como e por que** foram escolhidos exatamente esses 24 genes.

O gráfico mostra uma relação entre:
- **Eixo Y:** C-index OOF (Out-of-Fold) — a capacidade discriminativa do modelo de Cox elastic-net estimada por validação cruzada de 10 folds em SCAN-B
- **Eixo X:** df (degrees of freedom) — que no contexto do elastic-net com penalização L1/L2 representa o número efetivo de genes com pesos não-nulos no modelo

Cada ponto no gráfico representa um modelo com uma configuração de regularização λ específica — à medida que λ diminui (menos regularização), mais genes entram no modelo (df maior) e o C-index tende a aumentar. A pergunta crucial é: **a partir de qual df o ganho marginal em C-index se torna negligível?**

## Como ler este gráfico

- **Cada ponto:** Um modelo de Cox elastic-net com λ específico, avaliado por validação cruzada de 10 folds em SCAN-B.
- **Curva crescente:** O C-index geralmente sobe à medida que mais genes são incluídos (df cresce), mas com **retornos decrescentes** — a melhora marginal por gene adicional diminui.
- **Linha horizontal tracejada:** Marca o C-index do modelo de referência PAM50 completo (todos os 50 genes). Esta é a "linha de comparação" — modelos abaixo dela são inferiores ao PAM50.
- **Linha horizontal de não-inferioridade:** Marca o C-index PAM50 − ΔC (onde ΔC = 0,010 é o critério de não-inferioridade pré-definido nos parâmetros congelados). Qualquer modelo com C-index acima dessa linha é "não-inferior ao PAM50."
- **Ponto marcado (selecionado):** O ponto correspondente ao modelo CorePAM selecionado — o **menor df** cujo C-index OOF está acima da linha de não-inferioridade.
- **Barra de erro em cada ponto:** IC 95% do C-index OOF (variação entre os 10 folds).

## Achados principais

- **C-index PAM50 completo (referência):** Valor estimado em SCAN-B (OOF)
- **Critério de não-inferioridade:** ΔC = 0,010 (parâmetro congelado desde o protocolo)
- **df selecionado:** Corresponde ao CorePAM de 24 genes
- **C-index OOF do CorePAM:** C_PAM50 − ΔC ≤ C_CorePAM ≤ C_PAM50 (dentro da zona de não-inferioridade)
- **Genes eliminados:** 26 genes do PAM50 original com peso zero no modelo regularizado

A fronteira de Pareto mostra que com 24 genes (o df selecionado), o C-index OOF está dentro da margem de não-inferioridade definida. Adicionar mais genes traz ganhos marginais de C-index inferiores a 0,010 — portanto, os genes extras não justificam a complexidade adicional do painel.

## O que isso significa?

Este gráfico é a prova de conceito da redução: **24 genes são suficientes para capturar o sinal prognóstico dos 50 genes do PAM50**, dentro da margem de erro aceitável de 1 ponto percentual de C-index.

O critério de não-inferioridade (ΔC = 0,010) é o parâmetro mais importante do estudo. Ele representa a resposta à pergunta: "quanto pior o modelo menor pode ser comparado ao PAM50 completo antes de considerarmos a redução inadequada?" A escolha de 0,010 é conservadora e clínica: 1 ponto percentual de C-index é geralmente considerado abaixo da diferença mínima clinicamente relevante em estudos de prognóstico de câncer de mama.

A fronteira de Pareto também ilustra um princípio fundamental de modelagem preditiva: a **navalha de Occam estatística** — modelos mais simples que são estatisticamente equivalentes a modelos mais complexos devem ser preferidos por serem mais robustos, interpretáveis e generalizáveis.

O número de genes (24) **não foi pré-especificado** — emergiu do critério de não-inferioridade aplicado à fronteira Pareto. Isso é uma diferença filosófica importante em relação ao artigo anterior (que usava 40 genes por critério de cobertura de plataforma): aqui, o tamanho do painel é **determinado pelos dados e pelo critério estatístico**, não por uma escolha subjetiva prévia.

## Contexto no estudo

Esta figura suplementar é a justificativa metodológica para todo o estudo. Sem ela, o leitor poderia perguntar: "por que 24 genes? por que não 20 ou 30?" O gráfico de Pareto responde visualmente: **24 é o menor número que satisfaz o critério de não-inferioridade pré-definido.**

Para a tese de doutorado, este gráfico é um dos mais importantes — ele representa a inovação metodológica central (derivação statistically-principled do tamanho do painel) e é o que diferencia este trabalho de abordagens anteriores que escolhiam genes por critérios qualitativos.

## Pontos de atenção

- O C-index OOF é uma estimativa da performance esperada em dados externos, mas ainda pode ser otimista se houver algum grau de overfit na seleção de λ — a validação real está nas coortes externas (TCGA, METABRIC, GSE20685).
- A fronteira de Pareto foi gerada exclusivamente em SCAN-B — portanto, o "24 genes" reflete a informação disponível nessa coorte específica. Em outra coorte de treinamento, o número poderia ser diferente.
- O critério ΔC = 0,010 é razoável mas não único — outros valores (0,005 ou 0,020) levariam a painéis de tamanho diferente. A escolha de 0,010 foi pré-especificada no `analysis_freeze.csv` antes da análise.
- O elastic-net (α = 0,5) combina penalização L1 (shrinkage zero) e L2 (shrinkage contínuo) — a escolha de α influencia quais genes são zerados e quais retêm pesos. Com α = 1 (Lasso puro), possivelmente menos genes seriam selecionados; com α = 0 (Ridge), nenhum seria zerado.
- A interpretação da "fronteira" de Pareto pressupõe que o caminho de regularização (λ sequence) cobre adequadamente o espaço de df de 1 a 50 — verificar que o gráfico inclui pontos em toda essa faixa.
