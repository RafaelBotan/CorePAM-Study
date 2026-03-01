# Interpretação: FigS — Escore CorePAM por Status de ER, por Coorte

## O que é este gráfico?

Este gráfico suplementar mostra a **distribuição do escore CorePAM estratificada pelo status de ER (receptor de estrogênio)** — positivo (ER+) versus negativo (ER−) — em cada coorte. Tipicamente apresentado como boxplots, violin plots, ou histogramas sobrepostos, com um painel por coorte.

O status de ER é o marcador clínico mais importante em câncer de mama — divide os tumores em dois grandes grupos biológicos e define o tratamento. Tumores ER+ são tipicamente Luminais (A ou B) e têm melhor prognóstico; tumores ER− incluem Basal-like (triplo-negativos) e HER2-enriched, com comportamento mais agressivo.

## Como ler este gráfico

- **Cada painel:** Uma coorte (SCAN-B, TCGA-BRCA, METABRIC, GSE20685).
- **Eixo X:** Grupos (ER+ vs. ER−, ou vice-versa).
- **Eixo Y:** Valor do escore CorePAM.
- **Boxplot:** Mediana (linha central), IQR (caixa), whiskers (1,5×IQR), e outliers individuais.
- **Violino (se usado):** Mostra a distribuição completa além do boxplot — um violino largo indica muitos valores naquele ponto.
- **Separação esperada:** Escores CorePAM **mais altos em ER−** (porque ER− inclui mais tumores basal-like, que expressam mais EXO1, KRT5, CENPF, etc.) e **mais baixos em ER+** (que expressam mais ESR1, PGR, BCL2, MLPH, etc.).

## Achados principais

O resultado esperado e biologicamente necessário:

- **ER− tem escore CorePAM significativamente maior que ER+** em todas as coortes.
- A diferença entre os grupos (ER+ vs. ER−) deve ser **consistente entre coortes** — mesmo em plataformas diferentes (RNA-seq vs. microarray).
- O **METABRIC** pode mostrar maior sobreposição entre ER+ e ER− porque a classificação de ER era menos padronizada em eras mais antigas (corte diferente, método IHC diferente).
- **GSE20685** sendo sul-coreana pode ter distribuição ligeiramente diferente por proporção de subtipos.

Estatisticamente, espera-se:
- Teste t ou Wilcoxon com p < 0,001 em todas as coortes (exceto possivelmente GSE20685 por N menor).
- Diferença de mediana entre ER+ e ER− de ~0,5 a 1,5 pontos na escala do escore.

## O que isso significa?

Esta análise valida a **direção biológica** do escore CorePAM. Se o escore fosse invertido (ER+ com escore alto e ER− com escore baixo), indicaria um erro na direção dos pesos. A associação ER− ↔ escore alto é a "verificação de sanidade" biológica do modelo.

Mais profundamente, esta figura ilustra **por que** o escore CorePAM funciona:
- Os genes com **peso positivo** (EXO1, PTTG1, FOXC1, PHGDH, MYC, KRT5, CXXC5, CENPF, FGFR4, ERBB2) são mais expressos em tumores ER− (Basal-like e HER2-enriched).
- Os genes com **peso negativo** (NAT1, BLVRA, ACTR3B, MIA, MYBL2, MDM2, SFRP1, GPR160, PGR, KRT17, BCL2, MLPH, GRB7, ESR1) são mais expressos em tumores ER+ (Luminais).
- Portanto, tumores ER− acumulam mais termos positivos e menos negativos → escore mais alto → maior risco previsto. Tumores ER+ fazem o oposto.

A separação do escore por ER confirma que o painel CorePAM é, de certa forma, um **proxy molecular do status de ER** — mas não completamente. Dentro de cada grupo (ER+ alto vs. ER+ baixo escore; ER− alto vs. ER− baixo escore), ainda há variação do escore que é prognóstica. Isso é o que torna o CorePAM útil além do ER clínico.

## Contexto no estudo

Esta figura é a ponte entre a biologia molecular e a clínica. O status de ER é determinado por imunohistoquímica (IHC) na prática clínica — um teste simples, barato e universalmente disponível. O escore CorePAM, por outro lado, requer medição de expressão gênica de 24 genes. Esta figura mostra que:

1. **O CorePAM captura o sinal de ER** (confirmado pela separação dos grupos).
2. **O CorePAM vai além do ER** (há variação dentro de cada grupo de ER, e é justamente essa variação que fornece o valor incremental visto no ΔC-index).

O modelo CORE-A (que inclui ER como covariável) e o ΔC-index elevado mesmo após ajuste por ER confirmam esse ponto.

## Pontos de atenção

- A separação entre ER+ e ER− pelo escore CorePAM confirma a direção biológica, mas também levanta uma questão: **o CorePAM não é apenas "outro teste de ER mais caro?"** A resposta está no ΔC-index dentro de subgrupos ER+ e ER−, que mostraria que o CorePAM agrega informação dentro de cada grupo separadamente.
- A **qualidade do status de ER** varia entre coortes. Em SCAN-B (coorte moderna), ER é determinado por IHC com cortes padronizados. Em METABRIC e GSE20685 (coortes históricas), os critérios podem ter variado ao longo do tempo.
- A **proporção de ER+ e ER−** varia entre coortes e populações. SCAN-B tem proporção típica europeia (~75% ER+); GSE20685 (sul-coreana) pode ter proporção ligeiramente diferente.
- Os genes do CorePAM com peso negativo incluem **ESR1** (o gene para receptor de estrogênio) — portanto, o escore "captura" parcialmente o ER molecular, não apenas o ER clínico/IHC. Isso pode criar alguma circularidade se ER molecular e clínico são usados juntos, mas é biologicamente coerente.
