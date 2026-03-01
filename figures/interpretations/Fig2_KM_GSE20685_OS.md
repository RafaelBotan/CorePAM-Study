# Interpretação: Fig2 — Curva de Kaplan-Meier, GSE20685, Sobrevida Global (OS)

## O que é este gráfico?

Este gráfico de Kaplan-Meier mostra a probabilidade de sobrevivência global (Overall Survival, OS) ao longo do tempo na coorte GSE20685 (depositada no GEO — Gene Expression Omnibus), estratificada pelo escore CorePAM: **alto risco** (High, escore ≥ mediana) versus **baixo risco** (Low, escore < mediana).

O GSE20685 é a segunda coorte de validação em microarray, mas usa uma plataforma completamente diferente do METABRIC: o microarray **Affymetrix HGU133A**. Esta plataforma tem cobertura de genes mais restrita, e por isso apenas 44 dos 50 genes originais do PAM50 estão disponíveis nela (os 6 faltantes: ANLN, CXXC5, GPR160, NUF2, TMEM45B, UBE2T). Destes 44 disponíveis, todos os 24 genes do CorePAM que precisam estar presentes o estão — o que valida a aplicabilidade do painel nesta plataforma.

Esta coorte é particularmente interessante por sua origem sul-coreana (coletada em Seul), representando uma população asiática — o que adiciona diversidade étnica e geográfica à validação.

## Como ler este gráfico

- **Eixo X (horizontal):** Tempo em meses. O eixo chega a ~200 meses (aproximadamente 16 anos), refletindo o longo seguimento desta coorte.
- **Eixo Y (vertical):** Probabilidade de sobrevivência global acumulada (0 a 1).
- **Linha superior (azul — "Low risk"):** Pacientes com escore CorePAM abaixo da mediana; perfil luminal, melhor sobrevida.
- **Linha inferior (vermelha — "High risk"):** Pacientes com escore acima da mediana; perfil mais agressivo, pior sobrevida.
- **Ticks nas curvas:** Censuras — pacientes sem evento ao final do seguimento.
- **Tabela de risco:** A coorte começa com N=327, e o número cai progressivamente. Ao redor de 100 meses, ainda há algumas dezenas em risco em cada grupo.
- **IC sombreado:** A área sombreada ao redor de cada curva representa o intervalo de confiança pontual de Greenwood — quanto mais estreito, mais precisa é a estimativa naquele ponto.

## Achados principais

- **N total:** 327 pacientes
- **Eventos (óbitos):** 83 (25,4%)
- **Tempo de seguimento mediano (censurado):** 112,8 meses (aproximadamente 9,4 anos)
- **Hazard Ratio univariado:** HR = 1,40 (IC 95%: 1,12–1,76)
- **p-valor (log-rank):** p = 0,003
- **C-index:** 0,623 (IC 95%: 0,562–0,683)
- **Delta C-index (CORE-A → CORE-A + CorePAM):** ΔC = +0,093 (IC 95%: 0,021–0,175)

## O que isso significa?

O GSE20685 é a menor coorte de validação do estudo (N=327), mas ainda assim apresenta resultados estatisticamente significativos e clinicamente relevantes. Um HR de 1,40 (p = 0,003) com C-index de 0,623 demonstra que o CorePAM consegue discriminar prognóstico mesmo com:

1. **Plataforma completamente diferente** do treinamento (Affymetrix HGU133A vs. RNA-seq Illumina em SCAN-B),
2. **População de origem diferente** (sul-coreana vs. sueca e norte-americana),
3. **N menor**, que limita o poder estatístico.

O valor do HR (1,40) é comparável ao do METABRIC (1,41), o que é notável dada a diferença de tamanho entre as coortes. Isso sugere que, nas coortes com longo seguimento (FU mediano >100 meses em ambas), o efeito do CorePAM se manifesta com magnitude similar.

O ganho incremental de C-index (ΔC = +0,093) ao adicionar CorePAM ao modelo clínico CORE-A é substancial: em uma coorte de 327 pacientes com 83 eventos, um ganho de ~9 pontos percentuais em discriminação é clinicamente significativo.

A separação das curvas KM é visualmente clara e se inicia precocemente (primeiros 24–36 meses), mantendo-se ao longo de todo o período de seguimento disponível. O grupo de alto risco acumula eventos de forma consistentemente mais rápida.

## Contexto no estudo

O GSE20685 cumpre um papel único na estrutura de validação: é a evidência de que o CorePAM funciona em populações asiáticas e em plataformas de microarray mais antigas (HGU133A), amplamente usadas em estudos de câncer de mama da era 2005–2015. Muitos estudos clássicos de câncer de mama usaram exatamente esta plataforma (Affymetrix HGU133A/HGU133Plus2), de modo que a compatibilidade do CorePAM com ela é relevante para comparações históricas.

Junto com o METABRIC, o GSE20685 estabelece a validação cruzada de plataforma (cross-platform validation) — um dos critérios mais rigorosos para robustez de um painel molecular. A concordância de resultados entre RNA-seq (TCGA-BRCA) e dois microarrays diferentes (METABRIC e GSE20685) é a assinatura de um painel verdadeiramente robusto.

## Pontos de atenção

- O **N=327 é o menor do estudo**, e o intervalo de confiança do HR (1,12–1,76) é correspondentemente mais largo do que em METABRIC. A imprecisão estatística é maior.
- O **IC 95% do C-index** (0,562–0,683) é bastante amplo — um reflexo direto do tamanho amostral menor.
- O **ΔC de +0,093** tem IC amplo (0,021–0,175), sugerindo que, embora o ganho médio seja substancial, há incerteza considerável em torno dessa estimativa.
- A coorte GSE20685 representa uma população sul-coreana, que pode ter distribuição de subtipos moleculares ligeiramente diferente das populações ocidentais (menor frequência de Luminal B?). Isso deveria ser explorado em análises futuras.
- Os 6 genes PAM50 ausentes no HGU133A não fazem parte do CorePAM de 24 genes — o que confirma que a derivação do CorePAM foi, inadvertidamente, robusta à disponibilidade limitada desta plataforma. No entanto, isso deve ser reportado como limitação de cobertura da plataforma, não como design prospectivo.
