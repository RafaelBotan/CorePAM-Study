# Interpretação: Fig2 — Curva de Kaplan-Meier, METABRIC, Sobrevida Específica por Doença (DSS)

## O que é este gráfico?

Este gráfico de Kaplan-Meier mostra a probabilidade de sobrevivência acumulada ao longo do tempo na coorte METABRIC (Molecular Taxonomy of Breast Cancer International Consortium), estratificada em dois grupos pelo escore CorePAM: **alto risco** (High, escore ≥ mediana) e **baixo risco** (Low, escore < mediana).

Diferentemente das outras coortes, o **endpoint primário do METABRIC é DSS (Disease-Specific Survival)**, ou seja, sobrevida específica por doença — contabilizando apenas óbitos causados pelo câncer de mama como eventos, e censurando óbitos por outras causas. Essa escolha metodológica é deliberada e justificada: o METABRIC possui seguimento longo (mediana de 159 meses — mais de 13 anos!), e em períodos tão extensos, mortes por causas não relacionadas ao câncer representariam ruído considerável se incluídas como eventos. O DSS é, portanto, o endpoint mais informativo e metodologicamente adequado para esta coorte.

Esta figura representa a validação mais robusta do CorePAM no estudo, tanto pelo tamanho amostral (N=1.978) quanto pelo altíssimo número de eventos (646) e pelo longo período de seguimento.

## Como ler este gráfico

- **Eixo X (horizontal):** Tempo em meses. Note que o eixo se estende até mais de 200 meses (aproximadamente 17 anos) — isso reflete a riqueza única do seguimento METABRIC.
- **Eixo Y (vertical):** Probabilidade de sobrevivência específica por doença (0 a 1).
- **Linha superior (azul — "Low risk"):** Pacientes com escore CorePAM abaixo da mediana; perfil mais luminal, prognóstico mais favorável.
- **Linha inferior (vermelha — "High risk"):** Pacientes com escore acima da mediana; perfil mais basal-like, evolução mais agressiva.
- **Ticks nas curvas:** Pacientes censurados (saíram do estudo sem evento relacionado ao câncer).
- **Tabela de risco:** Fundamental aqui — ao longo de mais de 15 anos, o número de pacientes em risco reduz substancialmente. Observe como, mesmo em pontos tardios (>150 meses), ainda há pacientes sendo monitorados.
- **p-valor (log-rank) e HR:** Estatísticas globais comparando as duas curvas.

## Achados principais

- **N total:** 1.978 pacientes
- **Eventos DSS:** 646 (32,6%) — proporção de eventos muito superior às demais coortes
- **Tempo de seguimento mediano (censurado):** 159 meses (aproximadamente 13,2 anos)
- **Hazard Ratio univariado (DSS):** HR = 1,41 (IC 95%: 1,31–1,52)
- **p-valor:** p = 3,45 × 10⁻¹⁹ (altamente significativo)
- **C-index:** 0,638 (IC 95%: 0,615–0,659)
- **Delta C-index (CORE-A vs. CORE-A + CorePAM):** ΔC = +0,142 (IC 95%: 0,097–0,163) — o maior ganho incremental do estudo
- **Sensibilidade OS (análise secundária):** HR_OS = 1,23 (IC 95%: 1,16–1,30)

## O que isso significa?

O METABRIC é, sem dúvida, a validação mais poderosa deste estudo. Com 646 eventos e quase 160 meses de seguimento, a precisão estatística é excepcional. Um HR de 1,41 (p = 3,45 × 10⁻¹⁹) é extremamente significativo — a probabilidade de esse resultado ser acaso é de aproximadamente 0,000000000000000000345 (1 em 2,9 × 10¹⁸).

O que torna o resultado ainda mais impressionante é que o METABRIC usa **microarray Illumina**, enquanto o CorePAM foi derivado em **RNA-seq** (SCAN-B). A transferência do escore de uma plataforma para outra sem qualquer recalibração dos coeficientes — apenas com Z-score intra-coorte — e ainda assim obter HR de 1,41 com separação clara das curvas KM demonstra uma robustez tecnológica notável do painel.

O ΔC-index de +0,142 no METABRIC é o maior ganho incremental observado em todo o estudo: adicionar o escore CorePAM ao modelo clínico (CORE-A = idade + ER) melhora a capacidade discriminativa em 14,2 pontos percentuais. Isso sugere que o CorePAM captura informação molecular que simplesmente não está disponível nas variáveis clínicas convencionais — especialmente relevante em uma coorte tratada em era pré-molecular como o METABRIC.

As curvas KM mostram separação precoce (já nos primeiros 24–36 meses) e essa separação se mantém ao longo de toda a janela de seguimento, inclusive além de 150 meses. Esse perfil de separação sustentada indica que o escore não apenas identifica risco de curto prazo, mas captura diferenças prognósticas de longo prazo — o que é exatamente o que se espera de um classificador baseado em subtipo molecular.

## Contexto no estudo

O METABRIC tem papel duplo neste estudo: é a validação microarray principal (com DSS como endpoint primário) e a coorte que fornece o maior número de eventos, ancorando a estabilidade da meta-análise. A análise de sensibilidade com OS (em vez de DSS) é apresentada na FigS5 como validação da escolha do endpoint.

A separação tão clara das curvas KM no METABRIC, contrastando com o HR mais modesto do TCGA-BRCA, reflete a complementaridade das coortes: o TCGA-BRCA é jovem (curto FU), enquanto o METABRIC é maduro (longo FU). Juntos, eles cobrem espectros temporais complementares do prognóstico de mama.

## Pontos de atenção

- O METABRIC foi coletado em era pré-trastuzumabe e pré-CDK4/6 inibidores para a maioria dos casos, o que pode superestimar o impacto prognóstico do escore em comparação com populações tratadas com terapias modernas.
- A escolha do **DSS como endpoint primário** (em vez de OS) é justificada pelo longo seguimento, mas é importante notar que 57,8% dos pacientes morreram por qualquer causa — a diferença entre OS e DSS é relevante e explorada na análise de sensibilidade.
- A plataforma microarray Illumina captura menos genes por sonda que o RNA-seq, mas todos os 24 genes do CorePAM estavam disponíveis no METABRIC — o que não acontece com o GSE20685 (HGU133A, 44/50 genes PAM50 disponíveis).
- A separação das curvas KM por mediana (50/50) é conveniente para visualização, mas as curvas por quartis (FigS) mostram que há gradiente contínuo — o escore não é simplesmente binário.
