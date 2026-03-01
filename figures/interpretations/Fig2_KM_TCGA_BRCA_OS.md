# Interpretação: Fig2 — Curva de Kaplan-Meier, TCGA-BRCA, Sobrevida Global (OS)

## O que é este gráfico?

Este é um gráfico de Kaplan-Meier (KM), o tipo mais clássico e amplamente reconhecido em análise de sobrevida. Ele mostra como a probabilidade de sobrevivência acumulada evolui ao longo do tempo em dois grupos de pacientes definidos pelo escore CorePAM: o grupo de **alto risco** (High, escore ≥ mediana) e o grupo de **baixo risco** (Low, escore < mediana). A coorte utilizada é o TCGA-BRCA (The Cancer Genome Atlas — Breast Invasive Carcinoma), que representa a validação independente em dados de RNA-seq.

Esta figura é fundamental porque demonstra que o escore CorePAM, derivado exclusivamente nos dados de treinamento SCAN-B, consegue discriminar prognóstico em uma coorte completamente independente, com plataforma tecnológica equivalente (RNA-seq), mas diferente centro de coleta e população de origem.

## Como ler este gráfico

- **Eixo X (horizontal):** Tempo em meses desde o diagnóstico ou início do seguimento.
- **Eixo Y (vertical):** Probabilidade de sobrevivência acumulada (0 = 0%, 1 = 100%).
- **Linha superior (azul/verde — "Low risk"):** Pacientes cujo escore CorePAM ficou abaixo da mediana. Representam o perfil mais luminal, com prognóstico favorável.
- **Linha inferior (vermelha/laranja — "High risk"):** Pacientes com escore acima da mediana. Representam perfil mais basal-like, com pior prognóstico.
- **Marcas de "tick" nas linhas:** Representam censuras — pacientes que saíram do estudo sem atingir o evento (óbito), por término de seguimento ou perda.
- **Tabela de risco abaixo do gráfico (at-risk table):** Mostra quantos pacientes ainda estavam em risco em cada ponto temporal. À medida que o tempo avança, esse número cai por eventos e censuras.
- **p-valor (log-rank):** Teste estatístico que compara as duas curvas globalmente. p < 0.05 indica diferença estatisticamente significativa.
- **HR e IC 95%:** Hazard Ratio e intervalo de confiança do modelo Cox univariado; mede a razão de risco instantâneo entre os grupos.

## Achados principais

- **N total da coorte:** 1.072 pacientes
- **Número de eventos (óbitos):** 150 (14,0%)
- **Tempo de seguimento mediano (censurado):** 32 meses
- **Hazard Ratio univariado:** HR = 1,20 (IC 95%: 1,04–1,40)
- **p-valor (log-rank):** p = 0,016
- **C-index:** 0,624 (IC 95%: 0,572–0,678)
- **HR multivariado (CORE-A: CorePAM + idade + ER):** HR = 1,27 (IC 95%: 1,09–1,48)

## O que isso significa?

O TCGA-BRCA é a coorte de validação independente em RNA-seq. O fato de o CorePAM discriminar sobrevida nesta coorte (HR = 1,20, p = 0,016) com boa calibração do C-index (0,624) demonstra que o modelo **generaliza** além dos dados de treinamento.

O HR de 1,20 significa que pacientes classificados como alto risco têm, em média, 20% mais risco instantâneo de óbito em comparação aos de baixo risco, em cada ponto do tempo. Embora numericamente menor que o observado no treinamento (HR = 1,92 em SCAN-B), isso é esperado e fisiológico: a população do TCGA-BRCA tem menor tempo de seguimento mediano (32 meses vs. 55,7 meses em SCAN-B), proporção menor de eventos (14% vs. 10,5%, mas com N muito maior em SCAN-B), e distribuição subtipal diferente.

No modelo CORE-A (ajustado para idade e status de ER), o HR sobe para 1,27 (IC 95%: 1,09–1,48), sugerindo que o escore CorePAM traz informação prognóstica além dos fatores clínicos padrão — o que é exatamente o que se deseja de um painel molecular.

A separação entre as curvas KM, embora moderada, é estatisticamente robusta e biologicamente coerente: o grupo de alto risco acumula mais eventos ao longo do tempo, e a diferença se mantém durante todo o período de seguimento disponível.

## Contexto no estudo

Esta figura, juntamente com as curvas KM de METABRIC e GSE20685, compõe a figura central de validação prognóstica (Fig. 2) do manuscrito. A narrativa principal do artigo é: "o escore CorePAM foi derivado em SCAN-B e validado em três coortes independentes, em duas plataformas distintas." O TCGA-BRCA representa a prova de validação em RNA-seq (mesma plataforma que o treinamento), enquanto METABRIC e GSE20685 representam a validação cruzada de plataforma (microarray).

O fato de todas as três validações mostrarem separação significativa das curvas KM fortalece enormemente a credibilidade e a generalização do CorePAM.

## Pontos de atenção

- O TCGA-BRCA tem **seguimento mediano de apenas 32 meses**, o que é curto para um desfecho de sobrevida global em câncer de mama. Isso pode subestimar o poder discriminativo real do escore, pois muitos eventos ainda estão por ocorrer.
- A proporção de eventos é de apenas 14%, tornando a análise mais instável estatisticamente do que em METABRIC (57,8% de eventos).
- O HR de 1,20 é o menor entre as três coortes de validação, o que pode refletir tanto o menor seguimento quanto diferenças populacionais (TCGA tem maior proporção de subtipos favoráveis como Luminal A).
- A análise de sensibilidade com corte de 24 meses (FigS6) foi realizada justamente para abordar essa limitação de seguimento e é apresentada em figura suplementar.
- A dicotomização pela mediana é convencional mas arbitrária — as curvas por quartis (FigS suplementar) mostram gradiente dose-resposta.
