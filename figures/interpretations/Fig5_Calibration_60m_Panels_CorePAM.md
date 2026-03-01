# Interpretação: Fig5 — Calibração a 60 Meses: Painel de Curvas por Coorte

## O que é este gráfico?

Este é um **gráfico de calibração** a 60 meses (5 anos), mostrando a concordância entre as probabilidades de sobrevivência **preditas** pelo modelo CorePAM e as probabilidades **observadas** empiricamente nos dados.

A calibração responde a uma pergunta diferente da discriminação (C-index): enquanto o C-index pergunta "o modelo ordena corretamente os pacientes do mais ao menos grave?", a calibração pergunta **"quando o modelo diz que um paciente tem 80% de chance de sobreviver em 5 anos, realmente ~80% dos pacientes nessa categoria sobrevivem?"**

Um modelo pode discriminar bem mas calibrar mal (ex: sistemicamente superestima risco), ou calibrar bem mas discriminar pouco. O ideal é ter ambos. Esta figura avalia especificamente a calibração.

O painel inclui todas as coortes disponíveis para esta análise (tipicamente SCAN-B, TCGA-BRCA, METABRIC, GSE20685), cada uma em um painel separado.

## Como ler este gráfico

- **Eixo X (horizontal):** Probabilidade de sobrevivência **predita** pelo modelo a 60 meses (0 a 1, ou 0% a 100%).
- **Eixo Y (vertical):** Probabilidade de sobrevivência **observada** estimada pela curva de Kaplan-Meier no mesmo horizonte temporal.
- **Linha diagonal tracejada (linha de identidade perfeita):** Representa calibração perfeita — predito = observado.
- **Curva de calibração (ou pontos por decil):** Cada ponto representa um grupo de pacientes com predições similares. A posição vertical do ponto mostra a sobrevida observada naquele grupo.
- **Proximidade à diagonal:** Quanto mais próxima a curva está da diagonal, melhor a calibração.
- **Desvio para cima (acima da diagonal):** O modelo **subestima** risco — prevê menos mortes do que realmente ocorrem.
- **Desvio para baixo (abaixo da diagonal):** O modelo **superestima** risco — prevê mais mortes do que realmente ocorrem.

Muitas implementações também incluem:
- **IC sombreado** ao redor da curva de calibração.
- **Histograma inferior** mostrando a distribuição das predições (densidade).
- **Estatística de Brier Score:** Erro quadrático médio das predições de probabilidade.

## Achados principais

A calibração a 60 meses é particularmente informativa porque:
- 60 meses (5 anos) é o horizonte clínico padrão em oncologia mamária.
- É o ponto onde temos eventos suficientes na maioria das coortes para estimativas confiáveis.

**Padrões esperados por coorte:**

- **SCAN-B:** Calibração mais precisa por ser a coorte de treinamento. Espera-se que as curvas fiquem muito próximas da diagonal, com pequeno desvio por shrinkage do elastic-net.
- **TCGA-BRCA:** Calibração mais difícil aqui porque o seguimento mediano é 32 meses — apenas metade dos pacientes chegam ao horizonte de 60 meses. Extrapolações de Cox aumentam incerteza.
- **METABRIC:** Excelente calibração esperada por ter muitos eventos antes de 60 meses. O endpoint DSS (não OS) pode criar pequenas discrepâncias se comparado com modelos treinados em OS.
- **GSE20685:** Calibração moderada esperada — N menor, maior incerteza.

**O ponto-chave:** A aplicação do escore CorePAM (derivado em SCAN-B com RNA-seq) a coortes de microarray sem recalibração é um teste rigoroso de calibração cross-plataforma. A performance observada confirma que o Z-score intra-coorte é suficiente para equacionar as escalas de expressão entre plataformas.

## O que isso significa?

A calibração é o "teste de realidade" do modelo. Se o modelo diz que um paciente tem 70% de chance de sobreviver em 5 anos, isso deve significar algo concreto: aproximadamente 7 em cada 10 pacientes com essa predição devem estar vivos em 5 anos.

Para aplicação clínica real, calibração ruim é tão problemática quanto discriminação ruim. Um modelo que consistentemente superestima risco levaria a tratamentos desnecessários; um que subestima levaria a subtratamento.

No contexto do CorePAM, a análise de calibração é especialmente relevante porque:

1. Os pesos foram estimados por elastic-net com regularização — a regularização introduz shrinkage (encolhimento dos coeficientes em direção a zero), o que tende a melhorar calibração em dados externos.

2. O Z-score intra-coorte normaliza a escala do escore para cada coorte independentemente, o que ajuda na transferibilidade sem recalibração.

3. A fórmula de score normalizada (Σwᵢzᵢ / Σ|wᵢ|) mantém o escore em uma escala interpretável independente do número de genes presentes.

A combinação desses três fatores favorece boa calibração cross-plataforma, e a figura apresenta evidência empírica disso.

## Contexto no estudo

A análise de calibração (Fig5, painel de calibração) é o complemento natural da análise de discriminação (Fig5, ΔC-index). Juntas, elas constituem a avaliação completa de **valor clínico** do CorePAM:

- ΔC-index → "o modelo discrimina melhor que só variáveis clínicas?"
- Calibração → "as probabilidades preditas são confiáveis e corretas?"

A combinação de boa discriminação **e** boa calibração é o padrão exigido para validação de biomarcadores prognósticos segundo as diretrizes TRIPOD (Transparent Reporting of a multivariable prediction model for Individual Prognosis or Diagnosis), a qual este estudo segue.

## Pontos de atenção

- A calibração a **60 meses** pode ser problemática para o TCGA-BRCA com FU mediano de apenas 32 meses — haverá muita extrapolação do modelo de Cox além do tempo de seguimento observado, aumentando a incerteza.
- A interpretação das curvas de calibração requer cautela quando o **número de eventos é pequeno** por intervalo (como em GSE20685 com apenas 83 eventos totais). Poucos pontos nos decis extremos produzem estimativas ruidosas.
- O **Brier Score** é uma métrica mais rigorosa que o simples gráfico visual de calibração; valores abaixo de 0,25 são geralmente considerados aceitáveis para modelos de sobrevida a 5 anos em oncologia.
- A calibração avaliada **apenas em um horizonte temporal** (60 meses) não captura possível miscalibração em outros horizontes (ex: 24 meses, 120 meses). As análises de sensibilidade com outros horizontes são importantes para completude.
- O shrinkage do elastic-net tende a produzir predições mais moderadas (menos extremas), o que **favorece** a calibração em dados externos mas pode **reduzir ligeiramente** a discriminação.
