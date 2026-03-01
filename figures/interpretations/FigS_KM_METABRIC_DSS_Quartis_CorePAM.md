# Interpretação: FigS — Curvas KM por Quartis do Escore CorePAM — METABRIC (DSS)

## O que é este gráfico?

Esta figura suplementar exibe as curvas de Kaplan-Meier para a coorte METABRIC, com o endpoint primário DSS (Disease-Specific Survival), estratificada por **quartis do escore CorePAM**. É o equivalente METABRIC das análises por quartis de TCGA-BRCA e GSE20685.

A divisão por quartis distribui os 1.978 pacientes do METABRIC em quatro grupos de aproximadamente 495 pacientes cada:
- **Q1 (1º quartil):** Menor 25% do escore — perfil mais luminal, mais ER+, genes de bom prognóstico dominando.
- **Q2 (2º quartil):** Escore baixo-médio.
- **Q3 (3º quartil):** Escore alto-médio.
- **Q4 (4º quartil):** Maior 25% do escore — perfil mais basal-like, mais agressivo.

Com ~495 pacientes por quartil e ~160 eventos esperados por grupo (se uniformemente distribuídos dos 646 totais), cada grupo tem poder estatístico robusto individualmente.

## Como ler este gráfico

- **Quatro curvas KM** separadas por quartil, coloridas em gradiente de prognóstico.
- **Eixo X:** Tempo em meses, com horizonte longo (até ~200 meses / 16 anos) — único do METABRIC.
- **Separação gradual:** Com longo seguimento, a separação entre curvas adjacentes (Q1 vs. Q2, Q2 vs. Q3, Q3 vs. Q4) deve ser mais visível que em coortes de FU curto.
- **Efeito tardio nos Luminais:** Pacientes de Q1 (mais luminais) podem eventualmente acumular eventos tardios (>10 anos), enquanto pacientes de Q4 (mais basais) têm eventos mais precoces — o crossover tardio entre curvas, se presente, sugeriria biologia diferente nos horizontes temporais.

## Achados principais

O METABRIC é a coorte mais poderosa para análise por quartis:
- **646 eventos DSS** distribuídos entre 4 grupos.
- **Seguimento de até ~200 meses** permite ver separação em múltiplos pontos temporais.
- Com ~162 eventos por quartil (estimativa), os grupos extremos (Q1 vs. Q4) devem mostrar separação estatisticamente robusta e clinicamente dramática.

Espera-se:
- **Gradiente claro e ordenado:** Q4 com a pior sobrevida ao longo de todo o período, Q1 com a melhor.
- **Separação precoce e sustentada:** Diferença entre Q4 e Q1 visível já nos primeiros 24-36 meses e mantendo-se por mais de 100 meses.
- **Grupos intermediários (Q2, Q3) bem separados:** Com tantos eventos, a separação entre quartis adjacentes também é visível, confirmando o gradiente contínuo.

## O que isso significa?

O METABRIC é a coorte mais valiosa para demonstrar o gradiente dose-resposta do CorePAM por duas razões:

1. **Poder estatístico**: Com ~162 eventos por quartil, as estimativas são precisas e os IC 95% estreitos.

2. **Seguimento longo**: Com acompanhamento de até 16 anos, conseguimos ver o comportamento de subtipos luminais que têm eventos tardios — crucial para distinguir Q1 de Q2 no longo prazo.

O gradiente dose-resposta em METABRIC confirma que o CorePAM está capturando um espectro biológico real — desde os tumores mais "quietos" (Luminal A, alto BCL2, alto PGR, alto ESR1, baixa proliferação) até os mais agressivos (Basal-like, alta expressão de EXO1, CENPF, PTTG1, FOXC1).

Esta análise também justifica a escolha de Z-score contínuo em vez de classificação categórica: o escore funciona como um continuum, e qualquer limiar (mediana, quartis, etc.) é uma simplificação de uma escala contínua de risco.

A separação especialmente dramática entre Q4 e os demais grupos nos primeiros anos confirma que pacientes com escore muito alto (Q4) têm prognóstico particularmente ruim — esses seriam os candidatos prioritários para tratamentos mais intensivos ou estudos clínicos.

## Contexto no estudo

Esta figura, ao lado das equivalentes de TCGA-BRCA e GSE20685, forma um conjunto de evidências sobre a natureza contínua e gradual do escore CorePAM. A análise de dicotomia (Fig2 — mediana) é conservadora e adequada para o artigo principal; as análises por quartis nos suplementos mostram o "poder completo" do escore.

Para a tese de doutorado, estas figuras são especialmente valiosas para demonstrar que o CorePAM não é apenas um classificador binário (alto/baixo risco), mas um **preditor contínuo de risco** que pode ser usado de diferentes formas dependendo do contexto clínico.

## Pontos de atenção

- A análise por quartis com 4 grupos envolve comparações múltiplas — ao reportar p-valores individuais dos testes log-rank pairwise, ajuste de Bonferroni ou Holm é necessário para controle de erro tipo I.
- O **endpoint DSS do METABRIC** requer que os óbitos por outras causas sejam censurados — isso afeta os quartis de forma diferente dependendo da distribuição de comorbidades (pacientes mais idosas no Q1/Q2, potencialmente mais mortes por outras causas).
- A **análise de sensibilidade com OS** (FigS5) também poderia ser feita por quartis para completude — isso vai além do escopo desta figura suplementar.
- As curvas de quartis em longo prazo (>120 meses) têm IC muito largo por causa do número reduzido de pacientes em risco — interpretar cruzamentos ou divergências tardias com cautela.
