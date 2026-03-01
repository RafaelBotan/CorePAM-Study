# Interpretação: Fig3 — Curva de Kaplan-Meier, SCAN-B, Sobrevida Global (OS) — Coorte de Treinamento

## O que é este gráfico?

Este gráfico de Kaplan-Meier mostra a sobrevida global (OS) na coorte **SCAN-B** (Sweden Cancerome Analysis Network — Breast), que é a **coorte de treinamento** do estudo. As 3.069 pacientes foram estratificadas em dois grupos pelo escore CorePAM: **alto risco** (High) e **baixo risco** (Low), usando a mediana como ponto de corte.

É importante destacar o papel especial desta figura: ela mostra o desempenho do CorePAM nos **próprios dados de onde foi derivado**. Isso não é uma validação — é a demonstração de que o modelo aprendeu algo biologicamente real nos dados de treinamento. A validação genuína está nas Fig2 (TCGA-BRCA, METABRIC, GSE20685). No entanto, a Fig3 é igualmente necessária para contextualizar o desempenho nas validações: sabemos que o modelo nunca vai performar melhor que nos dados de treinamento, portanto o HR de SCAN-B serve como "teto teórico" para as validações.

## Como ler este gráfico

- **Eixo X (horizontal):** Tempo em meses desde o início do seguimento. SCAN-B tem seguimento razoável (mediana ~55 meses), mas não tão longo quanto METABRIC ou GSE20685.
- **Eixo Y (vertical):** Probabilidade de sobrevivência global acumulada.
- **Linha superior (azul — "Low risk"):** Pacientes com escore < mediana; perfil luminal, prognóstico favorável. A curva permanece elevada (próxima de 1,0) por muito mais tempo.
- **Linha inferior (vermelha — "High risk"):** Pacientes com escore ≥ mediana; perfil basal-like ou agressivo, maior mortalidade.
- **Separação precoce e sustentada:** A diferença entre as curvas começa cedo e se alarga progressivamente — um padrão típico de marcadores que capturam risco de forma contínua.
- **Tabela de risco:** Começa com ~1.534 pacientes em cada grupo (metade de 3.069 por dicotomização pela mediana) e vai diminuindo com eventos e censuras.
- **Área sombreada:** IC 95% pontual de cada curva — estreito aqui por causa do grande N.

## Achados principais

- **N total:** 3.069 pacientes (maior coorte do estudo)
- **Eventos (óbitos):** 322 (10,5%)
- **Tempo de seguimento mediano (censurado):** 55,7 meses (aproximadamente 4,6 anos)
- **Hazard Ratio univariado:** HR = 1,92 (IC 95%: 1,74–2,13)
- **p-valor:** p = 3,56 × 10⁻³⁶ (extraordinariamente significativo)
- **C-index:** 0,698 (IC 95%: 0,668–0,729) — o maior C-index do estudo
- **HR multivariado CORE-A (CorePAM + idade + ER):** HR = 1,86 (IC 95%: 1,64–2,12)
- **Delta C-index (CORE-A → CORE-A + CorePAM):** ΔC = +0,030 (IC 95%: 0,014–0,049)

## O que isso significa?

Os resultados do SCAN-B são os mais impressionantes numericamente do estudo. Um HR de 1,92 significa que, a cada ponto no tempo, pacientes de alto risco têm **92% mais risco instantâneo de morte** em comparação com os de baixo risco. O p-valor de 3,56 × 10⁻³⁶ corresponde a uma evidência estatística de magnitude astronômica.

O C-index de 0,698 é o melhor do estudo — próximo de 0,70 — indicando que o modelo distingue corretamente o prognóstico em aproximadamente 70% dos pares de pacientes comparáveis. Para referência, um lançamento de moeda seria 0,50; um modelo perfeito seria 1,0. C-index de 0,70 representa excelente discriminação para biomarcadores prognósticos em câncer de mama.

No entanto, o ponto mais fascinante é o ΔC relativamente **modesto** (+0,030) no SCAN-B — aparentemente contraditório com o HR tão alto. A explicação é importante: **o modelo CORE-A já é muito bom em SCAN-B** (C_COREA = 0,750), porque SCAN-B tem informações clínicas ricas (ER bem documentado, idade bem distribuída). Adicionar o CorePAM a um modelo clínico já forte produz ganho incremental menor. Já em METABRIC (onde CORE-A tem C = 0,494 — pouco melhor que chance), o CorePAM agrega +0,142 pontos de C-index.

O HR de 1,86 no modelo CORE-A (vs. 1,92 univariado) sugere que o escore CorePAM é **minimamente confundido** pelas variáveis clínicas — sua associação com sobrevida é quase completamente independente de idade e ER.

## Contexto no estudo

O SCAN-B serve como âncora metodológica do estudo. Todo o processo de derivação do CorePAM — seleção de genes, estimativa de pesos pelo elastic-net, validação cruzada de 10 folds — aconteceu exclusivamente em SCAN-B. A Fig3 mostra que o modelo resultante não é apenas estatisticamente significativo nos dados de treinamento, mas tem separação biológica real: as curvas divergem precocemente e a divergência se mantém ao longo de quase 6 anos de seguimento.

A decisão de usar **todas as 3.069 amostras do SCAN-B como treinamento** (sem split interno) é metodologicamente justificada pela validação cruzada de 10 folds usada na derivação: o C-index reportado no treinamento vem do processo OOF (Out-of-Fold), não de avaliação nos próprios dados de ajuste. Isso elimina o overfitting na seleção do modelo, embora a Fig3 em si use o modelo full-data aplicado a toda a coorte.

## Pontos de atenção

- **Esta figura não é uma validação independente** — é o desempenho nos dados de treinamento. O desempenho real e imparcial vem das curvas KM das coortes de validação (Fig2).
- O **seguimento mediano de 55,7 meses é relativamente curto** para câncer de mama, onde eventos tardios (>10 anos) são comuns, especialmente em Luminal A. Isso pode subestimar a diferença real entre os grupos.
- O **ΔC de +0,030** no SCAN-B, embora pequeno, é estatisticamente significativo (IC 95%: 0,014–0,049 — não inclui zero). Isso demonstra que o CorePAM agrega informação além dos fatores clínicos mesmo na coorte de treinamento.
- O HR de 1,92 é certamente **otimista por inflation de treinamento** — o verdadeiro efeito de generalizacao é melhor estimado pelos HRs das validações (1,20–1,41).
- A proporção de eventos (10,5%) é baixa para o tamanho da coorte, o que é típico de SCAN-B: coorte contemporânea, com tratamento moderno, principalmente ER+. Eventos são relativamente raros a 5 anos.
