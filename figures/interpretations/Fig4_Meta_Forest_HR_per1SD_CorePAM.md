# Interpretação: Fig4 — Forest Plot da Meta-análise: HR por 1 DP do Escore CorePAM

## O que é este gráfico?

Este é um **forest plot** (gráfico de floresta), o gráfico canônico de meta-análise. Ele mostra o efeito do escore CorePAM sobre sobrevida em cada coorte individualmente, e depois combina esses efeitos em uma estimativa global pela meta-análise de efeitos aleatórios (random-effects meta-analysis).

O que está sendo medido aqui é o **HR por 1 desvio padrão (DP) do escore CorePAM** — ou seja, o risco relativo de morte para cada aumento de 1 DP no escore. Usar 1 DP (em vez de dicotomia alta/baixa) é mais informativo porque trata o escore como contínuo, evitando perda de informação pelo corte na mediana.

Este gráfico é possivelmente a figura mais importante do manuscrito do ponto de vista de síntese de evidência. Ele responde à pergunta central: **"Considerando todas as coortes juntas, o CorePAM tem associação prognóstica robusta e replicável?"**

## Como ler este gráfico

- **Cada linha horizontal** representa uma coorte. O quadrado no centro da linha é a estimativa pontual do HR; a linha horizontal é o IC 95%.
- **Tamanho do quadrado:** Proporcional ao peso da coorte na meta-análise (coortes maiores/mais precisas têm peso maior).
- **Posição vertical:** As coortes são listadas individualmente acima, com o sumário da meta-análise (diamante) na parte inferior.
- **Linha vertical em HR=1,0:** É a linha de nulidade. Se o quadrado e todo o IC estão à **direita** de 1,0, o escore está associado a **maior risco** (HR > 1).
- **Diamante na parte inferior:** Representa a estimativa combinada da meta-análise. A largura do diamante = IC 95% da estimativa combinada.
- **I² (heterogeneidade):** Mede a proporção de variação entre estudos que é devida a diferenças reais (e não ao acaso). I² alto (>50%) indica heterogeneidade substancial.
- **τ² (tau-quadrado):** Variância entre estudos nos modelos de efeitos aleatórios. Mede o quanto os verdadeiros efeitos variam entre coortes.

Duas versões deste forest plot existem:
1. **Univariada:** HR do escore CorePAM sozinho.
2. **CORE-A multivariada:** HR do escore CorePAM ajustado por idade e ER (covariáveis clínicas).

## Achados principais

**Meta-análise univariada (4 coortes):**
- SCAN-B: HR = 1,92 (IC 95%: 1,74–2,13)
- TCGA-BRCA: HR = 1,20 (IC 95%: 1,04–1,40)
- METABRIC: HR = 1,41 (IC 95%: 1,31–1,52)
- GSE20685: HR = 1,40 (IC 95%: 1,12–1,76)
- **Meta HR = 1,472 (IC 95%: 1,205–1,798)**
- **p-valor meta = 1,52 × 10⁻⁴**
- **I² = 90,6%** (alta heterogeneidade)
- **τ² = 0,036**

**Meta-análise CORE-A multivariada (3 coortes de validação):**
- TCGA-BRCA: HR_adj = 1,27 (IC 95%: 1,09–1,48)
- METABRIC: HR_adj = 1,41 (IC 95%: 1,31–1,52)
- GSE20685: HR_adj = ? (incluso na meta)
- **Meta HR_adj = 1,478 (IC 95%: 1,253–1,744)**
- **p-valor meta = 3,69 × 10⁻⁶**
- **I² = 83,6%**

## O que isso significa?

O resultado da meta-análise é cristalino: **em média, um aumento de 1 DP no escore CorePAM está associado a um aumento de 47,2% no risco de óbito** (Meta HR = 1,472), com evidência altamente significativa (p = 1,52 × 10⁻⁴). Esse efeito é consistente — todos os IC 95% individuais estão inteiramente à direita de 1,0 exceto o TCGA-BRCA, que tangencia levemente (1,04–1,40).

A **heterogeneidade alta** (I² = 90,6%) é esperada e clinicamente justificável por várias razões:
1. **Plataformas diferentes:** RNA-seq vs. microarray Illumina vs. Affymetrix — escala de sinal e ruído diferem.
2. **Populações diferentes:** Sueca (SCAN-B), americana (TCGA), europeia (METABRIC), sul-coreana (GSE20685).
3. **Endpoints diferentes:** OS em três coortes, DSS em METABRIC.
4. **Duração de seguimento muito diferente:** De 32 meses (TCGA) a 159 meses (METABRIC).
5. **Era de tratamento diferente:** SCAN-B e TCGA são modernos; METABRIC e GSE20685 incluem casos mais antigos.

Por isso, o modelo de **efeitos aleatórios** é o correto aqui — ele não assume que todas as coortes estimam o mesmo efeito verdadeiro, mas sim que cada coorte tem seu próprio HR verdadeiro, e esses HRs seguem uma distribuição. O τ² = 0,036 quantifica essa dispersão (relativamente modesta na escala log-HR).

O resultado multivariado (CORE-A) com Meta HR = 1,478 e I² = 83,6% é ainda mais importante: mostra que o escore CorePAM é **prognosticamente independente** das variáveis clínicas de uso rotineiro (idade e status de ER). Isso é o requisito mínimo para que um biomarcador molecular seja clinicamente útil.

## Contexto no estudo

O forest plot da meta-análise (Fig4) é a síntese definitiva do estudo. Após mostrar separação das curvas KM coorte a coorte (Fig2 e Fig3), este gráfico combina formalmente as evidências e produz uma estimativa de efeito global com toda a força estatística das 4 coortes somadas (N total > 6.400 pacientes, >1.200 eventos).

A meta-análise com I² alto confirma que os efeitos são **heterogêneos mas consistentemente positivos** — ou seja, o CorePAM tem associação prognóstica em todas as coortes, embora com magnitudes que variam por razões biológicas e tecnológicas legítimas. Isso é superior a um I² = 0% que poderia sugerir artificialmente que todos os estudos são idênticos.

## Pontos de atenção

- **I² = 90,6% é alto** — isso significa que a maioria da variação observada entre coortes é real (não acaso). A interpretação correta é: "há heterogeneidade substancial, portanto o efeito médio deve ser interpretado com cautela como uma estimativa pooled de efeitos variáveis."
- O **SCAN-B domina o peso** na meta-análise por ter o IC 95% mais estreito (N=3069). Isso pode puxar a estimativa combinada em direção ao seu HR mais alto (1,92).
- A **ausência de SCAN-B na meta-análise multivariada CORE-A** (já que não é coorte de validação) faz sentido metodológico: SCAN-B é onde os pesos foram estimados, incluí-lo na meta de validação seria conservadorismo incorreto.
- O número relativamente pequeno de estudos (4) limita a análise de viés de publicação (funnel plot). Com apenas 4 estudos, o poder para detectar assimetria é baixo.
