# Interpretação: FigS — Forest Plot: HR nas Coortes de Validação (Análises Secundárias)

## O que é este gráfico?

Este forest plot suplementar apresenta os Hazard Ratios (HRs) do escore CorePAM em **análises secundárias e subgrupos** nas coortes de validação. Complementa o forest plot principal (Fig4), que mostra os resultados primários por coorte.

Dependendo do conteúdo específico, este forest plot suplementar pode incluir:

1. **HRs por subtipo molecular** (Luminal A, Luminal B, HER2-enriched, Basal-like) — mostrando se o CorePAM funciona dentro de cada subtipo separadamente.
2. **HRs por faixa etária** (<50 anos vs. ≥50 anos) — avaliando se o efeito difere por idade.
3. **HRs por status de ER** (ER+ vs. ER−) — testando se o escore agrega valor além do ER dentro de cada grupo.
4. **HRs por coorte individual** com modelo univariado E multivariado lado a lado.
5. **Análises de interação** testando se o efeito do CorePAM é modificado por covariáveis clínicas.

O objetivo é demonstrar **consistência do efeito** e identificar possíveis **subgrupos com maior ou menor benefício** do escore.

## Como ler este gráfico

- **Estrutura similar ao forest plot principal (Fig4):** Linhas horizontais com quadrados (estimativas pontuais) e barras de erro (IC 95%).
- **Linha de nulidade em HR=1:** Quadrados à direita = maior risco com escore alto.
- **Subgrupos aninhados:** Se o gráfico mostrar subgrupos dentro de cada coorte, há linhas hierárquicas com níveis diferentes.
- **p de interação:** Testa se o HR em subgrupo A é significativamente diferente do HR em subgrupo B. Um p de interação alto (>0,10) indica que não há evidência de modificação do efeito.

## Achados principais

Para as coortes de validação (TCGA-BRCA, METABRIC, GSE20685), os resultados dos subgrupos esperados:

**Por status de ER:**
- **ER+:** O escore CorePAM deve ainda mostrar HR > 1 dentro dos ER+, confirmando que captura variação prognóstica dentro do grupo luminal (não apenas discrimina ER+ de ER−).
- **ER−:** Similarmente, dentro dos ER−, deve haver gradiente de risco.
- **p de interação:** Tipicamente não significativo — o efeito do CorePAM é similar em ER+ e ER−.

**Por subtipo PAM50 (se disponível nas coortes):**
- **Luminal A vs. B:** O CorePAM pode ter maior poder discriminativo em Luminal B (que é mais heterogêneo) do que em Luminal A (mais uniforme).
- **Basal-like:** Com menos variação intra-subtipo, o escore pode ter menor poder discriminativo.

**Comparação univariado vs. multivariado (CORE-A):**
- Os HRs ligeiramente reduzidos no modelo multivariado (vs. univariado) indicam que ER e idade confundem parcialmente o efeito — o que é esperado. O IC continuando acima de 1,0 no multivariado confirma independência prognóstica.

## O que isso significa?

O forest plot de subgrupos aborda uma questão crítica: **"o CorePAM funciona apenas para discriminar luminal de basal, ou tem valor dentro de cada subtipo?"**

Se o HR fosse significativo apenas para a comparação ER+ vs. ER− (ou Luminal vs. Basal), isso sugeriria que o CorePAM é apenas um "proxy" da informação já disponível pelo ER. Porém, se o HR permanece significativo **dentro** dos subgrupos ER+ e ER−, demonstra que o escore captura heterogeneidade prognóstica além da simples classificação ER.

Isso é particularmente relevante para **Luminal B**: pacientes ER+ de alto risco que se beneficiariam de quimioterapia adicional (além da hormonoterapia) são difíceis de identificar clinicamente. O CorePAM, ao estratificar risco dentro do grupo ER+, pode ajudar exatamente nessa decisão.

A consistência do efeito entre subgrupos (sem interação significativa) é também um critério de qualidade: indica que o escore não funciona apenas em condições específicas, mas genericamente em toda a população estudada.

## Contexto no estudo

O forest plot de subgrupos (FigS) é complementar ao forest plot principal (Fig4). Enquanto Fig4 mostra resultados por coorte (síntese cross-coorte), FigS mostra resultados por subgrupo (síntese cross-biologia).

Para a narrativa do artigo, estas figuras juntas constroem o argumento: "o CorePAM tem efeito prognóstico robusto, consistente entre coortes, plataformas, populações e subgrupos moleculares." Essa consistência multinível é o padrão-ouro para validação de biomarcadores prognósticos.

## Pontos de atenção

- **Análises de subgrupos têm poder reduzido**: Ao dividir cada coorte em subgrupos (ex: ER+ e ER−), o número de eventos por subgrupo cai, reduzindo o poder estatístico. Um HR não significativo em um subgrupo pequeno não é evidência de ausência de efeito.
- **Multiplicidade de testes**: Com muitos subgrupos testados, o risco de falso-positivo aumenta. As análises de subgrupos pré-especificadas (como ER+ vs. ER−) têm mais validade que subgrupos exploratórios.
- **p de interação é geralmente conservador**: Mesmo quando o HR parece diferente entre subgrupos visualmente, o p de interação raramente atinge significância estatística com N moderado — isso não significa que não há diferença real.
- **Definição de subtipos varia por coorte**: Subtipo PAM50 pode ser determinado por IHC clínica ou por expressão gênica (molecular). Usar definições inconsistentes entre coortes compromete comparações de subgrupos.
- **O forest plot suplementar é exploratório** — as conclusões principais do estudo baseiam-se na análise por coorte (Fig4), não nos subgrupos. Este gráfico gera hipóteses para estudos futuros.
