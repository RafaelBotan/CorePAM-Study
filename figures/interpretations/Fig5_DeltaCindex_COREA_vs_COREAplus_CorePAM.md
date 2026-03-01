# Interpretação: Fig5 — Delta C-index: CORE-A vs. CORE-A + CorePAM

## O que é este gráfico?

Este gráfico mostra o **ganho incremental de capacidade discriminativa** (Delta C-index, ΔC) obtido ao adicionar o escore CorePAM a um modelo clínico de base — o modelo CORE-A, que inclui apenas **idade** e **status de ER** (receptor de estrogênio).

O C-index (índice de concordância de Harrell) é a métrica padrão de discriminação em análise de sobrevida. Varia de 0,50 (sem discriminação, equivalente a jogar moeda) a 1,0 (discriminação perfeita). O **ΔC = C(CORE-A + CorePAM) − C(CORE-A)** mede especificamente quanto o CorePAM adiciona além do que os clínicos já sabem.

Este gráfico é a resposta à pergunta clínica mais importante do estudo: **"O CorePAM agrega valor real a um modelo clínico simples?"**

## Como ler este gráfico

O gráfico tipicamente aparece como um **gráfico de barras com barras de erro** (IC 95%), com uma barra por coorte. Cada barra mostra:

- **Altura da barra:** O valor de ΔC (ganho em C-index).
- **Barras de erro:** IC 95% do ΔC (estimado por bootstrap).
- **Linha de referência em ΔC = 0:** Se o IC inteiramente acima de zero, o ganho é estatisticamente significativo.
- **Cores diferentes** representam cada coorte (SCAN-B, TCGA-BRCA, METABRIC, GSE20685).

Também pode ser apresentado como **tabela de C-index** comparando CORE-A vs. CORE-A + CorePAM side-by-side, com o ΔC calculado.

## Achados principais

| Coorte | C (CORE-A) | C (CORE-A + CorePAM) | ΔC | IC 95% do ΔC |
|---|---|---|---|---|
| SCAN-B (treino) | 0,750 | 0,780 | +0,030 | (0,014 – 0,049) |
| TCGA-BRCA (valid.) | 0,638 | 0,676 | +0,038 | (0,011 – 0,086) |
| METABRIC (valid.) | 0,494 | 0,636 | **+0,142** | (0,097 – 0,163) |
| GSE20685 (valid.) | 0,536 | 0,628 | +0,093 | (0,021 – 0,175) |

Todos os IC 95% estão **acima de zero**, confirmando ganho incremental significativo em todas as coortes.

## O que isso significa?

Os resultados do ΔC-index são os dados mais diretamente clínicos do estudo. Vamos interpretar cada coorte:

**SCAN-B (ΔC = +0,030):** O ganho é modesto — mas C_COREA já era 0,750! Isso significa que o modelo clínico em SCAN-B é muito informativo (age e ER têm boa distribuição e são bem coletados). Mesmo assim, adicionar CorePAM eleva o C-index em 3 pontos percentuais, com IC claramente positivo. Em oncologia, um ΔC de +0,030 sobre um modelo já bom é clinicamente relevante.

**TCGA-BRCA (ΔC = +0,038):** Ganho um pouco maior que em SCAN-B, mas o IC é mais largo (0,011–0,086) por causa do menor N e seguimento curto. O C_COREA de 0,638 é mais baixo que em SCAN-B, indicando que as variáveis clínicas disponíveis no TCGA são menos informativas para esta coorte. O CorePAM eleva para 0,676 — um salto considerável.

**METABRIC (ΔC = +0,142):** Este é o resultado mais espetacular. O modelo CORE-A no METABRIC tem C = 0,494 — **praticamente equivalente a jogar moeda!** Isso ocorre porque o METABRIC tem variáveis clínicas limitadas ou codificadas de forma diferente (ER medido por imunohistoquímica em era diferente, com cortes diferentes). Adicionar o escore CorePAM eleva o C-index para 0,636 — um ganho de 14,2 pontos percentuais. Isso é extraordinário e demonstra que o CorePAM "resgatou" toda a capacidade discriminativa que as variáveis clínicas isoladas não conseguiam fornecer.

**GSE20685 (ΔC = +0,093):** Ganho substancial de ~9 pontos percentuais, a partir de um C_COREA de 0,536 (também baixo — variáveis clínicas pobres nesta coorte). O IC amplo (0,021–0,175) reflete a incerteza inerente ao N=327, mas o limite inferior (0,021) ainda está bem acima de zero.

**O padrão geral é revelador:** O CorePAM agrega mais valor onde as variáveis clínicas são menos informativas (METABRIC, GSE20685) e menos onde elas já são muito informativas (SCAN-B). Isso é exatamente o comportamento esperado de um bom biomarcador molecular: **complementa, não substitui** a informação clínica.

## Contexto no estudo

A Fig5 (ΔC-index) junto com a Fig5 de calibração formam o argumento de **valor clínico incremental** do CorePAM. Enquanto as curvas KM (Fig2-3) e o forest plot (Fig4) mostram que há associação prognóstica, a Fig5 quantifica **o quanto** essa associação é útil além do que o clínico já sabe.

Este tipo de análise — análise incremental — é cada vez mais exigida por revisores e agências de saúde para biomarcadores moleculares. Não basta mostrar que um marcador "funciona"; é preciso demonstrar que ele funciona **melhor** (ou complementa) o que já existe. O CorePAM passa nesse teste em todas as 4 coortes.

## Pontos de atenção

- O **modelo CORE-A é deliberadamente simples** (apenas idade e ER) para maximizar a generalização e a aplicabilidade prática. Modelos mais complexos com mais covariáveis clínicas (grau tumoral, tamanho, linfonodos) poderiam elevar o C_COREA e, consequentemente, reduzir o ΔC observável.
- Os IC 95% são estimados por **bootstrap**, que é um método não-paramétrico robusto, mas a precisão depende do número de replicações de bootstrap e do tamanho amostral.
- O ΔC de +0,038 no TCGA-BRCA tem IC inferior de 0,011 — muito próximo de zero — indicando que o ganho real pode ser pequeno nesta coorte especificamente (influenciado pelo seguimento curto de 32 meses).
- O METABRIC tem o maior ΔC mas também o maior IC inferior (0,097) — ou seja, o ganho mínimo plausível ainda é de ~10 pontos percentuais, o que é uma garantia estatística forte.
- A **Decision Curve Analysis (DCA)**, que avalia benefício clínico líquido em diferentes limiares de risco, complementa esta análise e é apresentada na figura de calibração.
