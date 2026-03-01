# Interpretação: FigS — C-index por Coorte: Comparação de Modelos

## O que é este gráfico?

Esta figura suplementar apresenta o **C-index (índice de concordância de Harrell)** comparando diferentes modelos em cada coorte, tipicamente em formato de gráfico de barras ou pontos com IC 95%. Os modelos comparados são:

1. **CORE-A sozinho** (apenas variáveis clínicas: idade + ER)
2. **CorePAM score sozinho** (apenas o escore molecular)
3. **CORE-A + CorePAM** (modelo completo)

Eventualmente também inclui:
4. **PAM50 completo** (50 genes — para referência)
5. **Modelo nulo** (apenas intercepto, sem preditores — C-index = 0,50)

O objetivo é mostrar side-by-side o poder discriminativo de cada abordagem, tornando visual o ganho incremental do CorePAM.

## Como ler este gráfico

- **Cada grupo de barras:** Uma coorte (SCAN-B, TCGA-BRCA, METABRIC, GSE20685).
- **Barras dentro de cada grupo:** Diferentes modelos (CORE-A, CorePAM, CORE-A+CorePAM).
- **Eixo Y:** C-index (0,50 = nenhuma discriminação; 1,0 = discriminação perfeita).
- **Barras de erro:** IC 95% (bootstrap) de cada C-index.
- **Linha de referência em 0,50:** O "piso" de discriminação aleatória.
- **Comparação de sobreposição dos IC:** Se os IC de dois modelos não se sobrepõem, a diferença é estatisticamente significativa.

## Achados principais

| Coorte | C (CORE-A) | C (CorePAM) | C (CORE-A + CorePAM) |
|---|---|---|---|
| SCAN-B | 0,750 | 0,698 | 0,780 |
| TCGA-BRCA | 0,638 | 0,624 | 0,676 |
| METABRIC | 0,494 | 0,638 | 0,636 |
| GSE20685 | 0,536 | 0,623 | 0,628 |

Padrões notáveis:
- **SCAN-B:** CORE-A > CorePAM sozinho, mas CORE-A + CorePAM é o maior. As variáveis clínicas são muito informativas aqui.
- **METABRIC:** CORE-A < CorePAM sozinho! O escore molecular supera as variáveis clínicas. O CORE-A com C=0,494 é praticamente aleatório — a informação clínica disponível no METABRIC (para esta análise específica) é quase sem valor discriminativo. O CorePAM "resgata" a discriminação.
- **GSE20685:** CORE-A < CorePAM sozinho. Similar ao METABRIC — variáveis clínicas têm poder discriminativo limitado, e o escore molecular supera.
- **TCGA-BRCA:** CorePAM e CORE-A têm C-index próximos (0,624 vs. 0,638), e a combinação melhora para 0,676.

## O que isso significa?

Esta figura conta uma história fascinante sobre **quando os marcadores moleculares são mais valiosos**:

Em coortes onde as variáveis clínicas já são muito informativas (SCAN-B — coorte moderna, bem anotada, com ER e idade confiáveis), o ganho marginal do escore molecular é menor. O clínico experiente com ER e idade já consegue discriminar bem os pacientes de SCAN-B.

Em coortes onde as variáveis clínicas são menos informativas (METABRIC, GSE20685 — coortes históricas, possivelmente com dados clínicos menos padronizados, ou onde a heterogeneidade de subtipos não é capturada pelo ER clínico), o escore molecular **supera** completamente as variáveis clínicas.

Isso tem implicação clínica direta: o CorePAM seria **mais valioso** em contextos onde a informação clínica é limitada ou menos confiável — que é exatamente o cenário em muitos serviços de oncologia do mundo real, especialmente em países de baixa e média renda.

A mensagem do METABRIC (C_COREA = 0,494) é especialmente poderosa: sem o CorePAM, as variáveis clínicas disponíveis não conseguem discriminar prognóstico melhor do que chance. **Com o CorePAM, C sobe para 0,636** — transformando uma análise que seria estatisticamente inútil em uma com discriminação razoável.

## Contexto no estudo

Esta figura suplementar é o complemento visual da tabela de ΔC-index (Fig5 principal). Enquanto o ΔC-index mostra o **ganho**, esta figura mostra os **valores absolutos** — permitindo que o leitor julgue tanto o ponto de partida (CORE-A) quanto o destino (modelo completo).

A figura também indiretamente justifica o modelo CORE-A como baseline: se tivéssemos usado um modelo clínico mais complexo (com grau, tamanho, linfonodos), o baseline seria mais alto e o ΔC aparente menor. A escolha de um modelo clínico simples (apenas idade + ER) é conservadora e generalizável — mas o trade-off é um baseline potencialmente mais baixo que o obtido na prática clínica com mais informações.

## Pontos de atenção

- O **C-index de CORE-A** varia substancialmente entre coortes (0,494 a 0,750). Isso reflete **diferenças na qualidade dos dados clínicos** e na relevância das variáveis clínicas em cada população — não apenas diferenças biológicas.
- A **interpretação do C-index** depende do endpoint: C-index para DSS (METABRIC) não é diretamente comparável ao C-index para OS (outras coortes).
- **C-index = 0,636 no METABRIC** para o modelo CorePAM sozinho pode parecer modesto, mas é notável considerando que é obtido sem nenhuma variável clínica e com pesos estimados em outra plataforma (RNA-seq).
- Os IC 95% dos C-indexes individuais são estimados por bootstrap — a comparação formal entre modelos dentro da mesma coorte usa o **ΔC bootstrapado** (apresentado na Fig5 principal), que é mais preciso para comparações pareadas.
- Um C-index de 0,70 em oncologia mamária é geralmente considerado bom — o padrão de referência de modelos clínicos como Oncotype DX RS tem C-index em torno de 0,60-0,65 em populações comparáveis.
