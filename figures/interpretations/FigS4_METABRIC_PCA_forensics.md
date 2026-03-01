# Interpretação: FigS4 — PCA Forense do METABRIC

## O que é este gráfico?

Este gráfico suplementar apresenta uma **Análise de Componentes Principais (PCA)** forense dos dados de expressão do METABRIC, com foco nos genes do CorePAM. O adjetivo "forense" é proposital: esta análise foi conduzida para **investigar** e **compreender** a estrutura dos dados METABRIC, identificar potenciais subgrupos, outliers, ou fontes de variação técnica que poderiam confundir a análise de sobrevida.

A PCA é uma técnica de redução de dimensionalidade que transforma um conjunto de variáveis correlacionadas (expressões de múltiplos genes) em um número menor de variáveis não correlacionadas (componentes principais, PCs), preservando o máximo de variância dos dados originais. Os primeiros PCs capturam as principais fontes de variação.

## Como ler este gráfico

- **Cada ponto:** Representa uma amostra (paciente) do METABRIC.
- **Eixo X (PC1):** Primeiro componente principal — captura a maior proporção de variância. Frequentemente corresponde ao "gradiente luminal ↔ basal" em dados de câncer de mama.
- **Eixo Y (PC2):** Segundo componente principal — segunda maior fonte de variação ortogonal ao PC1. Pode corresponder a gradiente HER2+ vs. outros, ou a efeitos de batch técnico.
- **Cores dos pontos:** Codificam informações biológicas ou técnicas — tipicamente:
  - Subtipo molecular PAM50 (Luminal A, Luminal B, HER2-enriched, Basal-like, Normal-like)
  - Status de ER
  - Quartil do escore CorePAM
  - Batch de microarray (se disponível)
- **Percentual de variância explicada (em cada eixo):** Informa quanta variação cada PC captura — ex: "PC1: 42.3%" significa que o primeiro componente captura 42,3% de toda a variância dos dados.
- **Clusters visuais:** Grupos de pontos agrupados refletem pacientes biologicamente similares. Separação clara entre subtipos valida a qualidade dos dados.
- **Outliers:** Pontos muito distantes do cluster central merecem investigação (amostras potencialmente corrompidas, trocadas, ou biologicamente extremas).

## Achados principais

A PCA forense do METABRIC tipicamente revela:

1. **Estrutura clara de subtipos:** No espaço dos genes CorePAM, os subtipos PAM50 se separam de forma razoavelmente clara — os Basais formam um cluster separado, enquanto os Luminais se agrupam do outro lado do PC1.

2. **Gradiente luminal-basal no PC1:** O PC1 geralmente captura o eixo principal de biologia de câncer de mama — alta expressão de genes proliferativos/basais (EXO1, PTTG1, FOXC1, KRT5, CENPF) vs. alta expressão de genes luminais (ESR1, PGR, BCL2, MLPH, NAT1).

3. **Ausência de batch effects óbvios:** Se os pontos fossem coloridos por lote de microarray e formassem clusters separados, indicaria confundimento técnico. A ausência dessa estrutura é desejável.

4. **Verificação de outliers:** A análise identificou se há amostras fora do espaço normal dos dados METABRIC que deveriam ser excluídas ou investigadas.

5. **Concordância do escore CorePAM com PC1:** Espera-se que o quartil do escore CorePAM se alinha bem com a posição no PC1 — confirmando que o escore captura o principal eixo de variação molecular.

## O que isso significa?

A PCA forense do METABRIC serve múltiplos propósitos:

**Validação de qualidade dos dados:** Dados de microarray Illumina de uma coorte coletada ao longo de décadas podem ter heterogeneidade técnica (diferentes lotes de hibridização, diferentes operadores). A PCA permite visualizar se há estrutura técnica confundindo a biológica.

**Confirmação da biologia:** A separação dos subtipos moleculares no espaço PCA dos genes CorePAM confirma que o painel de 24 genes captura os principais eixos de variação biológica do câncer de mama — o mesmo sinal que os 50 genes do PAM50 capturavam com maior redundância.

**Diagnóstico de outliers:** A identificação e tratamento de outliers é particularmente importante no METABRIC, que inclui amostras de várias décadas e centros. Amostras com perfil completamente atípico (ex: possível contaminação ou troca de amostra) seriam excluídas antes das análises de sobrevida.

**Justificativa metodológica:** A análise forense documenta o raciocínio por trás de decisões de processamento de dados — essencial para reprodutibilidade.

## Contexto no estudo

O METABRIC é a maior coorte de validação e a mais completa em termos de seguimento. Qualquer problema não detectado nos dados METABRIC teria grande impacto nos resultados gerais do estudo. A PCA forense (FigS4) é o "laudo pericial" que atesta a qualidade dos dados antes de eles entrarem nas análises de sobrevida.

Esta figura suplementar é especialmente valorizada em revisão por pares: mostra que os autores não apenas usaram os dados, mas os **entenderam** e **auditaram** antes de usá-los.

## Pontos de atenção

- A PCA com apenas 24 genes (CorePAM) captura menos variância total que uma PCA com todos os genes do transcriptoma — os PCs são mais específicos para o sinal do CorePAM.
- **Separação visual na PCA não implica separação estatística** nas análises de sobrevida — são análises complementares, não equivalentes.
- A interpretação dos PCs como "gradiente luminal-basal" é inferida pela associação com marcadores conhecidos (ESR1, KRT5, etc.) — não é automática nem garantida.
- Amostras identificadas como outliers na PCA devem ser investigadas antes de serem excluídas — a exclusão arbitrária pode introduzir viés de seleção.
- A PCA é sensível à **escala** dos dados — dados não normalizados podem ter PC1 dominado por genes de alta expressão absoluta em vez de variância biológica. A normalização com Z-score intra-coorte é crucial antes desta análise.
