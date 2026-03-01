# Interpretação: FigS3 — Correlações Off-Diagonal entre Coortes

## O que é este gráfico?

Este gráfico suplementar apresenta a análise de **correlação de Pearson off-diagonal** dos escores CorePAM entre coortes, verificando a **ausência de contaminação de dados** (data leakage) entre a coorte de treinamento (SCAN-B) e as coortes de validação (TCGA-BRCA, METABRIC, GSE20685).

Em termos práticos, o gráfico mostra matrizes de correlação (ou heatmaps de correlação) dos perfis de expressão dos genes CorePAM entre pares de coortes, com foco especial nas correlações **fora da diagonal** — ou seja, correlações entre amostras de coortes diferentes.

A regra metodológica central do estudo é: **zero pooling entre coortes, zero batch correction inter-coorte**. Esta figura é a evidência visual de que essa regra foi respeitada e que não há dependência espúria entre os conjuntos de dados.

## Como ler este gráfico

- **Matriz de correlação:** Cada célula (i,j) representa a correlação de Pearson entre o escore (ou perfil de expressão) da amostra i e da amostra j.
- **Diagonal:** Representa correlações de cada amostra consigo mesma — sempre 1,0 por definição.
- **Off-diagonal intra-coorte:** Correlações entre amostras da mesma coorte — esperadas serem moderadas a altas por compartilharem plataforma e origem.
- **Off-diagonal inter-coorte (o que importa):** Correlações entre amostras de coortes diferentes — devem ser baixas e sem estrutura sistemática se não há leakage.
- **Escala de cores:** Tipicamente vermelho = alta correlação, azul = baixa correlação (ou vice-versa).
- **Blocos visíveis:** Estrutura em blocos ao longo da diagonal indica separação clara entre coortes (desejável). Ausência de blocos fora da diagonal confirma independência.

## Achados principais

A análise off-diagonal está desenhada especificamente para detectar três tipos de problemas:

1. **Contaminação direta:** Amostras da coorte de validação que, acidentalmente, foram incluídas no treinamento (ou vice-versa) — se isso ocorresse, veríamos correlações anormalmente altas entre pares específicos inter-coorte.

2. **Batch effects compartilhados:** Se as mesmas amostras foram processadas em conjunto ou se há sobreposição de lotes, correlações inter-coorte seriam sistematicamente infladas.

3. **Sobreposição de pacientes:** Se a mesma paciente aparece em duas coortes com identificadores diferentes (especialmente TCGA e METABRIC, que têm populações overlapping de alguns centros europeus).

O resultado esperado e desejado: **correlações off-diagonal inter-coorte próximas a zero, sem estrutura** — confirmando que cada coorte é estatisticamente independente das demais.

## O que isso significa?

A análise off-diagonal é um controle de qualidade fundamental de **integridade metodológica**. Em estudos com validação externa, o pecado capital é o vazamento de dados entre treinamento e validação — mesmo pequenos vazamentos inflam artificialmente as métricas de validação e tornam as estimativas de generalização otimistas e não-confiáveis.

No CorePAM, as salvaguardas contra leakage incluem:
- Z-score intra-coorte estritamente separado (cada coorte normaliza com suas próprias médias e desvios padrão)
- Nenhum batch correction inter-coorte
- Pesos dos genes estimados exclusivamente em SCAN-B
- Análise de sobrevida por coorte realizada de forma separada

Esta figura serve como **auditoria pública** dessas salvaguardas — qualquer revisor ou leitor pode verificar que não há estrutura anômala nas correlações inter-coorte.

## Contexto no estudo

O FigS3 aparece como figura suplementar não porque seja menos importante metodologicamente, mas porque é uma análise de controle de qualidade que não contribui diretamente para os resultados biológicos. Sua posição nos suplementos reflete a convenção de que figuras de QC (quality control) ficam separadas das figuras de resultado.

Para a defesa de tese, esta figura é especialmente valorizada por bancas que conhecem os problemas de contaminação de dados em estudos de expressão gênica — ela demonstra rigor metodológico explícito.

## Pontos de atenção

- A análise off-diagonal testa correlação de **Pearson** — um método paramétrico sensível a outliers e à distribuição dos dados. Se necessário, Spearman é mais robusto.
- Correlações próximas de zero entre coortes **não garantem** ausência de todo tipo de leakage — por exemplo, se o leakage for em forma de seleção de features (usar informação das validações para selecionar genes), isso não seria detectado por correlações de expressão.
- A análise confirma independência das **expressões dos genes CorePAM** entre coortes, mas não testa independência de variáveis clínicas (que poderiam ter sobreposição de centros, por exemplo).
- Uma análise PCA mostrando separação clara entre clusters de coortes (FigS4 para METABRIC) complementa esta análise de correlação.
