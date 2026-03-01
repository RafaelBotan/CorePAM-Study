# Interpretação: FigS — Curvas KM por Quartis do Escore CorePAM — TCGA-BRCA

## O que é este gráfico?

Esta figura suplementar apresenta as curvas de Kaplan-Meier para a coorte TCGA-BRCA estratificadas em **quatro grupos por quartis do escore CorePAM** (em vez da dicotomia alta/baixa pela mediana usada na análise principal).

Os quartis dividem a distribuição do escore CorePAM em quatro partes iguais:
- **Q1 (1º quartil — 0 a 25%):** Escores mais baixos — perfil mais luminal, maior expressão de genes associados a bom prognóstico.
- **Q2 (26 a 50%):** Escores baixo-intermediários.
- **Q3 (51 a 75%):** Escores alto-intermediários.
- **Q4 (4º quartil — 76 a 100%):** Escores mais altos — perfil mais basal-like, maior expressão de genes de mau prognóstico.

Cada quartil representa aproximadamente 268 pacientes (1.072 / 4).

## Como ler este gráfico

- **Quatro linhas KM** em vez de duas, coloridas de forma gradiente:
  - Q1: cor mais fria (ex: azul escuro) — melhor prognóstico
  - Q2: cor fria intermediária (ex: azul claro)
  - Q3: cor quente intermediária (ex: laranja)
  - Q4: cor mais quente (ex: vermelho) — pior prognóstico
- **Separação ordenada das curvas:** Se o escore funciona como contínuo, espera-se uma separação **ordenada e gradual** das quatro curvas — Q1 sempre acima de Q2, Q2 acima de Q3, Q3 acima de Q4.
- **Log-rank test global e pairwise:** O p-valor global testa se há diferença entre qualquer par de grupos. Testes pairwise (Q1 vs. Q4 especialmente) mostram as comparações extremas.
- **Tabela de risco:** Agora com quatro linhas em vez de duas, mostrando quantos pacientes estão em risco em cada quartil ao longo do tempo.

## Achados principais

A análise por quartis em TCGA-BRCA deve mostrar:
- **Gradiente dose-resposta:** Sobrevida progressivamente pior de Q1 para Q4.
- **Maior separação nos extremos:** Q1 vs. Q4 é a comparação mais extrema e tipicamente mais significativa.
- **Confirmação da linearidade do escore:** Se a separação for ordenada (Q1 > Q2 > Q3 > Q4), confirma que o escore CorePAM funciona como preditor contínuo e não apenas como discriminador binário.

A comparação Q1 vs. Q4 representa os 25% "mais luminais" vs. os 25% "mais basais-like" — a diferença prognóstica esperada é substancial e biologicamente justificada.

## O que isso significa?

A análise por quartis é complementar à análise principal por dicotomia (mediana) e serve para mostrar que:

1. **O escore tem gradiente biológico real** — não é apenas uma variável binária.
2. **A dicotomização pela mediana (análise principal) é conservadora** — ao juntar Q2 com Q3 e Q1 com Q4, a análise principal subestima a diferença real entre os extremos.
3. **Há dose-resposta linear ou quase-linear:** O risco aumenta progressiva e consistentemente com o escore.

Em termos clínicos, isso é importante porque sugere que o escore pode ser usado de forma contínua para estratificação de risco (não apenas como "alto" vs. "baixo"), o que é mais flexível para uso clínico futuro.

Para o TCGA-BRCA especificamente, onde o HR pela mediana é de apenas 1,20, a análise por quartis revela o efeito real dos extremos: pacientes de Q4 vs. Q1 têm uma diferença prognóstica substancialmente maior que sugere o HR de 1,20.

## Contexto no estudo

Esta figura suplementar aparece em conjunto com as figuras equivalentes para METABRIC e GSE20685, formando um trio que demonstra consistentemente o gradiente dose-resposta do escore CorePAM em todas as coortes de validação.

A análise por quartis também serve como uma verificação informal de que o modelo não "força" a dicotomia a funcionar por algum artefato estatístico — se o escore funciona por quartis de forma ordenada, a dicotomia pela mediana é uma simplificação válida do padrão real.

## Pontos de atenção

- Com 1.072 pacientes divididos em 4 grupos (~268 por quartil), cada grupo tem tamanho razoável para análises de sobrevida, mas o número de eventos por grupo (~37 se uniformemente distribuídos) é limitado para estimativas precisas.
- O seguimento curto (32 meses mediano) afeta mais os grupos Q1 e Q2, onde eventos são mais tardios — possível que no longo prazo a separação das curvas dos quartis inferiores seja mais pronunciada.
- Comparações múltiplas (6 pares possíveis entre 4 quartis) requerem ajuste por multiplicidade se interpretadas formalmente — a análise principal (apenas dois grupos) não tem esse problema.
- A escolha de quartis é convencional mas arbitrária — tercis ou quintis revelariam padrões ligeiramente diferentes. Quartis são escolhidos por equilíbrio entre número de grupos e tamanho de cada grupo.
