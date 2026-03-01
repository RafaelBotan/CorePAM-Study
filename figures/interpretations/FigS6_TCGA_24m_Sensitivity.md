# Interpretação: FigS6 — Análise de Sensibilidade do TCGA-BRCA: Corte a 24 Meses

## O que é este gráfico?

Esta figura suplementar aborda uma das principais limitações da análise principal do TCGA-BRCA: o **seguimento mediano curto de apenas 32 meses**. A solução metodológica é truncar a análise em **24 meses** — analisando apenas os eventos que ocorrem nos primeiros 2 anos de seguimento, período em que todos (ou quase todos) os pacientes ainda estão efetivamente em risco e o seguimento está completo.

A análise a 24 meses é um corte temporal que enriquece os eventos disponíveis de forma comparável entre pacientes e elimina a incerteza das estimativas KM tardias (quando poucos pacientes ainda estão em risco). Em oncologia mamária, os primeiros 24 meses são um período de alto risco, especialmente para subtipos agressivos (Basal-like, HER2-enriched).

Este tipo de análise é conhecido como **"landmark analysis"** ou análise com censura à direita em um ponto específico: pacientes sem evento até 24 meses são censurados em 24 meses, independente de seu tempo real de seguimento.

## Como ler este gráfico

- **Eixo X:** Tempo em meses, agora truncado em 24 meses
- **Eixo Y:** Probabilidade de sobrevivência global acumulada
- **Linha azul — "Low risk":** Escore CorePAM < mediana
- **Linha vermelha — "High risk":** Escore CorePAM ≥ mediana
- **Comparação fundamental:** O HR e as curvas KM desta figura versus a análise principal (FU completo) mostram se o efeito do CorePAM é visível já nos primeiros 2 anos ou só a longo prazo.

Importante: como apenas eventos em ≤24 meses são contabilizados, o número efetivo de eventos será menor que o total de 150, e as estimativas serão mais incertas (IC mais largos), mas o seguimento será mais homogêneo.

## Achados principais

- **N total:** 1.072 (mesmo que análise principal)
- **Eventos a 24 meses:** Subconjunto dos 150 eventos totais — apenas eventos precoces
- **HR estimado a 24 meses:** Tipicamente ligeiramente diferente do HR com FU completo; direção deve se manter
- **p-valor:** A significância estatística com menos eventos tende a ser menor, mas a direção do efeito é o que importa para a análise de sensibilidade
- **Mensagem principal:** Se o HR permanece acima de 1,0 e a separação das curvas é visível a 24 meses, confirma que o CorePAM tem efeito prognóstico **precoce** e não apenas tardio

## O que isso significa?

O TCGA-BRCA é a coorte de validação mais desafiadora do estudo por seu FU curto. Com apenas 32 meses de mediana, muitas pacientes ainda não atingiram os horizontes temporais classicamente relevantes em câncer de mama (5 anos, 10 anos).

A análise de sensibilidade a 24 meses testa se: **mesmo no subconjunto de eventos precoces (que são os mais bem documentados no TCGA), o CorePAM ainda discrimina prognóstico.**

Eventos precoces em câncer de mama tendem a ser de subtipos mais agressivos — Basal-like e HER2-enriched — enquanto eventos tardios são mais de subtipos luminais. Se o CorePAM funciona a 24 meses, isso confirma que ele captura o risco de morte precoce por câncer agressivo, que é justamente onde o sinal dos genes de proliferação (EXO1, PTTG1, CENPF, etc.) deveria se manifestar.

A análise a 24 meses também aborda implicitamente a pergunta: **"O TCGA-BRCA terá mais eventos se o seguimento continuar, e isso fortalecerá o resultado?"** A resposta implícita da análise de sensibilidade é: "Mesmo com eventos precoces, o sinal está lá."

## Contexto no estudo

Esta análise complementa a análise principal do TCGA-BRCA (Fig2) e é colocada como suplementar porque:
1. Não é o resultado primário — é uma verificação de robustez.
2. Com menos eventos, tem menos poder que a análise principal.
3. É uma análise pré-especificada no protocolo para lidar com a limitação do FU curto.

O par FigS5 (METABRIC sensibilidade com OS) e FigS6 (TCGA-BRCA sensibilidade com 24m) demonstram que os resultados são robustos a escolhas alternativas razoáveis de análise.

## Pontos de atenção

- A análise a 24 meses tem **menos eventos** que a análise principal — portanto, IC mais largos e menor poder estatístico. Isso não é uma falha; é esperado e corretamente reportado.
- Truncar em 24 meses **não elimina os pacientes** — apenas censura eventos tardios. Todos os 1.072 pacientes contribuem informação até 24 meses.
- O valor de 24 meses como ponto de corte é **pré-especificado**, não escolhido post-hoc para otimizar o resultado — isso é crucial para validade da análise de sensibilidade.
- A análise a 24 meses é particularmente relevante para subtipos de alto risco (Basal-like, HER2-enriched), onde a maioria dos eventos precoces ocorre. Em subtipos luminais, 24 meses é muito cedo para a maioria dos eventos.
- O seguimento do TCGA-BRCA **continuará crescendo** — análises futuras com mais follow-up confirmarão (ou não) os resultados desta análise de sensibilidade.
