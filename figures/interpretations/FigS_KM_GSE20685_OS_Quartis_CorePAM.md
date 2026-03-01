# Interpretação: FigS — Curvas KM por Quartis do Escore CorePAM — GSE20685

## O que é este gráfico?

Esta figura suplementar apresenta as curvas de Kaplan-Meier para a coorte GSE20685 estratificadas por **quartis do escore CorePAM**, usando sobrevida global (OS) como endpoint. É a análise por quartis da menor coorte de validação do estudo.

Com N=327 pacientes divididos em quatro grupos de aproximadamente 82 pacientes cada e 83 eventos totais (proporção de ~20 eventos por quartil em média), esta é a análise com menor poder estatístico do conjunto. No entanto, a interpretação do gradiente biológico e da consistência com as outras coortes é o que importa aqui.

## Como ler este gráfico

- **Quatro curvas KM** representando os quartis Q1 (menor escore, melhor prognóstico previsto) a Q4 (maior escore, pior prognóstico previsto).
- **Eixo X:** Tempo em meses, com horizonte longo (>100 meses / >8 anos) — refletindo o longo seguimento desta coorte sul-coreana.
- **Número de pacientes em risco por grupo:** Com apenas ~82 por quartil, as curvas serão progressivamente menos precisas à medida que o tempo avança.
- **IC sombreado:** As faixas de confiança serão substancialmente mais largas que nas coortes maiores — isso é esperado e correto estatisticamente.

## Achados principais

Com N=327 e 83 eventos:
- Aproximadamente **20-21 eventos por quartil** em média (se uniformemente distribuídos, o que não ocorre — Q4 terá mais eventos e Q1 menos).
- O poder estatístico para separar quartis adjacentes é limitado — é possível que Q2 vs. Q3 não atinja significância, mas Q1 vs. Q4 deve ser significativo.
- **Gradiente esperado:** Curva Q1 acima, Q4 abaixo, com Q2 e Q3 intermediárias.
- **FU mediano de ~112 meses** permite observar a separação em perspectiva de longo prazo, o que é um ponto positivo dado o pequeno N.

## O que isso significa?

A análise por quartis do GSE20685 testa o gradiente de risco na menor e mais desafiante das três coortes de validação. Mesmo com N pequeno, observar uma separação ordenada (Q1 > Q2 > Q3 > Q4) seria evidência poderosa de que o CorePAM captura variação biológica real nesta coorte — independentemente da limitação de potência.

Para uma coorte de 327 pacientes, o gradiente dose-resposta em um marcador prognóstico pode ser difícil de demonstrar estatisticamente ao nível de quartis individuais, mas se a hierarquia das curvas (Q1 sempre acima de Q4) for visualmente clara e mantida ao longo de todo o seguimento, isso é biologicamente informativo.

O GSE20685 representa uma população sul-coreana com câncer de mama, potencialmente com distribuição de subtipos ligeiramente diferente da europeia/americana. Se o gradiente do CorePAM se mantém nesta população, isso amplia a generalização transcultural do painel.

## Contexto no estudo

As três figuras de quartis suplementares (TCGA-BRCA, METABRIC, GSE20685) devem ser interpretadas **em conjunto** como um padrão de evidência consistente: em todas as coortes, independentemente do N, plataforma, ou população, o escore CorePAM produz um gradiente ordenado de risco. Isso é a melhor evidência de que o sinal prognóstico é real e robusto.

A figura do GSE20685 por quartis é a mais incerta individualmente, mas ganha força pelo contexto das outras duas.

## Pontos de atenção

- Com ~20 eventos por quartil, os **testes log-rank pairwise** entre quartis adjacentes (Q1 vs. Q2 ou Q3 vs. Q4) têm poder muito baixo. Resultados não significativos nesses pares específicos não devem ser interpretados como ausência de efeito.
- A comparação **Q1 vs. Q4** é a única comparação com poder razoável (~20 vs. ~20 eventos dos grupos extremos — ainda limitado).
- A **variabilidade das curvas** nos tempos tardios (>100 meses) será grande — as curvas podem se cruzar ou apresentar comportamentos irregulares por acaso estatístico, não por biologia real.
- O METABRIC e o TCGA-BRCA por quartis são mais informativos individualmente — o GSE20685 contribui principalmente para o padrão de consistência cross-coorte.
- Dado o N pequeno, os IC 95% das curvas de Greenwood serão largos — especialmente para Q1 e Q4, que tendem a ter as estimativas mais incertas por serem os grupos "extremos" com comportamentos eventualmente divergentes.
