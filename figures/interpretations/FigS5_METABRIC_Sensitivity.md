# Interpretação: FigS5 — Análise de Sensibilidade do METABRIC: OS em vez de DSS

## O que é este gráfico?

Esta figura suplementar apresenta a **análise de sensibilidade do endpoint** para a coorte METABRIC, mostrando curvas de Kaplan-Meier e estatísticas de sobrevida usando **Sobrevida Global (OS — Overall Survival)** como endpoint alternativo, em vez do endpoint primário **DSS (Disease-Specific Survival)**.

Esta análise existe porque o endpoint primário do METABRIC foi deliberadamente definido como DSS (óbitos apenas por câncer de mama, com censura de óbitos por outras causas). Porém, nem todos os leitores, revisores ou bases de dados registram a causa de morte com precisão. Ao demonstrar que os resultados também são significativos quando se usa OS (que inclui todos os óbitos independentemente de causa), confirmamos que a escolha do endpoint não é "conveniente" para os resultados.

## Como ler este gráfico

O formato é idêntico à figura principal KM do METABRIC (Fig2), mas com endpoint OS em vez de DSS:
- **Eixo X:** Tempo em meses
- **Eixo Y:** Probabilidade de sobrevivência global (OS)
- **Linha azul — "Low risk":** Escore CorePAM < mediana
- **Linha vermelha — "High risk":** Escore CorePAM ≥ mediana
- **HR, IC 95%, p-valor:** Estatísticas da análise Cox com OS como evento

A comparação lado a lado com a figura primária (DSS) permite avaliar:
- Se a separação das curvas se mantém com OS
- Se o HR muda substancialmente (sugerindo confundimento por causas competitivas de morte)
- Se o p-valor permanece significativo

## Achados principais

- **Endpoint:** OS (todos os óbitos, independente de causa)
- **N total:** 1.978 (mesmo que análise primária DSS)
- **Eventos OS:** maior que 646 (DSS), pois inclui óbitos por outras causas — METABRIC tem FU de ~160 meses, então mortes por outras causas são numerosas
- **HR_OS (sensibilidade):** HR = 1,23 (IC 95%: 1,16–1,30)
- **Comparação com HR_DSS primário:** HR_DSS = 1,41 (IC 95%: 1,31–1,52)

O HR com OS (1,23) é **menor** que com DSS (1,41), o que é **esperado e biologicamente correto**: ao incluir óbitos por outras causas (que o CorePAM — um marcador de câncer de mama — não prevê), "dilui-se" o sinal, reduzindo o HR aparente. Esse padrão confirma que o CorePAM é específico para prognóstico de câncer de mama, não um marcador geral de fragilidade.

## O que isso significa?

Esta análise de sensibilidade é metodologicamente elegante porque demonstra **robustez à escolha do endpoint**. O raciocínio completo é:

1. O endpoint primário (DSS) é mais preciso biologicamente — mede exatamente o que queremos (morte por câncer de mama).
2. Porém, a classificação de causa de morte é subjetiva e pode ter erros, especialmente em coortes históricas.
3. Se o modelo funciona também com OS (mais objetivo — morte é morte), confirmamos que os resultados não dependem de uma classificação de causa de morte potencialmente imprecisa.
4. O fato de o HR_OS ser menor que HR_DSS confirma que o CorePAM é mais específico para morte por câncer do que para mortalidade geral — exatamente o que se espera.

O intervalo de confiança do HR_OS (1,16–1,30) ainda está muito longe de 1,0, com p-valor extremamente significativo, demonstrando robustez.

## Contexto no estudo

Esta análise faz parte das **análises de sensibilidade pré-especificadas** do protocolo do estudo. Em pesquisa clínica rigorosa, as análises de sensibilidade são planejadas antes da coleta/análise dos dados e testam se as conclusões principais mudam quando suposições metodológicas são alteradas de forma razoável.

Para o METABRIC especificamente, as análises de sensibilidade incluem:
1. OS vs. DSS (esta figura — FigS5)
2. Quartis vs. mediana para estratificação KM (FigS suplementar de quartis)
3. Análise restrita a subgrupos (ER+, ER−) se relevante

O fato de a análise primária (DSS) mostrar HR = 1,41 e a sensibilidade (OS) mostrar HR = 1,23 — ambos significativos, ambos na mesma direção — fortalece a credibilidade do resultado principal.

## Pontos de atenção

- A diferença entre HR_DSS = 1,41 e HR_OS = 1,23 é explicada pelo fenômeno de **riscos competitivos (competing risks)**: em uma coorte com seguimento de 13+ anos, muitas pacientes morrem de outras causas (cardiovascular, outras neoplasias, causas relacionadas à idade), e o CorePAM não prevê essas mortes. Com OS, esses eventos "competem" com o sinal do câncer de mama.
- Para análises mais rigorosas de riscos competitivos, modelos de Fine-Gray ou análise de incidência cumulativa seriam o padrão — isso vai além do escopo desta análise de sensibilidade mas é uma limitação reconhecível.
- A coorte METABRIC inclui pacientes de diferentes eras de tratamento, e a causa de morte pode ter sido registrada com diferentes padrões em diferentes décadas e centros — isso é uma limitação inerente da coorte, não do estudo.
