# Interpretação: FigS — Gráfico Lollipop: Pesos dos Genes do CorePAM

## O que é este gráfico?

O gráfico **lollipop** (literalmente, "pirulito") é uma versão elegante e legível do gráfico de barras para mostrar os **pesos (coeficientes) de cada gene no modelo CorePAM**. Cada gene do CorePAM é representado por um "palito" horizontal com um círculo na ponta, onde:
- O **comprimento do palito** representa a **magnitude do peso** (quanto aquele gene contribui para o escore).
- A **posição do círculo** (à esquerda ou direita do zero) representa o **sinal do peso** (positivo = associado a mau prognóstico; negativo = associado a bom prognóstico).

Este gráfico é a representação visual mais direta e intuitiva dos 24 genes do CorePAM e seus respectivos pesos derivados pelo elastic-net.

## Como ler este gráfico

- **Eixo X (horizontal):** Valor do coeficiente/peso do modelo. O zero está no centro.
  - **Direita do zero (valores positivos):** Genes com peso positivo — maior expressão → maior escore → maior risco previsto.
  - **Esquerda do zero (valores negativos):** Genes com peso negativo — maior expressão → menor escore → menor risco previsto.
- **Eixo Y (vertical):** Lista de genes CorePAM, tipicamente ordenados do maior peso positivo ao maior peso negativo (ou por magnitude absoluta).
- **Círculo no topo de cada palito:** O ponto exato do peso; pode ser colorido por categoria biológica.
- **Cores dos círculos:** Tipicamente distinguem genes por pathway/categoria:
  - Vermelho/laranja: genes associados a mau prognóstico (basal-like, proliferação)
  - Azul/verde: genes associados a bom prognóstico (luminal, hormônio-receptores)
  - Ou coloridos por subtipo PAM50 de origem (ex: genes "Luminal A", "Basal-like", "HER2-enriched")

## Achados principais

**Genes com peso POSITIVO** (expressão alta = maior risco):
- **EXO1:** Exonuclease 1 — reparo de DNA, proliferação
- **PTTG1:** Securing/PTTG1 — regulação do ciclo celular, anáfase
- **FOXC1:** Fator de transcrição — marcador basal-like
- **PHGDH:** Fosfoglicerato desidrogenase — biossíntese de serina, metabolismo
- **MYC:** Proto-oncogene c-Myc — proliferação, transcripção
- **KRT5:** Citoqueratina 5 — marcador basal-epitelial
- **CXXC5:** Proteína de dedo de zinco — regulação transcricional
- **CENPF:** Proteína centromérica F — divisão celular
- **FGFR4:** Receptor de FGF tipo 4 — sinalização de crescimento
- **ERBB2:** HER2/neu — receptor tirosina-quinase, sinalização

**Genes com peso NEGATIVO** (expressão alta = menor risco):
- **ESR1:** Receptor de estrogênio α — marcador luminal fundamental
- **PGR:** Receptor de progesterona — marcador luminal
- **BCL2:** Proto-oncogene bcl-2 — antiapoptótico, marcador luminal
- **MLPH:** Melanofillina — associado ao Luminal A
- **GRB7:** Proteína adaptadora de sinalização
- **NAT1:** N-acetiltransferase 1 — metabolismo, marcador luminal
- **BLVRA:** Biliverdina redutase A — antioxidante, luminal
- **ACTR3B:** Subunidade do complexo Arp2/3
- **MIA:** Inibidor de migração de melanócitos
- **MYBL2:** Fator de transcrição
- **MDM2:** Regulador negativo de p53
- **SFRP1:** Proteína secretada frizzled-related 1 — antagonista de Wnt
- **GPR160:** Receptor acoplado a proteína G
- **KRT17:** Citoqueratina 17

## O que isso significa?

O gráfico lollipop revela a **biologia capturada pelo CorePAM**:

**Polo de alto risco (pesos positivos):** Dominado por genes de **proliferação** (PTTG1, CENPF, EXO1 — todos envolvidos em ciclo celular e divisão), genes de **fenótipo basal-like** (KRT5, FOXC1) e o **oncogene HER2** (ERBB2). A inclusão de ERBB2 com peso positivo explica por que tumores HER2-enriched tendem a ter escores mais altos.

**Polo de baixo risco (pesos negativos):** Dominado pelos **receptores hormonais** (ESR1, PGR) e seus genes co-expressos (BCL2, MLPH, NAT1) — o assinatura clássica do subtipo Luminal A. Genes como SFRP1 (antagonista de Wnt) e MDM2 (regulador de p53) adicionam informação além da simples dicotomia ER+/ER−.

A **magnitude dos pesos** indica a importância relativa de cada gene. Genes com palitos mais longos (maiores pesos absolutos) têm maior influência na discriminação prognóstica. Em um painel regulado por elastic-net (α=0,5), os pesos são shrunk mas não zerados (ao contrário do Lasso puro), garantindo que todos os 24 genes contribuam de forma ponderada.

A presença de genes como **MYC** (proto-oncogene clássico) e **PHGDH** (reprogramação metabólica) com pesos positivos indica que o CorePAM captura dimensões da biologia do câncer além da simples proliferação — incluindo alterações metabólicas associadas a fenótipos agressivos.

## Contexto no estudo

Esta figura é a "identidade molecular" do CorePAM — ela mostra quem são os 24 genes e qual é o papel de cada um. Para um leitor especialista em biologia molecular de câncer de mama, cada gene conta uma história:
- **ESR1 com peso negativo:** Câncer luminal, receptor de estrogênio alto → melhor prognóstico.
- **KRT5 com peso positivo:** Citoqueratina basal-epitelial → tumor basal-like, pior prognóstico.
- **ERBB2 com peso positivo:** HER2-amplificação → tumor agressivo (porém tratável com trastuzumabe!).

O gráfico lollipop transforma coeficientes numéricos em narrativa biológica coerente — é literalmente a tradução visual de "o que o CorePAM detecta".

## Pontos de atenção

- Os pesos mostrados são os **coeficientes do elastic-net no caminho de regularização** para o λ selecionado. Eles incluem shrinkage e não são equivalentes aos coeficientes de um modelo Cox não-penalizado — não devem ser interpretados como "HRs por unidade de Z-score".
- A **escala relativa** entre os pesos é informativa (gene A tem maior peso que gene B), mas a escala absoluta é arbitrária (depende de λ).
- Os pesos foram estimados em **SCAN-B (RNA-seq)** — a transferência desses pesos para microarray implica que os genes têm a mesma direção de associação prognóstica em ambas as plataformas, o que é confirmado pelas validações mas não necessariamente pelos valores absolutos dos pesos.
- **MDM2 com peso negativo** pode surpreender (MDM2 é classicamente considerado oncogênico), mas em câncer de mama há evidência de que amplificação de MDM2 ocorre em tumores luminais com p53 selvagem — explicando a associação com prognóstico mais favorável neste contexto.
- A cor dos círculos, se usada por categoria biológica, é interpretativa e não derivada algoritmicamente — reflete o julgamento dos autores sobre a função conhecida de cada gene.
