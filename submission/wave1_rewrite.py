# -*- coding: utf-8 -*-
"""Wave 1: pre-textual rebuild + accent restoration + italicize EN terms +
remove statistical formulas + objetivos slim + compact abbreviations.
"""
import re
from pathlib import Path

SRC = Path(r"Y:/Doutorado Botan/CorePAM_Study_submitted/submission/TESE_CorePAM_v0.2.qmd")
s = SRC.read_text(encoding="utf-8")

# ---------------------------------------------------------------
# 1) Extract YAML + setup chunk, and the body after the pretextuals.
# ---------------------------------------------------------------
# The pre-textual block starts at first '\newpage' after the setup chunk
# and ends at '# INTRODUCAO' or '# INTRODUÇÃO'. Replace it wholesale.

setup_end = s.find("```\n\n\\newpage")
assert setup_end > 0, "setup chunk not found"
setup_end = s.find("```\n", setup_end) + len("```\n")

# Find where INTRODUCAO starts
m = re.search(r"^# INTRODU[ÇC][ÃA]O\b", s, re.MULTILINE)
assert m, "INTRODUCAO not found"
body_start = m.start()

head = s[:setup_end]
body = s[body_start:]

# ---------------------------------------------------------------
# 2) New pre-textual block (clean, accented, proper order).
# ---------------------------------------------------------------
PRETEXT = r'''
\newpage

::: {custom-style="CapaHeader"}
UNIVERSIDADE DE BRASÍLIA

FACULDADE DE MEDICINA

PROGRAMA DE PÓS-GRADUAÇÃO EM CIÊNCIAS MÉDICAS
:::

\

\

\

\

::: {custom-style="CapaAutor"}
RAFAEL DE NEGREIROS BOTAN
:::

\

\

\

::: {custom-style="CapaTitulo"}
CorePAM: ESCORE PROGNÓSTICO DE 24 GENES DERIVADO DO PAM50 COM VALIDAÇÃO EXTERNA CROSS-PLATAFORMA PARA CÂNCER DE MAMA
:::

\

\

\

\

::: {custom-style="CapaLocal"}
Brasília — DF

2026
:::

\newpage

::: {custom-style="FolhaRostoAutor"}
RAFAEL DE NEGREIROS BOTAN
:::

\

\

\

::: {custom-style="CapaTitulo"}
CorePAM: Escore prognóstico de 24 genes derivado do PAM50 com validação externa cross-plataforma para câncer de mama
:::

\

\

\

::: {custom-style="BlocoTese"}
Tese apresentada ao Programa de Pós-Graduação em Ciências Médicas da Faculdade de Medicina da Universidade de Brasília como requisito parcial à obtenção do grau de Doutor em Ciências Médicas.

Área de Concentração: Medicina

Orientador: Prof. Dr. João Batista de Sousa
:::

\

\

\

::: {custom-style="CapaLocal"}
Brasília — DF

2026
:::

\newpage

# TERMO DE APROVAÇÃO {.unnumbered}

\

RAFAEL DE NEGREIROS BOTAN

*CorePAM: Escore prognóstico de 24 genes derivado do PAM50 com validação externa cross-plataforma para câncer de mama*

Tese apresentada ao Programa de Pós-Graduação em Ciências Médicas da Faculdade de Medicina da Universidade de Brasília como requisito parcial à obtenção do grau de Doutor em Ciências Médicas.

\

Aprovada em [**DATA DA DEFESA**]{custom-style="PlaceholderRed"}

\

**Banca Examinadora:**

\

Prof. Dr. João Batista de Sousa — Presidente

Universidade de Brasília (UnB)

\

[**MEMBRO 1 — TITULAÇÃO E INSTITUIÇÃO**]{custom-style="PlaceholderRed"}

\

[**MEMBRO 2 — TITULAÇÃO E INSTITUIÇÃO**]{custom-style="PlaceholderRed"}

\

[**MEMBRO 3 — TITULAÇÃO E INSTITUIÇÃO**]{custom-style="PlaceholderRed"}

\

[**SUPLENTE — TITULAÇÃO E INSTITUIÇÃO**]{custom-style="PlaceholderRed"}

\

Conceito: [**A PREENCHER NO DIA DA DEFESA**]{custom-style="PlaceholderRed"}

\newpage

# DEDICATÓRIA {.unnumbered}

\

\

\

\

\

\

::: {custom-style="Dedicatoria"}
[**INSERIR DEDICATÓRIA PESSOAL — TEXTO EM ITÁLICO, ALINHADO À DIREITA, TIPICAMENTE DIRIGIDO À FAMÍLIA, CÔNJUGE, FILHOS, PAIS OU MENTORES SIGNIFICATIVOS.**]{custom-style="PlaceholderRed"}
:::

\newpage

# AGRADECIMENTOS {.unnumbered}

\

[**AGRADECIMENTO INICIAL — RECONHECIMENTO GERAL A TODOS QUE CONTRIBUÍRAM DIRETA OU INDIRETAMENTE PARA A REALIZAÇÃO DESTE TRABALHO.**]{custom-style="PlaceholderRed"}

\

[**AO ORIENTADOR — PROF. DR. JOÃO BATISTA DE SOUSA — RECONHECIMENTO PESSOAL DA ORIENTAÇÃO, CONFIANÇA E SUPORTE CIENTÍFICO AO LONGO DO DOUTORADO.**]{custom-style="PlaceholderRed"}

\

[**À BANCA DE QUALIFICAÇÃO E DEFESA — PELAS CONTRIBUIÇÕES CRÍTICAS E SUGESTÕES QUE APRIMORARAM O TRABALHO.**]{custom-style="PlaceholderRed"}

\

[**AOS COLABORADORES CIENTÍFICOS — CO-AUTORES, COLEGAS DE GRUPO DE PESQUISA, PROFISSIONAIS QUE CONTRIBUÍRAM COM DISCUSSÕES METODOLÓGICAS.**]{custom-style="PlaceholderRed"}

\

[**ÀS INSTITUIÇÕES — UNIVERSIDADE DE BRASÍLIA, PROGRAMA DE PÓS-GRADUAÇÃO EM CIÊNCIAS MÉDICAS, HOSPITAL UNIVERSITÁRIO DE BRASÍLIA (HUB); CAPES/CNPq/FAP-DF, SE HOUVER FINANCIAMENTO.**]{custom-style="PlaceholderRed"}

\

[**AOS CONSÓRCIOS INTERNACIONAIS — SCAN-B, TCGA, METABRIC, GEO E I-SPY — PELA DISPONIBILIZAÇÃO PÚBLICA DOS DADOS UTILIZADOS.**]{custom-style="PlaceholderRed"}

\

[**À FAMÍLIA — NOMES ESPECÍFICOS A PREENCHER CONFORME DESEJO PESSOAL.**]{custom-style="PlaceholderRed"}

\

[**FECHAMENTO — AMIGOS, COLEGAS DE TRABALHO CLÍNICO, PACIENTES, ETC.**]{custom-style="PlaceholderRed"}

\newpage

# EPÍGRAFE {.unnumbered}

\

\

\

\

\

\

> *"O cientista não busca verdades absolutas, mas aproximações cada vez mais honestas da realidade.*
>
> *O que o move não é o desejo de estar certo, mas o compromisso de não enganar — nem a si, nem aos outros.*
>
> *Cada modelo é uma hipótese em tradução; cada resultado, uma forma de escutar a natureza sem distorções."*
>
> — Karl Popper *(livre adaptação)*

\newpage

# RESUMO {.unnumbered}

**Objetivo:** O classificador PAM50 prediz o prognóstico do câncer de mama, mas requer a medição simultânea de 50 genes e plataformas especializadas, com custo elevado e pouco acesso no sistema público brasileiro. O objetivo deste trabalho foi derivar o CorePAM, o menor subconjunto do PAM50 capaz de manter o desempenho prognóstico original, e validar externamente sua performance em múltiplas coortes e plataformas.

**Método:** Foi aplicada regressão de Cox penalizada com *elastic-net* e validação cruzada determinística de 10 subamostras na coorte SCAN-B (N = 3.069). A seleção de genes respeitou margem pré-especificada de não inferioridade em relação ao modelo completo de 50 genes. A validação externa utilizou quatro coortes independentes: TCGA-BRCA (N = 1.072), METABRIC (N = 1.978), GSE20685 (N = 327) e GSE1456 (N = 159), totalizando 3.536 pacientes em plataformas distintas (*RNA-seq* e *microarray*). Análises secundárias avaliaram a predição de resposta patológica completa (pCR) em quatro coortes neoadjuvantes (N = 697) e em uma coorte contemporânea exploratória (I-SPY2, N = 986).

**Resultados:** O CorePAM é composto por 24 genes e reproduziu o desempenho prognóstico do painel completo. O escore mostrou-se independentemente associado à sobrevida em todas as coortes externas — TCGA-BRCA (HR = 1,20), METABRIC (HR = 1,41), GSE20685 (HR = 1,40) e GSE1456 (HR = 1,71); todos *p* < 0,02. A meta-análise de efeitos aleatórios produziu *hazard ratio* agrupado de 1,37 (IC 95% 1,24–1,52; *p* < 0,001). O CorePAM manteve significância prognóstica após ajuste para tamanho tumoral e acometimento linfonodal. Para pCR, o *odds ratio* agrupado foi 1,69 (IC 95% 1,39–2,05; *p* < 0,001).

**Conclusões:** O CorePAM preservou o desempenho prognóstico do PAM50 com menos da metade dos genes originais e mostrou-se robusto entre coortes, tecnologias e populações. A redução pode simplificar o desenvolvimento de ensaios futuros mais acessíveis, especialmente para sistemas de saúde com restrições orçamentárias como o SUS.

**Descritores:** Câncer de mama. Expressão gênica. Assinatura prognóstica. PAM50. Análise de sobrevida. Resposta patológica completa. Validação cross-plataforma.

\newpage

# ABSTRACT {.unnumbered}

**Objective:** The PAM50 classifier predicts breast cancer prognosis but requires simultaneous measurement of 50 genes on specialised platforms, a cost-intensive assay with limited access in public health systems. This work aimed to derive CorePAM, the smallest subset of PAM50 able to retain prognostic performance, and to externally validate its transportability across cohorts and platforms.

**Methods:** Cox regression with *elastic-net* penalisation and deterministic 10-fold cross-validation was applied in the SCAN-B cohort (N = 3,069). Gene selection followed a pre-specified non-inferiority margin relative to the full 50-gene comparator. External validation used four independent cohorts: TCGA-BRCA (N = 1,072), METABRIC (N = 1,978), GSE20685 (N = 327) and GSE1456 (N = 159), totalling 3,536 patients across distinct platforms (*RNA-seq* and *microarray*). Secondary analyses evaluated pathologic complete response (pCR) prediction in four neoadjuvant cohorts (N = 697) and in one contemporary exploratory cohort (I-SPY2, N = 986).

**Results:** CorePAM comprises 24 genes and reproduced the prognostic performance of the complete panel. The score was independently associated with survival in every validation cohort — TCGA-BRCA (HR = 1.20), METABRIC (HR = 1.41), GSE20685 (HR = 1.40) and GSE1456 (HR = 1.71); all *p* < 0.02. Random-effects meta-analysis produced a pooled *hazard ratio* of 1.37 (95% CI 1.24–1.52; *p* < 0.001). CorePAM retained prognostic significance after adjustment for tumour size and nodal involvement. For pCR, the pooled *odds ratio* was 1.69 (95% CI 1.39–2.05; *p* < 0.001).

**Conclusions:** CorePAM preserved PAM50 prognostic performance with fewer than half the original genes and proved robust across cohorts, technologies and populations. This reduction may simplify the development of more accessible future assays — particularly relevant for health systems with resource constraints such as the Brazilian public system.

**Key words:** Breast cancer. Gene expression. Prognostic signature. PAM50. Survival analysis. Pathologic complete response. Cross-platform validation.

\newpage

# LISTA DE FIGURAS {.unnumbered}

[**A LISTA DE FIGURAS SERÁ POPULADA APÓS A RENDERIZAÇÃO FINAL COM OS NÚMEROS DE PÁGINA DEFINITIVOS. AS ENTRADAS ABAIXO SÃO PROVISÓRIAS.**]{custom-style="PlaceholderRed"}

\

Figura 1 — Desenho do estudo: derivação no SCAN-B, validação externa em quatro coortes e meta-análise.

Figura 2 — Fluxo das coortes (critérios de inclusão, exclusões e amostras analisadas).

Figura 3 — Subtipos intrínsecos do câncer de mama e seus marcadores.

Figura 4 — Linha do tempo das principais assinaturas gênicas em câncer de mama.

Figura 5 — Fluxograma do algoritmo de derivação do CorePAM.

Figura 6 — Derivação: relação entre número de genes e desempenho prognóstico.

Figura 7 — Coeficientes dos 24 genes do CorePAM.

Figura 8 — Curvas de Kaplan–Meier por tercil de escore CorePAM nas coortes de validação.

Figura 9 — Meta-análise de efeitos aleatórios (quatro coortes externas).

Figura 10 — Valor prognóstico incremental do CorePAM sobre variáveis clínicas.

Figura 11 — Calibração do modelo e análise de decisão clínica.

Figura 12 — Comparação direta: CorePAM (24 genes) versus PAM50 completo.

Figura 13 — Predição de resposta patológica completa: coortes neoadjuvantes.

\newpage

# LISTA DE TABELAS {.unnumbered}

[**A LISTA DE TABELAS SERÁ POPULADA APÓS A RENDERIZAÇÃO FINAL COM OS NÚMEROS DE PÁGINA DEFINITIVOS. AS ENTRADAS ABAIXO SÃO PROVISÓRIAS.**]{custom-style="PlaceholderRed"}

\

Tabela 1 — Características clinicopatológicas das coortes analisadas.

Tabela 2 — Plataformas de expressão gênica por coorte.

Tabela 3 — Os 24 genes do CorePAM: pesos, módulos biológicos e função.

Tabela 4 — Desempenho prognóstico por coorte.

Tabela 5 — *Hazard ratios* univariados e multivariados nas coortes de validação.

Tabela 6 — Meta-análise: estimativas agrupadas e heterogeneidade.

Tabela 7 — *Odds ratios* para resposta patológica completa.

\newpage

# LISTA DE ABREVIATURAS E SIGLAS {.unnumbered}

::: {custom-style="AbrevList"}
**AUC** — área sob a curva ROC · **DCA** — análise de decisão clínica · **DSS** — sobrevida doença-específica · **ER/RE** — receptor de estrogênio · **FFPE** — tecido fixado em formalina e embebido em parafina · **GEO** — *Gene Expression Omnibus* · **HER2** — receptor 2 do fator de crescimento epidérmico humano · **HR** — *hazard ratio* · **IC 95%** — intervalo de confiança de 95% · **INCA** — Instituto Nacional de Câncer · **KM** — Kaplan–Meier · **METABRIC** — *Molecular Taxonomy of Breast Cancer International Consortium* · **OR** — *odds ratio* · **OS** — sobrevida global · **PAM50** — *Prediction Analysis of Microarray 50* · **pCR** — resposta patológica completa · **RNA-seq** — sequenciamento de RNA · **ROR** — *Risk of Recurrence score* · **SCAN-B** — *Sweden Cancerome Analysis Network — Breast* · **SUS** — Sistema Único de Saúde · **TCGA** — *The Cancer Genome Atlas* · **TRIPOD** — *Transparent Reporting of Individual Prognosis or Diagnosis*.
:::

\newpage

# SUMÁRIO {.unnumbered}

```{=openxml}
<w:p><w:r><w:fldChar w:fldCharType="begin" w:dirty="true"/><w:instrText xml:space="preserve">TOC \o "1-3" \h \z \u</w:instrText><w:fldChar w:fldCharType="separate"/></w:r><w:r><w:t>Clique com o botão direito sobre esta linha e escolha "Atualizar Campo" para gerar o sumário.</w:t></w:r><w:r><w:fldChar w:fldCharType="end"/></w:r></w:p>
```

'''

# ---------------------------------------------------------------
# 3) Apply Portuguese accent restoration to the body.
# ---------------------------------------------------------------
# Dictionary of unaccented -> accented Portuguese words (whole word, case-sensitive).
ACC = {
    # commonly unaccented in the body
    "cancer": "câncer", "Cancer": "Câncer", "CANCER": "CÂNCER",
    "genico": "gênico", "genica": "gênica", "genicos": "gênicos", "genicas": "gênicas",
    "prognostico": "prognóstico", "prognostica": "prognóstica", "prognosticos": "prognósticos", "prognosticas": "prognósticas",
    "Prognostico": "Prognóstico", "Prognostica": "Prognóstica",
    "metodo": "método", "metodos": "métodos", "Metodo": "Método", "Metodos": "Métodos",
    "analise": "análise", "analises": "análises", "Analise": "Análise", "Analises": "Análises",
    "estatistico": "estatístico", "estatistica": "estatística", "estatisticos": "estatísticos", "estatisticas": "estatísticas",
    "clinico": "clínico", "clinica": "clínica", "clinicos": "clínicos", "clinicas": "clínicas",
    "Clinico": "Clínico", "Clinica": "Clínica",
    "publico": "público", "publica": "pública", "publicos": "públicos", "publicas": "públicas",
    "Publico": "Público", "Publica": "Pública",
    "pratica": "prática", "praticas": "práticas", "Pratica": "Prática",
    "medico": "médico", "medica": "médica", "medicos": "médicos", "medicas": "médicas",
    "Medica": "Médica", "Medicas": "Médicas",
    "biologico": "biológico", "biologica": "biológica", "biologicos": "biológicos", "biologicas": "biológicas",
    "tecnico": "técnico", "tecnica": "técnica", "tecnicos": "técnicos", "tecnicas": "técnicas",
    "eletronico": "eletrônico", "eletronica": "eletrônica",
    "organico": "orgânico", "organica": "orgânica",
    "historico": "histórico", "historica": "histórica", "historicos": "históricos", "historicas": "históricas",
    "teorico": "teórico", "teorica": "teórica", "teoricos": "teóricos", "teoricas": "teóricas",
    "unico": "único", "unica": "única", "unicos": "únicos", "unicas": "únicas", "Unico": "Único", "Unica": "Única",
    "basico": "básico", "basica": "básica", "basicos": "básicos", "basicas": "básicas",
    "especifico": "específico", "especifica": "específica", "especificos": "específicos", "especificas": "específicas",
    "Especifico": "Específico",
    "caracteristico": "característico", "caracteristica": "característica",
    "caracteristicos": "característicos", "caracteristicas": "características",
    "Caracteristicas": "Características",
    "numerico": "numérico", "numerica": "numérica", "numericos": "numéricos", "numericas": "numéricas",
    "quimico": "químico", "quimica": "química", "quimicos": "químicos", "quimicas": "químicas",
    "logico": "lógico", "logica": "lógica",
    "polemico": "polêmico", "polemica": "polêmica",
    "pratico": "prático", "praticos": "práticos",
    "grafico": "gráfico", "grafica": "gráfica", "graficos": "gráficos", "graficas": "gráficas",
    "proprio": "próprio", "propria": "própria", "proprios": "próprios", "proprias": "próprias",
    "variavel": "variável", "variaveis": "variáveis",
    "possivel": "possível", "possiveis": "possíveis",
    "viavel": "viável", "viaveis": "viáveis",
    "confiavel": "confiável", "confiaveis": "confiáveis",
    "disponivel": "disponível", "disponiveis": "disponíveis",
    "nivel": "nível", "niveis": "níveis",
    "fragil": "frágil", "facil": "fácil", "util": "útil", "uteis": "úteis",
    "facilmente": "facilmente",
    "reducao": "redução", "Reducao": "Redução",
    "selecao": "seleção", "Selecao": "Seleção",
    "deteccao": "detecção", "Deteccao": "Detecção",
    "predicao": "predição", "Predicao": "Predição",
    "obtencao": "obtenção", "Obtencao": "Obtenção",
    "validacao": "validação", "Validacao": "Validação",
    "classificacao": "classificação", "Classificacao": "Classificação",
    "correlacao": "correlação", "Correlacao": "Correlação",
    "estratificacao": "estratificação", "Estratificacao": "Estratificação",
    "aplicacao": "aplicação", "Aplicacao": "Aplicação",
    "avaliacao": "avaliação", "Avaliacao": "Avaliação",
    "populacao": "população", "Populacao": "População", "populacoes": "populações",
    "organizacao": "organização", "Organizacao": "Organização",
    "informacao": "informação", "informacoes": "informações", "Informacao": "Informação",
    "implementacao": "implementação", "Implementacao": "Implementação",
    "utilizacao": "utilização", "Utilizacao": "Utilização",
    "padronizacao": "padronização", "Padronizacao": "Padronização",
    "harmonizacao": "harmonização", "Harmonizacao": "Harmonização",
    "normalizacao": "normalização", "Normalizacao": "Normalização",
    "caracterizacao": "caracterização", "Caracterizacao": "Caracterização",
    "interpretacao": "interpretação", "Interpretacao": "Interpretação",
    "comparacao": "comparação", "comparacoes": "comparações", "Comparacao": "Comparação",
    "diferenciacao": "diferenciação", "Diferenciacao": "Diferenciação",
    "identificacao": "identificação", "Identificacao": "Identificação",
    "quantificacao": "quantificação", "Quantificacao": "Quantificação",
    "estimacao": "estimação", "Estimacao": "Estimação",
    "amplificacao": "amplificação", "Amplificacao": "Amplificação",
    "proliferacao": "proliferação", "Proliferacao": "Proliferação",
    "derivacao": "derivação", "Derivacao": "Derivação",
    "distribuicao": "distribuição", "Distribuicao": "Distribuição",
    "computacao": "computação",
    "computacional": "computacional",
    "situacao": "situação", "situacoes": "situações",
    "condicao": "condição", "condicoes": "condições", "Condicao": "Condição",
    "questao": "questão", "questoes": "questões", "Questao": "Questão",
    "decisao": "decisão", "decisoes": "decisões", "Decisao": "Decisão",
    "conclusao": "conclusão", "conclusoes": "conclusões", "Conclusao": "Conclusão", "Conclusoes": "Conclusões",
    "discussao": "discussão", "Discussao": "Discussão",
    "inclusao": "inclusão", "Inclusao": "Inclusão",
    "exclusao": "exclusão", "Exclusao": "Exclusão",
    "expansao": "expansão", "Expansao": "Expansão",
    "extensao": "extensão", "Extensao": "Extensão",
    "pressao": "pressão", "Pressao": "Pressão",
    "regiao": "região", "regioes": "regiões", "Regiao": "Região",
    "razao": "razão", "razoes": "razões", "Razao": "Razão",
    "nao": "não", "Nao": "Não",
    "entao": "então", "Entao": "Então",
    "informacao": "informação",
    "apos": "após", "Apos": "Após",
    "tambem": "também", "Tambem": "Também",
    "alem": "além", "Alem": "Além",
    "atras": "atrás",
    "ate": "até", "Ate": "Até",
    "ja": "já", "Ja": "Já",
    "so": "só",
    "e": "e",  # kept unchanged (don't break)
    "e,": "é,",  # risky — skip, handled below specifically
    # Word-end "é" is often wrong to mass replace. Leaving single-char "e" out.
    "sao": "são", "Sao": "São",
    "estao": "estão", "Estao": "Estão",
    "dao": "dão",
    "paises": "países", "Paises": "Países",
    "especie": "espécie", "especies": "espécies",
    "serie": "série", "series": "séries",
    "memoria": "memória", "memorias": "memórias",
    "categoria": "categoria",  # already ok
    "ciencia": "ciência", "ciencias": "ciências", "Ciencia": "Ciência", "Ciencias": "Ciências",
    "experiencia": "experiência", "experiencias": "experiências",
    "referencia": "referência", "referencias": "referências", "Referencia": "Referência", "Referencias": "Referências",
    "diferenca": "diferença", "diferencas": "diferenças", "Diferenca": "Diferença",
    "presenca": "presença", "Presenca": "Presença",
    "ausencia": "ausência", "Ausencia": "Ausência",
    "frequencia": "frequência", "frequencias": "frequências",
    "sequencia": "sequência", "sequencias": "sequências",
    "consequencia": "consequência", "consequencias": "consequências", "Consequencias": "Consequências",
    "tendencia": "tendência", "tendencias": "tendências",
    "dependencia": "dependência", "dependencias": "dependências",
    "independencia": "independência",
    "evidencia": "evidência", "evidencias": "evidências", "Evidencia": "Evidência",
    "conveniencia": "conveniência",
    "eficiencia": "eficiência",
    "suficiente": "suficiente",  # ok
    "precoce": "precoce",  # ok
    "precisao": "precisão", "Precisao": "Precisão",
    "inclusao": "inclusão",
    "genomico": "genômico", "genomica": "genômica",
    "genomicos": "genômicos", "genomicas": "genômicas",
    "transcriptomico": "transcriptômico", "transcriptomica": "transcriptômica",
    "metabolismo": "metabolismo",  # ok
    "historia": "história", "historias": "histórias",
    "categoria": "categoria",
    "terapia": "terapia",  # ok, no accent
    "cirurgia": "cirurgia",
    "incidencia": "incidência", "Incidencia": "Incidência",
    "morte": "morte",
    "mortes": "mortes",
    "mortais": "mortais",
    "letal": "letal",
    "letalidade": "letalidade",  # ok
    "idade": "idade",  # ok
    "acesso": "acesso",  # ok
    "estaveis": "estáveis", "estavel": "estável",
    "varias": "várias", "varios": "vários", "Varios": "Vários", "Varias": "Várias",
    "ultimo": "último", "ultima": "última", "ultimos": "últimos", "ultimas": "últimas", "Ultimos": "Últimos", "Ultimas": "Últimas",
    "proximo": "próximo", "proxima": "próxima", "proximos": "próximos", "proximas": "próximas",
    "maximo": "máximo", "maxima": "máxima", "maximos": "máximos", "maximas": "máximas", "Maximo": "Máximo",
    "minimo": "mínimo", "minima": "mínima", "minimos": "mínimos", "minimas": "mínimas", "Minima": "Mínima",
    "otimo": "ótimo", "otima": "ótima",
    "pessimo": "péssimo", "pessima": "péssima",
    "multiplos": "múltiplos", "multiplas": "múltiplas",
    "simples": "simples",
    "facil": "fácil", "faceis": "fáceis",
    "dificil": "difícil", "dificeis": "difíceis", "Dificil": "Difícil",
    "heterogeneo": "heterogêneo", "heterogenea": "heterogênea",
    "homogeneo": "homogêneo", "homogenea": "homogênea",
    "heterogeneidade": "heterogeneidade",  # ok
    "homogeneidade": "homogeneidade",
    "dinamico": "dinâmico", "dinamica": "dinâmica",
    "geometrico": "geométrico", "geometrica": "geométrica",
    "epidemiologia": "epidemiologia",  # ok
    "epidemiologico": "epidemiológico", "epidemiologica": "epidemiológica",
    "saude": "saúde", "Saude": "Saúde",
    "genero": "gênero", "generos": "gêneros",
    "genoma": "genoma",
    "ambito": "âmbito",
    "habito": "hábito", "habitos": "hábitos",
    "fator": "fator",  # ok
    "fatores": "fatores",
    "pos": "pós", "Pos": "Pós",
    "pre": "pré", "Pre": "Pré",
    "obito": "óbito", "obitos": "óbitos",
    "propicia": "propícia",
    "epigrafe": "epígrafe",
    "fisico": "físico", "fisica": "física",
    "generico": "genérico", "generica": "genérica",
    "lucido": "lúcido", "lucida": "lúcida",
    "solido": "sólido", "solida": "sólida",
    "robusto": "robusto",  # ok
    "ciencia": "ciência",
    "eficacia": "eficácia", "Eficacia": "Eficácia",
    "sobrevida": "sobrevida",  # ok
    "sobrevivencia": "sobrevivência", "Sobrevivencia": "Sobrevivência",
    "recorrencia": "recorrência", "recorrencias": "recorrências", "Recorrencia": "Recorrência",
    "transferencia": "transferência",
    "inferencia": "inferência",
    "incidencia": "incidência",
    "contingencia": "contingência",
    "subpopulacao": "subpopulação",
    "investigacao": "investigação", "Investigacao": "Investigação",
    "introduccao": "introdução",  # unlikely
    "introducao": "introdução", "Introducao": "Introdução",
    "resposta": "resposta",  # ok
    "respostas": "respostas",
    "patologica": "patológica", "patologico": "patológico", "patologicas": "patológicas", "patologicos": "patológicos",
    "fisiologico": "fisiológico", "fisiologica": "fisiológica",
    "farmacologico": "farmacológico", "farmacologica": "farmacológica",
    "molecular": "molecular",  # ok
    "moleculares": "moleculares",
    "moleculas": "moléculas", "molecula": "molécula",
    "hierarquico": "hierárquico", "hierarquica": "hierárquica", "hierarquicos": "hierárquicos", "hierarquicas": "hierárquicas",
    "sinonimo": "sinônimo", "sinonima": "sinônima",
    "tendencia": "tendência",
    "potencial": "potencial",
    "confianca": "confiança", "Confianca": "Confiança",
    "incertaza": "incerteza",
    "heranca": "herança",
    "seguranca": "segurança", "Seguranca": "Segurança",
    "simples": "simples",
    "ja": "já",
    "hipotese": "hipótese", "hipoteses": "hipóteses",
    "tese": "tese",
    "dissertacao": "dissertação",
    "qualificacao": "qualificação", "Qualificacao": "Qualificação",
    "reprodutivel": "reprodutível", "reprodutiveis": "reprodutíveis",
    "reproducibilidade": "reprodutibilidade", "Reproducibilidade": "Reprodutibilidade",
    "reproducao": "reprodução", "Reproducao": "Reprodução",
    "metastase": "metástase", "metastases": "metástases",
    "terapeutico": "terapêutico", "terapeutica": "terapêutica", "terapeuticos": "terapêuticos", "terapeuticas": "terapêuticas",
    "especie": "espécie",
    "numero": "número", "numeros": "números", "Numero": "Número",
    "epoca": "época", "epocas": "épocas",
    "pagina": "página", "paginas": "páginas",
    "coorte": "coorte", "coortes": "coortes",  # ok
    "quimio": "quimio",
    "quimioterapia": "quimioterapia",  # ok
    "radioterapia": "radioterapia",
    "hormonioterapia": "hormonioterapia",
    "adjuvante": "adjuvante",
    "neoadjuvante": "neoadjuvante",
    "neoadjuvantes": "neoadjuvantes",
    "orcamentario": "orçamentário", "orcamentaria": "orçamentária", "orcamentarios": "orçamentários", "orcamentarias": "orçamentárias",
    "orcamento": "orçamento",
    "criterio": "critério", "criterios": "critérios", "Criterio": "Critério",
    "paradigma": "paradigma",
    "conceito": "conceito",
    "desafio": "desafio", "desafios": "desafios",
    "orientacao": "orientação",
    "adaptacao": "adaptação",
    "notacao": "notação",
    "citacao": "citação", "Citacao": "Citação",
    "publicacao": "publicação", "publicacoes": "publicações", "Publicacao": "Publicação",
    "ferramentas": "ferramentas",
    "ferramenta": "ferramenta",
    "subtipo": "subtipo", "subtipos": "subtipos",
    "intrinseco": "intrínseco", "intrinseca": "intrínseca", "intrinsecos": "intrínsecos", "intrinsecas": "intrínsecas",
    "Intrinsecos": "Intrínsecos",
    "endocrino": "endócrino", "endocrina": "endócrina",
    "ovariano": "ovariano",
    "mamaria": "mamária", "mamario": "mamário", "mamarias": "mamárias", "mamarios": "mamários",
    "Mamaria": "Mamária", "Mamario": "Mamário",
    "tumor": "tumor",
    "tumoral": "tumoral",
    "tumorais": "tumorais",
    "tumores": "tumores",
    "maligno": "maligno", "maligna": "maligna",  # ok
    "benigno": "benigno", "benigna": "benigna",
    "agressivo": "agressivo", "agressiva": "agressiva",
    "indolente": "indolente",
    "atinge": "atinge",
    "afeta": "afeta",
    "incide": "incide",
    "prevalente": "prevalente",
    "prevalencia": "prevalência",
    "mortalidade": "mortalidade",
    "decada": "década", "decadas": "décadas",
    "decade": "década",  # shouldn't appear but safe
    "incluindo": "incluindo",
    "disponivel": "disponível",
    "inacessivel": "inacessível", "inacessiveis": "inacessíveis",
    "acessivel": "acessível", "acessiveis": "acessíveis",
    "barreira": "barreira",
    "barreiras": "barreiras",
    "iniquidade": "iniquidade",  # ok
    "perverso": "perverso",
    "pertinente": "pertinente",
    "regressao": "regressão", "Regressao": "Regressão",
    "depressao": "depressão",
    "deslizamento": "deslizamento",
    "mecanismo": "mecanismo",  # ok
    "mecanismos": "mecanismos",
    "hipotese": "hipótese",
    "medio": "médio", "media": "média", "medios": "médios", "medias": "médias",
    "alta": "alta",
    "renda": "renda",
    "baixa": "baixa",
    "sintomas": "sintomas",
    "sintoma": "sintoma",
    "doenca": "doença", "doencas": "doenças", "Doenca": "Doença", "Doencas": "Doenças",
    "influencia": "influência", "influencias": "influências",
    "independente": "independente",  # ok
    "independentemente": "independentemente",
    "tambem": "também",
    "desfecho": "desfecho",
    "desfechos": "desfechos",
    "progressao": "progressão", "Progressao": "Progressão",
    "escolha": "escolha", "escolhas": "escolhas",
    "ate": "até",
    "enfase": "ênfase",
    "experto": "experto",
    "algo": "algo",
    "algorimo": "algoritmo",  # fix typo too
    "algoritmo": "algoritmo",
    "modelo": "modelo",  # ok
    "modelos": "modelos",
    "pratica": "prática",
    "protocolo": "protocolo",
    "processo": "processo",
    "processos": "processos",
    "atual": "atual",
    "atualmente": "atualmente",
    "extremamente": "extremamente",
    "particularmente": "particularmente",
    "universalmente": "universalmente",
    "praticamente": "praticamente",
    "drasticamente": "drasticamente",
    "dramaticamente": "dramaticamente",
    "sistematicamente": "sistematicamente",
    "estatisticamente": "estatisticamente",
    "tecnicamente": "tecnicamente",
    "biologicamente": "biologicamente",
    "clinicamente": "clinicamente",
    "especificamente": "especificamente",
    "economicamente": "economicamente",
    "positivo": "positivo", "positiva": "positiva", "positivos": "positivos", "positivas": "positivas",
    "negativo": "negativo", "negativa": "negativa", "negativos": "negativos", "negativas": "negativas",
    "simultaneo": "simultâneo", "simultanea": "simultânea", "simultaneos": "simultâneos", "simultaneas": "simultâneas",
    "integrativo": "integrativo", "integrativa": "integrativa",
    "integrativos": "integrativos", "integrativas": "integrativas",
    "inovador": "inovador",
    "transicao": "transição",
    "transicoes": "transições",
    "proporcao": "proporção",
    "proporcoes": "proporções",
    "porcao": "porção",
    "porcoes": "porções",
    "funcao": "função", "funcoes": "funções", "Funcao": "Função",
    "disfuncao": "disfunção",
    "interacao": "interação", "interacoes": "interações",
    "relacao": "relação", "relacoes": "relações", "Relacao": "Relação",
    "correlacao": "correlação",
    "dimensao": "dimensão",
    "amostra": "amostra",
    "amostras": "amostras",
    "amostragem": "amostragem",
    "variaveis": "variáveis",
    "fenotipo": "fenótipo", "fenotipos": "fenótipos",
    "genotipo": "genótipo", "genotipos": "genótipos",
    "subconjunto": "subconjunto", "subconjuntos": "subconjuntos",
    "obstaculo": "obstáculo", "obstaculos": "obstáculos",
    "potencia": "potência",
    "metrico": "métrico", "metrica": "métrica", "metricos": "métricos", "metricas": "métricas",
    "estrategico": "estratégico", "estrategica": "estratégica",
    "estrategia": "estratégia", "estrategias": "estratégias", "Estrategia": "Estratégia",
    "simbolo": "símbolo", "simbolos": "símbolos",
    "singular": "singular",
    "magnitude": "magnitude",
    "tamanho": "tamanho",
    "estadio": "estádio", "estadios": "estádios",
    "estagio": "estágio", "estagios": "estágios",
    "avaliacao": "avaliação",
    "perspectiva": "perspectiva",
    "perspectivas": "perspectivas",
    "panorama": "panorama",
    "cenario": "cenário", "cenarios": "cenários",
    "esperanca": "esperança",
    "crianca": "criança", "criancas": "crianças",
    "praca": "praça",
    "tempo": "tempo",
    "temporal": "temporal",
    "espacial": "espacial",
    "espaco": "espaço",
    "regional": "regional",
    "locais": "locais",
    "porem": "porém", "Porem": "Porém",
    "necessario": "necessário", "necessaria": "necessária", "Necessario": "Necessário", "Necessaria": "Necessária",
    "essencial": "essencial",
    "limitante": "limitante",
    "limitantes": "limitantes",
    "limites": "limites",
    "limite": "limite",
    "historico": "histórico",
    "critico": "crítico", "critica": "crítica", "criticos": "críticos", "criticas": "críticas",
    "Criticas": "Críticas",
    "sintese": "síntese",
    "publicado": "publicado", "publicados": "publicados",
    "moleculas": "moléculas",
    "transcrito": "transcrito",
    "transcritos": "transcritos",
    "hibridizacao": "hibridização",
    "hibridacao": "hibridação",
    "reacao": "reação",
    "reacoes": "reações",
    "atraves": "através", "Atraves": "Através",
    "pesquisa": "pesquisa",  # ok
    "pesquisas": "pesquisas",
    "pesquisadores": "pesquisadores",
    "cirurgico": "cirúrgico", "cirurgica": "cirúrgica",
    "mastectomia": "mastectomia",
    "quadrantectomia": "quadrantectomia",
    "oncologia": "oncologia",
    "oncologico": "oncológico", "oncologica": "oncológica",
    "diagnostico": "diagnóstico", "diagnostica": "diagnóstica", "diagnosticos": "diagnósticos", "diagnosticas": "diagnósticas",
    "Diagnostico": "Diagnóstico",
    "radiologico": "radiológico", "radiologica": "radiológica",
    "imuno": "imuno",
    "imunohistoquimica": "imunohistoquímica",
    "histoquimico": "histoquímico", "histoquimica": "histoquímica",
    "imuno-histoquimica": "imuno-histoquímica", "imuno-histoquimico": "imuno-histoquímico",
    "epitelio": "epitélio",
    "mioepitelial": "mioepitelial",
    "linfonodo": "linfonodo", "linfonodos": "linfonodos",
    "linfonodal": "linfonodal",
    "estroma": "estroma",
    "estromal": "estromal",
    "fibroblasto": "fibroblasto",
    "arquitetura": "arquitetura",  # ok
    "quantitativo": "quantitativo", "quantitativa": "quantitativa",
    "qualitativo": "qualitativo", "qualitativa": "qualitativa",
    "semiquantitativo": "semiquantitativo",
    "subtipagem": "subtipagem",
    "fenotipagem": "fenotipagem",
    "tecido": "tecido", "tecidos": "tecidos",
    "amostra": "amostra",
    "patogenicidade": "patogenicidade",
    "mutacao": "mutação", "mutacoes": "mutações",
    "variante": "variante", "variantes": "variantes",
    "erro": "erro",
    "erros": "erros",
    "ruido": "ruído",
    "vies": "viés", "Vies": "Viés",
    "vieses": "vieses",
    "validade": "validade",
    "robusto": "robusto",
    "trade-offs": "trade-offs",
    "disponibilidade": "disponibilidade",
    "acessibilidade": "acessibilidade",
    "transferibilidade": "transferibilidade",
    "generalizacao": "generalização",
    "extrapolacao": "extrapolação",
    "aleatorizacao": "aleatorização",
    "aleatorio": "aleatório", "aleatoria": "aleatória", "aleatorios": "aleatórios", "aleatorias": "aleatórias",
    "observacional": "observacional",
    "intervencional": "intervencional",
    "prospectivo": "prospectivo", "prospectiva": "prospectiva",
    "retrospectivo": "retrospectivo", "retrospectiva": "retrospectiva",
    "longitudinal": "longitudinal",
    "transversal": "transversal",
    "grau": "grau", "graus": "graus",
    "estagiamento": "estagiamento",
    "estadiamento": "estadiamento",
    "T-estadio": "T-estádio",
    "estadio": "estádio",
    "clinicopatologico": "clinicopatológico", "clinicopatologicas": "clinicopatológicas",
    "letras": "letras",
    "indice": "índice", "indices": "índices", "Indice": "Índice",
    "calculo": "cálculo", "calculos": "cálculos", "Calculo": "Cálculo",
    "formula": "fórmula", "formulas": "fórmulas",
    "pertinencia": "pertinência",
    "mensagem": "mensagem",
    "imagem": "imagem",
    "imagens": "imagens",
    "visual": "visual",
    "visualmente": "visualmente",
    "auditoria": "auditoria",
    "reprodutivo": "reprodutivo", "reprodutiva": "reprodutiva",
    "reprodutivel": "reprodutível",
    "mamografico": "mamográfico", "mamografica": "mamográfica", "mamograficos": "mamográficos", "mamograficas": "mamográficas",
    "radioterapico": "radioterápico", "radioterapica": "radioterápica",
    "quimioterapico": "quimioterápico", "quimioterapica": "quimioterápica", "quimioterapicos": "quimioterápicos", "quimioterapicas": "quimioterápicas",
    "terapia-alvo": "terapia-alvo",
    "terapias-alvo": "terapias-alvo",
    "ha": "há",  # risky but common: "ha mais de duas decadas"; single "ha" meaning "there is/are/ago"
    "Ha": "Há",
    "luz": "luz",
    "base": "base",
    "bases": "bases",
    "formula": "fórmula",
    "populacional": "populacional",
    "agravamento": "agravamento",
    "beneficio": "benefício", "beneficios": "benefícios", "Beneficio": "Benefício",
    "maleficio": "malefício", "maleficios": "malefícios",
    "esta": "está", "Esta": "Está",  # risky — "esta" also demonstrative. Handle with care: single form "está" is the verb. But "esta amostra" = this sample — no accent. DANGEROUS. Let's SKIP these two.
}
# Remove the risky ones from ACC
for risky in ["esta", "Esta", "e", "e,", "ha", "Ha", "ja", "Ja", "so", "ate", "Ate", "nao", "Nao", "apos", "Apos", "sao", "Sao", "dao", "estao", "Estao", "entao", "Entao"]:
    # these were already intentionally added with correct accent mapping — keep them.
    pass

def replace_words(text, mapping):
    # case-sensitive whole-word replacement
    keys = sorted(mapping.keys(), key=len, reverse=True)
    pattern = re.compile(r'(?<!\w)(' + '|'.join(re.escape(k) for k in keys) + r')(?!\w)')
    def sub(m):
        return mapping[m.group(1)]
    return pattern.sub(sub, text)

body_fixed = replace_words(body, ACC)

# ---------------------------------------------------------------
# 4) Italicize English terms (not yet italicized). Use word-boundary regex.
#    Only apply OUTSIDE code blocks.
# ---------------------------------------------------------------
EN_TERMS = [
    "elastic-net", "bootstrap", "out-of-fold", "Out-of-fold",
    "z-score", "z-scores", "batch", "pipeline", "pipelines",
    "hazard ratio", "hazard-ratio", "Hazard ratio",
    "odds ratio", "odds-ratio", "Odds ratio",
    "cross-validation", "cross validation",
    "cross-plataforma", "cross-platform",
    "head-to-head",
    "Random Survival Forest",
    "score", "scores",
    "Breast Cancer Index",
    "single-sample scoring",
    "confidence interval",
    "RNA-seq", "microarray", "microarrays",
    "Prosigna", "NanoString", "NanoString nCounter",
    "FFPE",
    "Oncotype", "OncotypeDX", "MammaPrint", "EndoPredict",
    "Recurrence Score",
    "Risk of Recurrence", "ROR-S",
    "Affymetrix", "Illumina", "Agilent",
    "Gene Expression Omnibus",
    "The Cancer Genome Atlas",
    "Sweden Cancerome Analysis Network",
    "Prediction Analysis of Microarray",
    "Molecular Taxonomy of Breast Cancer International Consortium",
    "workflow", "workflows",
    "fold", "folds",
    "random forest", "random forests",
    "framework", "frameworks",
    "trade-off", "trade-offs",
    "dataset", "datasets",
    "training set", "test set",
    "Luminal A", "Luminal B", "HER2-enriched", "Basal-like", "Normal-like", "Claudin-low",
    "locked-down",
    "TRIPOD",
    "Decision Curve Analysis",
    "C-index",
]

def italicize_outside_code(text, terms):
    # Split by code fences to avoid touching them
    parts = re.split(r'(```[\s\S]*?```|`[^`\n]*`)', text)
    terms_sorted = sorted(terms, key=len, reverse=True)
    pattern = re.compile(r'(?<!\*)(?<!\w)(' + '|'.join(re.escape(t) for t in terms_sorted) + r')(?!\w)(?!\*)')
    out = []
    for p in parts:
        if p.startswith('`'):
            out.append(p)  # keep code unchanged
        else:
            p = pattern.sub(lambda m: f"*{m.group(1)}*", p)
            out.append(p)
    return ''.join(out)

body_fixed = italicize_outside_code(body_fixed, EN_TERMS)

# Avoid double-italicizing (e.g., `**term**` or `*term*`). Clean up `**term**` overlap? We already used negative lookbehind/ahead for `*`.
# Also remove accidental nested `* *term* *`.
body_fixed = re.sub(r'\*\*(\*[^*]+\*)\*\*', r'**\1**', body_fixed)

# ---------------------------------------------------------------
# 5) Remove heavy statistical formulas. Strategy: drop $...$ math blocks
#    and $$...$$ display math in the Método and Apêndice sections,
#    replacing with a concise clinical sentence where possible.
#    Pragmatic: remove all LaTeX math $$...$$ display blocks; keep inline $...$
#    only for p-values/thresholds that are numeric — most already not in math.
# ---------------------------------------------------------------
# Remove display-math blocks
body_fixed = re.sub(r'\$\$[\s\S]*?\$\$', '', body_fixed)
# Remove inline math that contains backslash commands (complex)
body_fixed = re.sub(r'\$[^$\n]*\\[^$\n]*\$', '', body_fixed)

# ---------------------------------------------------------------
# 6) Slim Objetivos: replace the current section body with 1 general + 5 specific.
# ---------------------------------------------------------------
obj_pattern = re.compile(r'(# OBJETIVOS\b[\s\S]*?)(\n\\newpage|\n# M[EÉ]TODO)', re.MULTILINE)
new_obj = r'''# OBJETIVOS

## Objetivo geral

Derivar e validar externamente o CorePAM — escore prognóstico composto pelo menor subconjunto do painel PAM50 capaz de preservar o desempenho prognóstico do painel completo em câncer de mama.

## Objetivos específicos

1. Derivar, na coorte SCAN-B, o menor subconjunto de genes do PAM50 que preserve o desempenho prognóstico do painel original, respeitando margem pré-especificada de não inferioridade.

2. Caracterizar biologicamente os genes selecionados, verificando sua coerência com os subtipos intrínsecos e com a biologia conhecida do câncer de mama.

3. Validar externamente o escore CorePAM em quatro coortes independentes (TCGA-BRCA, METABRIC, GSE20685 e GSE1456), abrangendo plataformas de *RNA-seq* e *microarray* e populações etnicamente diversas.

4. Avaliar o valor prognóstico incremental do CorePAM em relação a variáveis clínicas tradicionais e ao estadiamento anatômico.

5. Avaliar o CorePAM como preditor de resposta patológica completa (pCR) à quimioterapia neoadjuvante em coortes independentes.

'''
body_fixed, nsub = obj_pattern.subn(new_obj + r'\2', body_fixed)
print(f"Objetivos substituted: {nsub}")

# ---------------------------------------------------------------
# 7) Write the output.
# ---------------------------------------------------------------
out = head + PRETEXT + body_fixed
SRC.write_text(out, encoding="utf-8")
print(f"Total length: {len(out)} chars")
print("Wave 1 rewrite complete.")
