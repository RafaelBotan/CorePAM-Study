# Registro de Correções — Tese CorePAM (v1.0 → v1.1)

**Autor:** Rafael de Negreiros Botan
**Documento:** controle de revisão da tese *"CorePAM: escore prognóstico de 24 genes derivado do PAM50 com validação externa cross-plataforma para câncer de mama"*
**Arquivo-fonte corrigido:** `TESE_CorePAM_v1.1_corrigida.qmd`
**Arquivo original preservado:** `TESE_CorePAM_FINAL_envio_Dr_Joao_impressao.qmd` (v1.0, intacto)
**Data da revisão:** 04/06/2026

---

## Sumário executivo

Revisão sistemática da tese (texto integral + auditoria figura-a-figura + checagem de coerência de valores + revisão metodológica). Resultado:

- **Correções de texto aplicadas na v1.1:** 11 (detalhadas na Parte A).
- **Itens de figura pendentes** (exigem regeneração da PNG ou confirmação no script — **não corrigíveis só no texto**): 7 principais + ajustes menores (Parte B).
- **Vulnerabilidades metodológicas** (não são erros; preparação para arguição): 6 (Parte C).

**Natureza das correções:** todas são de **redação, documentação ou apresentação**, mais **um fato externo** (política do SUS) que era factualmente incorreto e foi removido. **Nenhuma correção altera a metodologia, as análises estatísticas ou os resultados quantitativos.** As tabelas de resultados foram conferidas número a número e são internamente consistentes. A tese permanece **válida e robusta** — o artigo correspondente foi revisado por pares e **aceito no *Breast Cancer Research***.

---

## Parte A — Correções de texto aplicadas (v1.1)

| ID | Local (v1.0) | Antes | Depois | Tipo | Impacto |
|----|--------------|-------|--------|------|---------|
| **A1** | Discussão, "Coerência biológica" (4 bullets dos eixos) | Listava 7 genes **fora** do CorePAM nos eixos biológicos: `FOXA1` (luminal); `UBE2C, BIRC5, TYMS, CDC20, MKI67` (proliferação); `KRT14` (basal) | Bullets reescritos contendo **apenas** os 24 genes do CorePAM; adicionado bullet "Demais reguladores selecionados" (ACTR3B, MDM2, SFRP1, GPR160, PHGDH, CXXC5, FGFR4) | Redação (contradição com a Tabela de genes e com os Resultados) | Nenhum (a análise sempre usou 24 genes) |
| **A2** | Discussão, "Comparação com assinaturas estabelecidas" | "…permanece competitivo em toda a faixa — **nunca sendo pior que os dois competidores simultaneamente em nenhuma coorte**" | "…vantagem nas coortes RNA-seq…, e **desvantagem nas coortes microarray (METABRIC e GSE20685), em que ROR-S e OncotypeDX-RS discriminam ligeiramente melhor**. Nenhuma das três é universalmente superior…" | Redação (afirmação **refutada pela própria Tabela** — CorePAM é o menor C-index dos três em METABRIC e GSE20685) | Nenhum nos dados; deixa o texto honesto e **mais defensável** |
| **A3** | Métodos (cohort count) — 3 locais | "**nove** coortes" (×3: descrição dos blocos, "Coortes do estudo", legenda do fluxograma) | "**dez** coortes (cinco de sobrevida e cinco neoadjuvantes)" | Redação (inconsistência: 5 sobrevida + 5 pCR = 10; I-SPY2 exploratória estava ora contada ora não) | Nenhum |
| **A4** | Resultados, legenda + texto da KM do SCAN-B | "HR = **3,56** (IC 95% **2,77–4,58**)"; "risco de morte **mais de três vezes** maior" | "HR = **2,83** (IC 95% **2,22–3,62**)"; "**quase três vezes** maior" | Redação alinhada à figura (a figura `Fig3_KM_SCANB` mostra 2,83; o texto dizia 3,56) | Nenhum nos resultados principais (HR contínuo do SCAN-B = 1,92, inalterado). **⚠ Confirmar no script `07_survival_analysis_SCANB.R` qual valor é o canônico — ver F1** |
| **A5** | Discussão, subgrupo ER METABRIC | "(N = **1 978**, 646 eventos DSS)" | "(N = **1 936**, 646 eventos DSS)" | Redação (a análise de subgrupo usou 1.936; 1.978 é o N total da coorte) | Nenhum |
| **A6** | Apêndice A, pseudocódigo (3 cabeçalhos) | `05_derivation.R`, `06_scoring.R`, `08_meta_analysis.R` | `05_reduce_pam50_to_corepam_FINAL.R`, `06_zscore_and_score_<coorte>.R`, `08_meta_survival.R` | Documentação (nomes não batiam com a lista oficial de scripts) | Nenhum |
| **A7** | Introdução, acesso no SUS | "o **OncotypeDX foi incorporado ao SUS** por decisão da Conitec **em 2024** apenas para subgrupo restrito" | "as assinaturas genômicas (OncotypeDX, Prosigna, MammaPrint, EndoPredict) **não estão incorporadas ao SUS** para uso rotineiro; sua adoção **permanece em avaliação** pela Conitec, **inclusive por meio de consultas públicas** sobre testagem genômica" | **Fato externo FALSO** (confirmado pelo autor e por busca: não há incorporação do OncotypeDX em 2024) | Nenhum na ciência; corrige afirmação factualmente incorreta. **Recomendado** citar o processo/consulta específico da Conitec se quiser precisão documental |
| **A8** | Resultados pCR, texto | "I² = 0 %… **ausência de heterogeneidade** — **notável** dado que…" | "I² = 0 %, **compatível com homogeneidade** — **embora K = 4 limite a precisão** dessa estimativa…" | Redação (suavização de superinterpretação; I² com K=4 é impreciso) | Nenhum |
| **A9** | Resultados pCR, legenda do forest | "(OR por 1-DP, **K = 4 coortes primárias + I-SPY2 exploratória**)… consistência **notável**" | "(OR por 1-DP, **K = 4 coortes primárias**)… consistência (ressalva K=4)" | Redação (a figura `Fig_pCR1_Forest` mostra apenas K=4; I-SPY2 não aparece nela) | Nenhum |
| **A10** | Resultados, recapitulação dos genes | "Os 24 genes **recapitulam os quatro eixos**" (listava 17 genes) | "Dos 24, **17 mapeiam** nos quatro eixos…; os **demais sete** (ACTR3B, MDM2, SFRP1, GPR160, PHGDH, CXXC5, FGFR4) contribuem com funções correlatas" | Redação (precisão: 7/24 genes ficavam fora dos 4 eixos) | Nenhum |
| **A11** | Discussão, abertura "Coerência biológica" | "representando **todos** os quatro eixos…" | "cobrindo os quatro eixos…, **além de reguladores correlatos**" | Redação (consistência com A10) | Nenhum |

### A12 — Item documentado, NÃO auto-aplicado (decisão consciente)
**Separador de milhar inconsistente:** tabelas usam ponto (`1.978`, `3.069`) e a prosa usa espaço (`1 978`, `3 069`). Ambas as convenções são aceitáveis em português; o problema é a **mistura**. **Não foi corrigido à mão** para não arriscar introduzir erros de digitação em ~30 números. Recomenda-se padronização global por regex no momento da renderização (ex.: escolher uma convenção única e aplicar com `-replace` controlado), ou configurar a formatação numérica no Quarto. **Não é erro de conteúdo.**

---

## Parte B — Itens de FIGURA pendentes (requerem regeneração da PNG ou confirmação no script — não corrigíveis só no texto)

> Estes itens foram encontrados na auditoria visual figura-a-figura, conferindo os valores exibidos contra as tabelas. A maioria exige **regerar a figura** (rodar o script R correspondente) ou **confirmar o valor canônico no script**. Onde o texto pôde ser alinhado à figura, já foi feito (ver A4, A9).

| ID | Figura | Problema | Ação recomendada | Gravidade |
|----|--------|----------|------------------|-----------|
| **F1** | `Fig3_KM_SCANB` + `Fig2_KM_MultiPanel` (Painel A) | HR dicotomizado do SCAN-B = **2,83 (2,22–3,62)** na figura, mas o texto v1.0 dizia **3,56 (2,77–4,58)**. O texto foi alinhado para 2,83 (A4). | **Confirmar no script `07_survival_analysis_SCANB.R` qual é o valor verdadeiro** (CIs quase não se sobrepõem → podem ser análises diferentes). Se 3,56 for o correto, reverter A4 e regerar a figura. | **Alta** |
| **F2** | `FigS_Cindex_ByCohort` (legenda `fig-cindex-comparison`) | A legenda promete "**CorePAM vs ROR-S vs OncotypeDX-RS**", mas a figura mostra **apenas o CorePAM, em 3 coortes** (sem SCAN-B, sem os comparadores). Legenda ≠ figura. | **Regerar a figura** com a comparação dos 3 escores (os dados estão na `tbl-comparison`), ou apontar o `include_graphics` para o arquivo correto. | **Alta** |
| **F3** | `FigS3_Correlation_OffDiagonal` (legenda `fig-correlation`) | A legenda promete "**escore intra-coorte vs z-score congelado** (r ≥ 0,958…)", mas a figura é uma **matriz de correlação ENTRE coortes** (valores 0,979–0,999). Os r congelados (0,958/0,923/0,916) não aparecem. Legenda ≠ figura. | **Regerar/substituir** pela figura correta (frozen vs intra), ou corrigir a legenda para descrever a matriz real. | **Alta** |
| **F4** | `B6_Fluxograma_Coortes` | A figura se **autocontradiz**: subtítulo diz "**9 coortes públicas**", rodapé diz "**10 coortes públicas**"; e conta o SCAN-B como coorte de validação (5 coortes, N=6.605), enquanto a `Fig1` o exclui (4 coortes, 3.536). | **Regerar B6** com contagem única (alinhar a "dez coortes" do texto; e separar derivação de validação, evitando aparência de vazamento). | Média-alta |
| **F5** | `FigS_Bootstrap_GeneFreq` | O texto diz "**15 dos 24 genes ≥ 70%**", mas a leitura visual da figura sugere **~5–6 genes ≥ 70%** (maioria entre 35–55%). | **Verificar no script `27_bootstrap_gene_stability.R`** o número real. Se o texto estiver errado, corrigir o número. (ACTR3B/NAT1 ~100% e mínimo ~35,5% batem.) | Média (verificar) |
| **F6** | `Fig_pCR2_ROC`, `FigS_DCA_pCR` (GSE25066) | GSE25066: a figura DCA mostra **pCR 20,3%** e o ROC mostra **AUC 0,569**, mas a tabela diz **23,1%** e **0,576**. Possível N/subconjunto diferente nessas figuras. | **Verificar o N/subconjunto de GSE25066** usado nas figuras de pCR (o forest `Fig_pCR1` bate com a tabela: 182/23%). Harmonizar. | Média (verificar) |
| **F7** | `Fig2_KM_MultiPanel` (Painéis A e C) | Mostram **`p = 0.0e+00`** (SCAN-B e METABRIC). **p nunca é exatamente zero.** | **Regerar** exibindo `p < 0,001` (ou `p < 2×10⁻¹⁶`). | Média |

### Itens de figura menores (cosméticos / boa prática antes da impressão)
- **Títulos truncados** (clipping de renderização) em: `FigS3_Correlation`, `FigS4_METABRIC_PCA`, `Fig3_KM_SCANB`, subtítulo de `FigS_Cindex`. → aumentar largura/margem ao regerar.
- **Inconsistência de cor:** a KM avulsa do SCAN-B usa **azul** para baixo risco; as demais usam **verde**. → padronizar.
- **p-values divergentes** entre figura avulsa e multipainel: TCGA `6,0e-04` vs `5,0e-04`; GSE20685 `1,1e-03` vs `8,5e-04`. → harmonizar.
- **`B3_Timeline`:** a narrativa "redução progressiva no número de genes" é imprecisa (70→21→50→12→24 não é monotônica); o **Breast Cancer Index está ausente** do esquema. → ajustar legenda/figura se desejado.
- **`B2_Subtipos`:** percentuais somam 90% (Normal-like 10% omitido deliberadamente) — aceitável; uma nota de rodapé resolve.
- **`FigS_DropGene`:** especificar a coorte do "máximo |ΔC| = 0,044 (EXO1)" — em GSE1456/GSE20685 o NAT1 rivaliza/supera o EXO1.

> **Observação importante:** a espinha dorsal quantitativa das figuras está **correta** — todos os N e eventos batem; o forest da meta-análise (HR 1,37; I² 38,2%), o forest de pCR (OR 1,69; I² 0%), o Delta-C incremental (todos os IC excluem zero) e todas as contagens de genes (PAM50=50, CorePAM=24, OncotypeDX=21, MammaPrint=70, EndoPredict=12) conferem. Os problemas são de **apresentação/consistência**, não de resultado.

---

## Parte C — Vulnerabilidades metodológicas (NÃO são erros — preparação para arguição)

> A revisão metodológica **não encontrou erro de método**: o desenho (congelamento analítico, folds por SHA-256, validação externa cross-plataforma, não-inferioridade pré-especificada) é rigoroso e bem executado. Os pontos abaixo são **exposições defensáveis** que a banca pode sondar — a maioria já é reconhecida na própria tese.

- **M1 — Comparador "PAM50-full (49)".** O comparador é o ponto de **maior cardinalidade da mesma trajetória elastic-net** (menos penalizado), não um PAM50 otimizado ao seu próprio λ. Logo, "CorePAM supera o PAM50" repousa, em parte, sobre um comparador sub-regularizado (que generaliza pior por construção). *Defesa:* é exatamente o argumento viés-variância da tese; mencionar que o objetivo é parcimônia transportável, não vencer o PAM50 clínico.
- **M2 — Margem de não-inferioridade na borda.** O modelo selecionado fica a **0,009** do máximo, com margem **0,010** — praticamente no limite. A afirmação "margem < ½ EP do C-index OOF" não pôde ser verificada só no texto. *Defesa:* ter à mão o EP do C-index OOF (script 05) para mostrar que 0,010 < ½ EP.
- **M3 — SCAN-B (treino) em análises rotuladas como validação.** O SCAN-B aparece na tabela de valor incremental e (na figura B6) entre as "5 coortes de validação". *Defesa:* deixar explícito que o SCAN-B é referência interna; o N de validação externa é 3.536 (4 coortes), como na Fig1.
- **M4 — Calibração da TCGA frágil.** O slope 2,74 vem de pouquíssimos pontos numa faixa estreita de probabilidade (seguimento de 24 m). *Defesa:* já discutido; reforçar que discriminação ≠ calibração e que recalibração local é recomendada.
- **M5 — z-score intra-coorte.** Depende da coorte para pontuar uma paciente isolada. *Defesa:* a análise de z-score congelado (single-sample) já cobre o uso clínico real.
- **M6 — Motivação SUS × ausência de coorte brasileira.** A tese é motivada por equidade no SUS mas não testa nenhuma coorte brasileira. *Defesa:* é prova de princípio; validação no HUB/UnB é a prioridade #1 nas perspectivas.

---

## Parte D — Declaração de integridade

As correções registradas neste documento são de natureza **editorial (redação, documentação, apresentação)** e, num único caso (A7), de **correção de um fato externo factualmente incorreto**. **Nenhuma delas altera o desenho do estudo, as análises estatísticas, os coeficientes do modelo, os escores calculados ou os resultados reportados.** As tabelas de resultados permanecem inalteradas e foram verificadas como internamente consistentes. A validade científica da tese — corroborada pela aceitação do artigo correspondente no *Breast Cancer Research* — **não é afetada**. As correções **aumentam a precisão e a defensabilidade** do documento.

---

## Parte E — Erros adicionais identificados pelo autor

*(Reservado para os erros que o autor (Rafael) trará. Serão registrados aqui com o mesmo padrão: Local · Antes · Depois · Tipo · Impacto, e marcados se já aplicados na v1.1.)*

- [ ] (a preencher)

---

## Parte F — Verificacao executando os scripts R/Python (plena certeza dos dados)

Ambiente: R 4.5.3 + pacotes (arrow, survival, survminer, glmnet, pROC, ggplot2) e Python 3 + matplotlib. Pipeline totalmente reproduzivel localmente (dados em 01_Base_Pura_CorePAM/PROCESSED/, 88 CSVs em results/). Cada item foi COMPUTADO da fonte, nao estimado.

| Item | Verificacao executada | Resultado |
|---|---|---|
| **F1** (HR SCAN-B) | Recomputado o HR dicotomizado das 5 coortes (metodo exato do script 37). | **SCAN-B = 2,83 (2,22-3,62)** confirmado. Correcao A4 (3,56->2,83) esta CERTA. (TCGA 1,78; METABRIC 2,03; GSE20685 2,13; GSE1456 2,99.) |
| **F5** (bootstrap 15/24) | Lido bootstrap_gene_stability.csv. | FALSO ALARME - exatamente 15/24 genes >=70% (ACTR3B/NAT1=100%, min ERBB2=35,5%). Texto correto. |
| **F7** (p=0.0e+00) | Underflow 1-pchisq no script 37 -> trocado por pchisq(lower.tail=FALSE). | Corrigido e multipainel REGENERADO: SCAN-B/METABRIC agora "p < 2e-16" (reais 4,05e-18 e 9,35e-19). |
| **F2** (Cindex comparison) | Script 26 so gerava CSV. Criado scripts/39_fix_figures_F2_F3.R. | Nova figura FigS_Cindex_Comparison_PT/EN (3 escores por coorte, identicos a tabela). .qmd atualizado. |
| **F3** (frozen vs intra) | Script 30 so gerava CSV. Reproduzido scoring exato por paciente. | Nova figura FigS_Frozen_vs_Intra_PT/EN (scatter+diagonal+r). r reconfirmados: 1,000/0,999/0,958/0,923/0,916. .qmd atualizado. |
| **F4** (B6) | Corrigido generate_didactic_figures.py ("9"->"10" coortes). | B6 regenerado, consistente ("10 coortes publicas"). |

### D1 - ACHADO DE INTEGRIDADE DE DADOS (decisao do PI) - GSE25066 (bloco pCR)

Ao recomputar o F6: o analysis_ready.parquet do GSE25066 NO DISCO rende N=488 (pCR 99; 20,3%; OR 1,343; AUC 0,569). Mas a TABELA/forest/meta da tese usam o GSE25066 congelado em N=182 (pCR 42; 23,1%; OR 1,435; AUC 0,576). As outras 3 coortes pCR batem exatamente.

- O prep (20_prepare_pCR_GSE25066.R) so remove pCR NA dos ~508 de Hatzis 2011 -> ~488; NAO ha filtro para 182 no codigo atual. Parquet de 01/mar; CSV congelado (182) de 03/mar.
- Consequencia: 3 paineis suplementares que recomputam ao vivo (Fig_pCR2 ROC, Fig_pCR3 quartil, Fig_pCR4 distribuicao) exibem 488; tabela/forest/meta exibem 182. INCONSISTENCIA PRE-EXISTENTE na v1.0 / bundle aceito.
- A ciencia aceita permanece integra: OR agrupado de pCR (1,69) e conclusoes usam o GSE25066 congelado (182/OR 1,435).
- PENDENTE DE DECISAO DO PI: definir o N correto do GSE25066 (quase certamente 182 do artigo aceito) e (a) restaurar/refazer o dataset de 182 e regenerar os 3 paineis, OU (b) documentar o criterio de inclusao dos 182. NAO alterar tabela/meta/forest. NAO regenerar os paineis a partir dos 488.

### Itens menores pendentes (cosmeticos)
- FigS_DCA_pCR exibe GSE25066 em 20,3% (488) - mesmo problema do D1; aguardar decisao.
- Titulos truncados; cor azul vs verde na KM avulsa do SCAN-B; B3 "reducao progressiva"/ausencia do Breast Cancer Index; B2 soma 90% - cosmeticos, opcionais.

### D1 — RESOLVIDO (2026-06-05)
Causa-raiz traceada: extract_pcr_column pegava a coluna mesclada do GEOquery 'pathologic_response_pcr_rd:ch1' (99 pCR/488). A coluna crua characteristics_ch1.10 = 42 pCR/140 RD = 182 (Hatzis 2011 = audit_pcr_definition.csv = analise aceita). FIX aplicado: script 20_prepare_pCR_GSE25066.R fixado na coluna ch1.10 com guarda de auditoria (stop se != 182/42). Parquet regenerado (182), e Fig_pCR2/3/4 + FigS_DCA_pCR regeneradas. Verificado: N=182, pCR=42, rate=23,1%, OR=1,435, AUC=0,576 — identicos a tabela/forest/meta (que NAO foram alterados). Backup do parquet 488 em analysis_ready.parquet.bak488_20260605.
