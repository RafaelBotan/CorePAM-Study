# PRD — Inspeção completa TESE_CorePAM_v0.2.pdf

**Escopo:** Verificar página a página (102 páginas) se o PDF está adequado para entrega, cruzando com a lista de modificações acordadas (tasks #1-24 + correção Sérgio/João Batista).

**Método:**
1. Render todas as páginas em PNG a 130 DPI
2. Ler em lotes de 4 (ou 2 quando a página é densa)
3. Para cada página, anotar: seção, status (OK / issue / dúvida)
4. Consolidar lista de correções
5. Aplicar correções no qmd
6. Rebuild pipeline 4-pass
7. Re-render só as páginas afetadas e confirmar fix

---

## Checklist por categoria (do backlog #1-24)

- [ ] **#10 Acentuação** — procurar palavras sem acento por toda a tese
- [ ] **#11 Termo de aprovação** — sem linhas de assinatura, placeholders em vermelho
- [ ] **#12 Lista de abreviaturas compacta**
- [ ] **#13 Lista de figuras populada** com paginação real
- [ ] **#14 Inglês italicizado** — endpoint, baseline, downstream, cross-plataforma, etc.
- [ ] **#15 Objetivos enxutos** — 1-2 gerais + até 5 específicos
- [ ] **#16 Figuras em largura total** da página
- [ ] **#17 Vancouver numérico** — `[1]`, `[2-5]`, nunca `(Perou et al., 2000)`
- [ ] **#18 Títulos azul / subtítulos azul claro**
- [ ] **#19 H3s reduzidos** — ~22 H3s esperados
- [ ] **#20 Tabelas não quebradas** — entre páginas, dentro de células
- [ ] **#21 Sem fórmulas órfãs** — LaTeX que não renderiza
- [ ] **#22 Legendas ABNT** — `Figura N — Título` com espaço antes do em-dash
- [ ] **#23 Figuras em PT** — sem labels em inglês
- [ ] **#24 Explicações didáticas** — parágrafos explicativos antes/depois das figuras
- [ ] **#3 Dedicatória + Agradecimentos** — inseridos corretamente
- [ ] **#8 4 figuras didáticas** — B2, B3, B6, B9 inseridas nos locais certos
- [ ] **Orientador correto** — João Batista de Sousa na folha de rosto, termo, agradecimentos (Sérgio Arruda só como amigo)

## Findings log

Formato: `p{N}: {seção} | {status} | {obs}`

- p1: Capa | OK
- p2: Folha de rosto | OK (orientador = João Batista de Sousa)
- p3-4: Termo de aprovação | OK (sem linhas assinatura, placeholders vermelhos)
- p5-6: Dedicatória | OK (itálico, alinhado à direita)
- p7-8: Agradecimentos | OK (João Batista como orientador, Sérgio como amigo)
- **p9: EM BRANCO | ISSUE 2** — verificar \newpage desbalanceado
- p10: Epígrafe Popper | OK
- p11: Resumo | OK
- p12: Abstract | OK
- p13-16: Lista de figuras | OK (27 figs, paginação real)
- p17: Lista de tabelas | OK (18 tabelas)
- p18: Lista de abreviaturas | OK (compacta, inline)
- **p19: Sumário | ISSUE 3** — placeholder "Clique com botão direito" em vez do sumário populado
- p20: Introdução 1.2 | OK
- p21: Fig 1 Subtipos moleculares | OK
- p22: Fig 2 Timeline assinaturas | OK
- p23: §1.5 Plataformas expressão | OK (RNA-seq, Microarray italicizado)
- p24: §1.6 Coortes + §1.7 Pergunta central | OK
- p25: Final introdução | OK
- p26: §2 Objetivos | OK (1 geral + 5 específicos)
- p27: §3 Método 3.1 | OK
- p28: Fig 3 Desenho + 3.1.2/3.1.3 | OK
- p29: §3.2 Coortes | OK
- **p30: Tabelas 1 e 2 | ISSUE 4** — quebras em cabeçalhos/células estreitas ("Endpoi\nnt", "Validaçã\no", "3 0\n69"); "I-SPY2 exp.)" com parêntese solto
- p31: Fig 4 Fluxograma | OK
- p32: Tabela 3 Parâmetros | OK (sem quebras)
- p33: §3.3.1/3.3.2 Margem e SHA-256 | OK
- **p34: Tabela 4 Cobertura | ISSUE 4 (mesmo)** — "META\nBRIC", "GSE206\n85", "GSE145\n6"
- p35: Fig 5 B9 Derivação + §3.5.1 | OK
- p36: §3.5.1 continuação | OK
- p37: Etapas 1-3 derivação | OK
- p38: §3.6 + §3.7.1/3.7.2 Cox | OK
- p39: §3.7.3 KM até §3.7.6 DCA | OK
- p40: §3.7.7 + §3.8 Meta-análise + §3.9 H2H | OK
- p41: §4 Resultados 4.1 + Fig 6 Pareto | OK
- **p42: QUASE VAZIA | ISSUE 5** — apenas "excedentes." órfão no topo
- **p43: Tabela 5 genes | ISSUE 4 + ISSUE 6** — quebras ("BLVR\nA", "ACTR\n3B", "MYBL\n2", "PTTG\n1", "MDM\n2", "GPR16\n0", "FOXC\n1", "PHGD\n H", "CXXC\n5", "KRT1\n7", "CENP\n F", "FGFR\n4") + ACENTOS FALTANDO: "Regulacao de p53", "cromossomica", "Regulacao epigenetica"
- **p44: continuação Tabela 5 | ISSUE 6** — "ERBB\n2", "Transduccao" (sem ç), "sinal HER2" (ok)
- p45: Fig 7 Lollipop pesos | OK
- p46: Recap 4 eixos biológicos | OK
- p47: Fig 8 Bootstrap estabilidade | OK
- p48: Fig 9 Leave-one-out | OK
- p49: §4.2 Validação externa + Tabela 6 | ISSUE 4 (mesma quebra: "Coort\ne", "Endpo\nint", "SCAN\n-B*", "META\nBRIC", "GSE20\n685", "GSE14\n56**")
- p50: Fig 10 KM 5 painéis | OK
- p51: Fig 11 KM SCAN-B | OK
- p52: Fig 12 KM TCGA-BRCA | OK
- p53: Fig 13 KM METABRIC | OK
- p54: Fig 14 KM GSE20685 | OK
- p55: §4.3 Meta + Tabela 7 + Fig 15 forest | OK
- **p56: QUASE VAZIA | ISSUE 5** — órfão de legenda Fig 15
- **p57: Tabela 8 Delta-C | ISSUE 4** — quebras "Covariávei\ns", "C (COR\nE-A)", "CORE-\nA +\nCorePAM"
- p58: Fig 16 Forest Delta-C | OK
- p59: §4.4 + Fig 17 Calibração | OK
- **p60: Tabelas 9, 10, 11 header | ISSUE 4** — MUITAS quebras: "Covariáv\neis", "META\nBRIC", "Popul\nação", "HR\nCorePAM\n(IC 95%)"
- p61: Tabela 11 sens. congelado + texto | ISSUE 4 (menor)
- p62: Fig 18 Correlações PCA | OK (título "Correlações... quantile" cortado à direita)
- **p63: Tabelas 12 e 13 | ISSUE 4** — "Horizont\ne" quebra
- p64: Fig 19 DCA | OK
- p65: Fig 20 PCA forense | OK (título cortado direita)
- p66: §4.7 + Tabela 14 | OK (minor quebras)
- p67: Tabela 14 + Fig 21 Forest | OK
- **p68: §4.8 Tabela 15 | ISSUE 4** — "pCR\nn (%)", "42\n(23,1\n%)"
- p69: Tabela 15 + Fig 22 pCR forest | OK
- p70: Fig 23 ROC | OK
- p71: Fig 24 Densidade pCR | OK
- p72: Fig 25 Taxa pCR quartil | OK
- **p73: §4.9 + Tabela 16 | ISSUE 4** — quebras em "congelad\no"
- p74: Tabela 16 final + Fig 26 H2H | OK
- **p75: §4.10 + Tabela 17 | ISSUE 4** — "Coort\ne", "CoreP\nAM C", "ROR-\nS C", "MET\nABRI\nC", "GSE2\n0685"
- p76: Fig 27 H2H C-index | OK
- p77: §5 Discussão 5.1 | OK
- p78: §5.2 Inovação reprodutibilidade | OK
- p79: §5.3.1 Biologia | OK
- p80: §5.4 Transportabilidade | OK
- p81: §5.5 Amostra única | OK
- p82: §5.5 cont + §5.6 | OK
- p83: §5.6 + §5.8 pCR | OK
- p84: §5.10 Reflexões metodológicas | OK
- p85: §5.11 Limitações | OK
- p86: §5.12 + §5.13 | OK
- p87: §5.14 Perspectivas | OK
- p88: §6 Conclusão | OK
- p89: Refs 1-10 (Vancouver) | OK
- p90: Refs 11-18 | OK
- p91: Refs 19-28 | OK
- p92: Refs 29-34 | OK
- p93: Refs 35-50 | OK
- p94: Refs 51-63 | OK
- p95: Refs 64-74 | OK
- p96: Refs 75-76 | OK (ref 76 Botan 2026 "Joao Batista" sem til)
- **p97: Apêndice A | ISSUE 7** — "periodico" sem acento
- **p98: Apêndice B + Tabela 18 | ISSUE 4 + 7** — "Configuracao", "Revisao major"
- **p99: §6.1 Pseudocódigo | ISSUE 8 + 9** — justify esquisito no bloco de código; numbering "6.1" em Apêndice B (deveria ser B.1)
- **p100: Algoritmo 2 | ISSUE 8** — mesmo justify esquisito
- **p101: Algoritmos 3 e 4 | ISSUE 8** — mesmo
- **p102: Final versão | ISSUE 7** — "Versao 0.2", "impressao", "Joao Batista de Sousa" (sem tilde/acento), "Pós-Graduacao"

## Consolidated issues — final status (103 pp)

1. **ISSUE 2** p9 blank — **RESOLVIDO** (removido `\newpage` redundante antes de EPÍGRAFE; Heading1 já tem page-break-before)
2. **ISSUE 3** SUMÁRIO — **RESOLVIDO** (inject_figtab_lists.py populando 74 entradas)
3. **ISSUE 4** quebras de palavra em tabelas — **NÃO RESOLVIDO** (limitação de template: colunas estreitas quebram palavras longas; exigiria retuning por tabela)
4. **ISSUE 5** páginas quase vazias — **NÃO RESOLVIDO** (fluxo natural pós-figura; aceitável)
5. **ISSUE 6** acentos Tabela 5 — **RESOLVIDO** (Regulação/cromossômica/Transdução)
6. **ISSUE 7** acentos apêndices + créditos finais + ref 76 — **RESOLVIDO**
7. **ISSUE 8** pseudocódigo com justify esquisito — **NÃO RESOLVIDO** (estilo SourceCode herda justify; precisaria modificar reference.docx)
8. **ISSUE 9** numbering B.1 — **RESOLVIDO** (explícito `## B.1 ... {.unnumbered}`)
