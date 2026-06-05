# HANDOFF — Tese CorePAM v1.1 (continuação em nova sessão)

**Data:** 2026-06-05 · **Para:** próxima sessão Claude Code · **Projeto:** `Y:\Doutorado Botan\Estudos\CorePAM_Study_accepted`

Defesa: **22/06/2026 14h, Auditório FM-UnB**. A tese-fonte de trabalho é o `.qmd` v1.1; o original v1.0 está **intacto**.

---

## 0) CONTEXTO RÁPIDO (o que aconteceu nesta sessão)

Revisão pós-defesa exaustiva da tese CorePAM. Foram corrigidos ~11 erros de texto, ~7 de figura, **1 de integridade de dados (GSE25066)**, consolidadas interpretações duplicadas de figuras, renderizado o PDF e feita 1ª passada de revisão. Tudo commitado **localmente** (commit `2d1d687`), mas **NÃO pushado** (auth pendente).

**Arquivos canônicos:**
- Fonte corrigida: `submission/TESE_CorePAM_v1.1_corrigida.qmd`
- Render: `submission/TESE_CorePAM_v1.1_corrigida.docx` e `.pdf` (5,47 MB)
- Original preservado: `submission/TESE_CorePAM_FINAL_envio_Dr_Joao_impressao.qmd` (v1.0)
- Registro de todas as correções: `submission/TESE_CorePAM_REGISTRO_CORRECOES_v1.1.md` (Partes A–F)

---

## 1) JÁ FEITO E VERIFICADO ✅

- **Texto A1–A11** (ver registro Parte A): genes da Discussão só os 24; comparação ROR-S/ODX honesta (CorePAM perde em microarray); "dez coortes"; HR SCAN-B = **2,83** (computado da fonte); N subgrupo METABRIC 1.936; nomes de script; **afirmação falsa de incorporação no SUS REMOVIDA**; T14/T15 citadas; `x`→`×`; título de calibração; fechamentos 4.7/4.8; tabela de subgrupo unificada (T16); banca sem justificação espalhada.
- **Figuras F1–F7** (Parte F): `p=0.0e+00` do multipainel KM corrigido (`lower.tail=FALSE` no `scripts/37_R2_km_multipanel_v12.R`); **figuras criadas** `FigS_Cindex_Comparison` e `FigS_Frozen_vs_Intra` (via `scripts/39_fix_figures_F2_F3.R`); B6 "9→10 coortes" (`generate_didactic_figures.py`).
- **DADO GSE25066 (D1) — CRÍTICO:** restaurado a **N=182** (era 488 por drift do auto-detect do GEOquery, que pegava a coluna mesclada `pathologic_response_pcr_rd:ch1`). O prep `scripts/20_prepare_pCR_GSE25066.R` agora está **fixado** na coluna crua `characteristics_ch1.10` (42 pCR/140 RD) **com guarda de auditoria** (`stop` se ≠ 182/42). Figuras pCR regeneradas (ROC 0,576, taxa 23,1%, OR 1,435 — batem com tabela/forest/meta, que **NÃO foram alterados**). Backup do parquet 488: `01_Base_Pura_CorePAM/PROCESSED/pCR/GSE25066/analysis_ready.parquet.bak488_20260605`.
- **#3 interpretações duplicadas:** 12 pares genuinamente redundantes consolidados (4 KMs mescladas + pareto, multipainel, calibração, correlação, er-forest, pcr-forest, head-to-head, cindex), preservando o didático complementar.
- **Render + 1ª passada (conteúdo):** PDF gerado; verificada fidelidade total da pré-textual + introdução (págs. 1–22): cross-refs resolvem ("Figura N"), citações viram superscritos (vancouver), figuras e números corretos. Folha de aprovação corrigida (pág. 3 OK).
- **Git:** commit local `2d1d687` (874 arquivos; dados grandes em `01_Base_Pura_CorePAM/` estão gitignored).

---

## 2) O QUE FALTA (com instruções)

### 2.1 — PUSH do git (auth pendente)
O commit `2d1d687` está só local. Push falhou: "Password authentication is not supported".
- Remoto: `https://github.com/RafaelBotan/CorePAM-Study.git`, branch `main`.
- git só existe no **WSL** (git 2.43.0); não há git.exe no Windows. `gh.exe` (Windows) está autenticado como RafaelBotan (scope repo).
- **Tentado e falhou:** passar token do gh via WSLENV+credential.helper — o `$GHT` provavelmente não atravessou para o ambiente WSL.
- **Caminhos a tentar na próxima sessão:**
  1. Pedir ao usuário rodar `! ` (no prompt do Claude Code) um `git push` autenticado, OU
  2. No WSL, armazenar credencial: `wsl -e bash -lc 'git config --global credential.helper store'` e fazer um push interativo onde o usuário cola um PAT (usuário=RafaelBotan, senha=PAT com scope repo), OU
  3. Gerar PAT e exportar dentro do WSL: `wsl -e bash -lc 'export GH=<PAT>; git -C "/mnt/y/Doutorado Botan/Estudos/CorePAM_Study_accepted" push https://RafaelBotan:$GH@github.com/RafaelBotan/CorePAM-Study.git main'` (cuidado: não imprimir o PAT).
- **Regra:** não imprimir tokens no chat. Senhas/tokens do usuário podem estar em `C:\Users\oncol\Desktop\Senhas.docx`.

### 2.2 — Sumário, Lista de Figuras, Lista de Tabelas (estão como placeholder)
No `.qmd` há 3 seções placeholder ("gerar após a paginação final"): `SUMÁRIO`, `LISTA DE FIGURAS`, `LISTA DE TABELAS`. O render docx **não** as gera. Poppler/pdftotext **não está instalado** (não dá para extrair paginação facilmente).
- **Opção recomendada (mais simples e robusta):** abrir o `.docx` no **Word** → Referências → Inserir Sumário; e Inserir Índice de Ilustrações (uma vez para Figuras, filtro "Figura"; outra para Tabelas, filtro "Tabela"). O Word numera as páginas e atualiza com F9. Substituir os 3 placeholders.
- **Opção alternativa (automatizável):** macro LibreOffice Basic headless que insere `com.sun.star.text.ContentIndex` (ToC), `IllustrationIndex` (figuras) e `TableIndex` (tabelas), atualiza e exporta PDF. Mais trabalhoso; só se quiser 100% reprodutível.
- **Opção Quarto:** setar `toc: true` no YAML gera um Sumário (mas a posição no docx fica no topo, fora do padrão ABNT; LoF/LoT não saem). Útil só para o Sumário.

### 2.3 — Revisão página-a-página COMPLETA (1ª passada conteúdo + 2ª forma)
Feito só págs. 1–22 (pré-textual+intro). **Falta ler Métodos, Resultados, Discussão, Conclusão, Referências e Apêndices** página a página. O motor de render é consistente (o que renderizou nas 22 págs renderiza no resto), e todas as correções estão na fonte — então o risco é baixo, mas o usuário pediu explicitamente a leitura completa.
- **1ª passada (conteúdo):** conferir números das tabelas/figuras de Resultados (HR 2,83 SCAN-B, meta 1,37, pCR 0,576/23,1%, head-to-head, comparação genefu), cross-refs resolvidos (sem "??"), Referências renderizadas (bib `tese_references_FINAL_envio_Dr_Joao.bib`, sem citações quebradas), Apêndice A (pipeline) e Apêndice B (artigo BCR, se incluído).
- **2ª passada (forma/layout):** paginação, tabelas que não estouram a margem (especialmente `tbl-head-to-head` 8 colunas, `tbl-comparison`, a nova `tbl-er-subgroup` de 7 colunas), figuras centradas e legíveis, legendas completas, quebras de página, conformidade ABNT.
- **Forma já detectada:** instituições longas da banca (ex.: Hélio Carrara) ainda justificam a 1ª linha quando quebram (menor); B3 timeline tem subtítulo "redução progressiva" que é impreciso (opcional). 

---

## 3) TOOLCHAIN DE REFERÊNCIA (tudo testado nesta sessão)

```
R:        "C:\Program Files\R\R-4.5.3\bin\Rscript.exe"   (R 4.5.2 está QUEBRADO; usar 4.5.3)
Quarto:   "C:\Program Files\Quarto\bin\quarto.exe"        (1.8.27)
soffice:  "C:\Program Files\LibreOffice\program\soffice.exe"
git:      via WSL  ->  wsl -e git -C "/mnt/y/Doutorado Botan/Estudos/CorePAM_Study_accepted" ...
          (se WSL der "Catastrophic failure": wsl --shutdown ; aguardar 5s ; retentar)
```

**Render (rodar da pasta `submission`, para os caminhos relativos `../figures/...` resolverem):**
```powershell
$env:COREPAM_ROOT="Y:\Doutorado Botan\Estudos\CorePAM_Study_accepted"; Set-Location "$env:COREPAM_ROOT\submission"
& "C:\Program Files\Quarto\bin\quarto.exe" render "TESE_CorePAM_v1.1_corrigida.qmd" --to docx
# docx -> pdf (IMPORTANTE: fechar instâncias do LibreOffice e usar perfil isolado, senão falha silenciosa):
Get-Process soffice* -ea SilentlyContinue | Stop-Process -Force
& "C:\Program Files\LibreOffice\program\soffice.exe" --headless "-env:UserInstallation=file:///C:/Users/oncol/AppData/Local/Temp/lo_corepam_profile" --convert-to pdf --outdir "$env:COREPAM_ROOT\submission" "TESE_CorePAM_v1.1_corrigida.docx"
```

**Rodar scripts R do pipeline (da raiz do repo):**
```powershell
$env:COREPAM_ROOT="Y:\Doutorado Botan\Estudos\CorePAM_Study_accepted"; Set-Location $env:COREPAM_ROOT; $env:FORCE_RERUN="TRUE"
& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" "scripts\<nome>.R"
```

---

## 4) AVISOS CRÍTICOS ⚠️

1. **GSE25066 = 182, NÃO 488.** Se alguém rodar o prep sem a correção, ou se um sync restaurar o parquet 488, a guarda de auditoria no `20_prepare_pCR_GSE25066.R` vai dar `stop`. Os 182 são os corretos (Hatzis 2011, `audit_pcr_definition.csv`). **Nunca** regenerar figuras pCR a partir dos 488. **Nunca** alterar tabela/forest/meta de pCR (já corretos com 182).
2. **Não tocar no v1.0** (`..._envio_Dr_Joao_impressao.qmd`) — é o aceito/preservado.
3. Se re-rodar `scripts/39_fix_figures_F2_F3.R`, ele recria `FigS_Cindex_Comparison` e `FigS_Frozen_vs_Intra` (o `.qmd` já aponta para elas).
4. O sync (rclone/git) da `.claude` já causou corrupção antes (plugins). Cuidado com sync sobre `01_Base_Pura_CorePAM`.

---

## 5) PRÓXIMOS PASSOS SUGERIDOS (ordem)
1. Resolver o **push** (2.1).
2. Gerar **Sumário/LoF/LoT** no Word (2.2) e re-exportar PDF.
3. **Revisão página-a-página** completa (2.3), corrigindo o que aparecer no `.qmd` e re-renderizando.
4. Declarar a tese pronta.
5. Rodar **/checkknowledge** para consolidar memórias (caso GSE25066, pipeline de render, bugs de figura).
