# -*- coding: utf-8 -*-
# Constroi Sumario, Lista de Figuras, Lista de Tabelas.
# Pagina: ULTIMA ocorrencia no PDF (a legenda/heading REAL vem depois das listas).
# Legenda: vinda do .qmd (fig-cap / caption), completa e truncada limpa.
import sys, re
import pypdf

pdf_path, qmd_path = sys.argv[1], sys.argv[2]

# ---- 1) Captions do QMD, em ordem (Nth = numero N) ----
with open(qmd_path, encoding='utf-8') as f:
    qtext = f.read()

def strip_md(s):
    s = re.sub(r'[*_`]', '', s)
    s = s.replace('\\$', '$').replace('\\%', '%').replace('\\&', '&')
    s = re.sub(r'\s+', ' ', s).strip()
    return s

def short_cap(s, n=78):
    s = strip_md(s)
    # primeira sentenca se curta o suficiente
    cut = s.find('. ')
    if 30 <= cut <= n:
        return s[:cut]
    if len(s) > n:
        sp = s.rfind(' ', 0, n)
        s = (s[:sp] if sp > 30 else s[:n]).rstrip(' ,;:').rstrip()
        # remove preposicao/artigo solto no fim
        s = re.sub(r'\s+(do|da|de|e|em|na|no|com|para|por|que|os|as|um|uma)$', '', s, flags=re.I)
        s += '…'
    return s

fig_caps = re.findall(r'#\|\s*fig-cap:\s*"((?:[^"\\]|\\.)*)"', qtext)
# tabelas: linhas ": Texto {#tbl-...}"  (ignora as auto-listas, que nao tem {#tbl-})
tbl_caps = re.findall(r'^:\s*(.+?)\s*\{#tbl-[^}]*\}', qtext, flags=re.M)
fig_caps = [short_cap(c) for c in fig_caps]
tbl_caps = [short_cap(c) for c in tbl_caps]

# ---- 2) Paginas (ULTIMA ocorrencia) no PDF ----
r = pypdf.PdfReader(pdf_path)
CHAP = {"INTRODUÇÃO","OBJETIVOS","MÉTODO","RESULTADOS","DISCUSSÃO","CONCLUSÃO"}
POSTEXT = re.compile(r'^(REFER[ÊE]NCIAS BIBLIOGR[ÁA]FICAS|AP[ÊE]NDICE [A-Z])')
H_NUM = re.compile(r'^(\d+(?:\.\d+)*)\.?\s+(\S.{1,90})$')
FIG = re.compile(r'^Figura\s+(\d+)\b')
TBL = re.compile(r'^Tabela\s+(\d+)\b')

fig_pg={}; tbl_pg={}; head={}   # head[num]=(lvl,title,page)
for i,pg in enumerate(r.pages,1):
    for raw in (pg.extract_text() or "").split("\n"):
        line=raw.strip()
        if not line: continue
        m=FIG.match(line)
        if m: fig_pg[int(m.group(1))]=i; continue
        m=TBL.match(line)
        if m: tbl_pg[int(m.group(1))]=i; continue
        if POSTEXT.match(line):
            head[line[:22]]=(0, re.sub(r'\s+\d+$','',line).strip(), i); continue
        m=H_NUM.match(line)
        if m:
            num=m.group(1); title=re.sub(r'\s+\d+$','',m.group(2)).strip()
            chap=num.split('.')[0]
            if not chap.isdigit() or not (1<=int(chap)<=6): continue
            lvl=len(num.split('.'))
            if lvl==1:
                if title.split()[0].upper().rstrip(':') not in CHAP: continue
            else:
                if not title[0].isalpha(): continue
            head[num]=(lvl, f"{num} {title}", i)

# ordena headings: numericos por tupla, postextuais ao fim por pagina
def hkey(k):
    v=head[k]
    if v[0]==0: return (99,)+(v[2],)
    return tuple(int(x) for x in k.split('.'))
toc=[head[k] for k in sorted(head, key=hkey)]

def table(col1, rows):
    out=[f"| {col1} | Página |","|:---|---:|"]
    for label,page,bold in rows:
        out.append(f"| {'**'+label+'**' if bold else label} | {page} |")
    return "\n".join(out)

sumario = table("Seção", [(t,p,(lvl in (0,1))) for (lvl,t,p) in toc])
lof = table("Figura", [(f"Figura {n} — {fig_caps[n-1]}", fig_pg.get(n,'?'), False)
                        for n in range(1,len(fig_caps)+1)])
lot = table("Tabela", [(f"Tabela {n} — {tbl_caps[n-1]}", tbl_pg.get(n,'?'), False)
                        for n in range(1,len(tbl_caps)+1)])

def inject(q, marker, placeholder, content):
    block=f"<!-- {marker}_START -->\n\n{content}\n\n<!-- {marker}_END -->"
    pat=re.compile(re.escape(f"<!-- {marker}_START -->")+r".*?"+re.escape(f"<!-- {marker}_END -->"),re.S)
    return pat.sub(lambda m: block, q) if pat.search(q) else q.replace(placeholder, block, 1)

q=qtext
q=inject(q,"AUTO_SUMARIO","Placeholder técnico: gerar o sumário após a paginação final.",sumario)
q=inject(q,"AUTO_LOF","Placeholder técnico: gerar a lista de figuras após a paginação final.",lof)
q=inject(q,"AUTO_LOT","Placeholder técnico: gerar a lista de tabelas após a paginação final.",lot)
with open(qmd_path,'w',encoding='utf-8') as f: f.write(q)

print(f"figcaps={len(fig_caps)} tblcaps={len(tbl_caps)} | fig_pg={len(fig_pg)} tbl_pg={len(tbl_pg)} | toc={len(toc)}")
print("amostra LoF:", [f"F{n}:{fig_pg.get(n,'?')}" for n in range(1,8)])
print("amostra LoT:", [f"T{n}:{tbl_pg.get(n,'?')}" for n in range(1,8)])
print("INTRO/RESULT/DISC pg:", [(t,p) for (lvl,t,p) in toc if lvl==1][:6])
