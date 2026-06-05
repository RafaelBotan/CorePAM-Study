"""Populate LISTA DE FIGURAS and LISTA DE TABELAS statically from captions + PDF pagination.

Runs after fix_docx_postrender.py and the first LibreOffice PDF pass.
Reads captions from document.xml, matches each to its page in the PDF, then
replaces the placeholder paragraphs with a populated list (one paragraph per entry
with caption text left-aligned and page number right-aligned via tab stop).
"""
import zipfile, re, shutil, os, pypdf, html

DOCX = os.environ.get('DOCX', 'TESE_CorePAM_v0.2.docx')
PDF = os.environ.get('PDF', 'TESE_CorePAM_v0.2.pdf')
TMP = DOCX + '.listtmp'

with zipfile.ZipFile(DOCX) as z:
    entries = {n: z.read(n) for n in z.namelist()}
doc = entries['word/document.xml'].decode('utf-8')

# Extract captions by iterating <w:p> blocks
caps = []
pos = 0
while True:
    s1 = doc.find('<w:p>', pos); s2 = doc.find('<w:p ', pos)
    cands = [s for s in (s1, s2) if s >= 0]
    if not cands: break
    start = min(cands)
    end = doc.find('</w:p>', start)
    if end < 0: break
    block = doc[start:end+len('</w:p>')]
    if 'w:pStyle w:val="ImageCaption"' in block:
        txt = ''.join(re.findall(r'<w:t[^>]*>([^<]*)</w:t>', block))
        txt = html.unescape(txt).strip()
        caps.append(txt)
    pos = end + 1

def short_title(full):
    """ABNT: lista de figuras/tabelas deve mostrar só o título, não a legenda inteira.
    Corta no primeiro ponto final seguido de espaço (fim da frase-título)."""
    full = full.strip()
    # Procura o primeiro ". " — fim da frase-título; legenda explicativa vem depois.
    m = re.search(r'\.\s+[A-ZÁÉÍÓÚÂÊÔÃÕÇ]', full)
    if m:
        title = full[:m.start()].rstrip('.').strip()
    # Se não houver, tenta cortar no primeiro ponto final simples, caso exista.
    elif full.endswith('.'):
        title = full[:-1].strip()
    else:
        title = full
    # Listas longas no Word podem colar o número da página ao título quando o
    # texto chega no tab stop. Para a lista, removemos detalhes parentéticos de
    # títulos extensos; a legenda completa permanece no corpo da tese.
    if len(title) > 58:
        title = re.sub(r'\s*\([^)]{8,}\)', '', title)
    title = title.replace('das principais assinaturas gênicas', 'das assinaturas gênicas')
    title = title.replace(
        'Correlação entre escores intra-coorte e z-score congelado por coorte',
        'Correlação entre escores por coorte',
    )
    title = re.sub(r'\s+', ' ', title).strip()
    return title

figs = []  # (num, title)
tabs = []
for t in caps:
    mf = re.match(r'^Figura\s*(\d+)\s*[:\.\—\-–]\s*(.+)$', t)
    mt = re.match(r'^Tabela\s*(\d+)\s*[:\.\—\-–]\s*(.+)$', t)
    if mf: figs.append((int(mf.group(1)), short_title(mf.group(2))))
    elif mt: tabs.append((int(mt.group(1)), short_title(mt.group(2))))
figs.sort(); tabs.sort()
print(f'Captions parsed: {len(figs)} figuras, {len(tabs)} tabelas')

# Build page map from PDF
r = pypdf.PdfReader(PDF)
page_texts = []
for p in r.pages:
    try: page_texts.append(p.extract_text() or '')
    except: page_texts.append('')

def _body_start_hint():
    """Locate the first body page — i.e. the page after the SUMÁRIO/TOC block.
    Pre-textual pages (LISTA DE FIGURAS/TABELAS/SUMÁRIO) duplicate caption and
    heading strings and must be skipped or we always match them first.
    """
    # Prefer the real INTRODUÇÃO heading. In the sumário it appears on a line
    # with dotted leaders; in the body it appears as a standalone heading line.
    intro_pat = re.compile(r'^\s*1\.\s+INTRODU[ÇC][ÃA]O\s*$', re.I)
    for i, pt in enumerate(page_texts):
        pt_norm = pt.replace('\u00a0', ' ')
        for line in pt_norm.splitlines():
            if intro_pat.match(line.strip()):
                return i

    sumario_start = None
    for i, pt in enumerate(page_texts):
        if 'SUMÁRIO' in pt or 'SUMARIO' in pt:
            sumario_start = i
            break
    if sumario_start is None:
        return 0
    body_start = sumario_start + 1
    for j in range(sumario_start, len(page_texts)):
        leaders = len(re.findall(r'\.{3,}', page_texts[j]))
        if leaders >= 5:
            body_start = j + 1
        else:
            break
    return body_start

def find_page(kind, num):
    # match any of: ':' '—' '–' '-' after the number, with optional space or nbsp
    pat = re.compile(rf'{kind}[ \u00a0]{num}[ \u00a0]*[:\—\–\-]', re.U)
    # Skip pre-textual pages (LISTA DE FIGURAS/TABELAS/SUMÁRIO) which contain
    # identical caption strings. Fallback to global scan if not found.
    start = _body_start_hint()
    for i in range(start, len(page_texts)):
        if pat.search(page_texts[i]):
            return i + 1
    for i, pt in enumerate(page_texts):
        if pat.search(pt):
            return i + 1
    return None

# Preliminary page, which may not exist in TOC pages (pages 10-11 in current output) — skip them
skip_pages = set()  # we'll rely on find_page which picks first occurrence

fig_entries = [(n, t, find_page('Figura', n)) for n, t in figs]
tab_entries = [(n, t, find_page('Tabela', n)) for n, t in tabs]

# If a caption occurs both in its figure page AND in the list-of-figures (circular), first match wins.
# Since lists currently have placeholder text (no captions), first real match is the real one. OK.

for n, t, pg in fig_entries[:5]: print(f'  Figura {n} -> p{pg}')
for n, t, pg in tab_entries[:5]: print(f'  Tabela {n} -> p{pg}')

# Build OOXML list paragraphs
NS = 'xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"'

def esc(s):
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')

def build_entry(kind, num, title, page):
    page_str = str(page) if page else '—'
    line = f'{kind} {num} – {title}'
    line_esc = esc(line)
    # right tab at ~9000 twips (page width)
    return (
        '<w:p>'
        '<w:pPr>'
        '<w:tabs><w:tab w:val="right" w:leader="dot" w:pos="8640"/></w:tabs>'
        '<w:spacing w:after="0" w:line="360" w:lineRule="auto"/>'
        '<w:ind w:left="0" w:right="0"/>'
        '</w:pPr>'
        f'<w:r><w:t xml:space="preserve">{line_esc}</w:t></w:r>'
        '<w:r><w:tab/>'
        f'<w:t xml:space="preserve">{page_str}</w:t></w:r>'
        '</w:p>'
    )

# Wrap each injected list in XML comment markers so this script can be re-run
# on an already-populated docx (second pass, to refresh page numbers after the
# first injection changed the pagination).
def wrap(name, xml):
    return f'<!--{name}_START-->{xml}<!--{name}_END-->'

fig_xml = wrap('FIGLIST', ''.join(build_entry('Figura', n, t, pg) for n, t, pg in fig_entries))
tab_xml = wrap('TABLIST', ''.join(build_entry('Tabela', n, t, pg) for n, t, pg in tab_entries))

# Find and replace placeholder paragraphs.
# Placeholder text fragment to locate:
FIG_PH = 'gerar a lista de figuras'
TAB_PH = 'gerar a lista de tabelas'

def replace_placeholder(doc_xml, ph_text, inject_xml, marker_name):
    """Replace placeholder paragraph on first pass, or re-replace wrapped content on subsequent passes."""
    start_tag = f'<!--{marker_name}_START-->'
    end_tag = f'<!--{marker_name}_END-->'
    s = doc_xml.find(start_tag)
    if s >= 0:
        e = doc_xml.find(end_tag, s)
        if e >= 0:
            return doc_xml[:s] + inject_xml + doc_xml[e + len(end_tag):]
    idx = doc_xml.find(ph_text)
    if idx < 0:
        print(f'Placeholder NOT FOUND and no marker: "{ph_text}"')
        return doc_xml
    p_start = doc_xml.rfind('<w:p>', 0, idx)
    p_start2 = doc_xml.rfind('<w:p ', 0, idx)
    if p_start2 > p_start: p_start = p_start2
    p_end = doc_xml.find('</w:p>', idx) + len('</w:p>')
    return doc_xml[:p_start] + inject_xml + doc_xml[p_end:]

doc = replace_placeholder(doc, FIG_PH, fig_xml, 'FIGLIST')
doc = replace_placeholder(doc, TAB_PH, tab_xml, 'TABLIST')

# ---- SUMÁRIO (TOC from Heading1/2/3) ----
heads = []  # list of (level, text)
pos = 0
while True:
    s1 = doc.find('<w:p>', pos); s2 = doc.find('<w:p ', pos)
    cands = [s for s in (s1, s2) if s >= 0]
    if not cands: break
    start = min(cands)
    end = doc.find('</w:p>', start)
    if end < 0: break
    block = doc[start:end+len('</w:p>')]
    m = re.search(r'w:pStyle w:val="Heading(\d)"', block)
    if m:
        lvl = int(m.group(1))
        if lvl in (1, 2, 3):
            txt = ''.join(re.findall(r'<w:t[^>]*>([^<]*)</w:t>', block))
            txt = html.unescape(txt).strip()
            if txt: heads.append((lvl, txt))
    pos = end + 1

# Skip pre-textual entries that ABNT keeps out of Sumário
# (Sumário itself must not self-reference; we also drop elements that appear
#  before it in the document, but keep from RESUMO onwards.)
pretextual_drop = {
    'RESUMO', 'ABSTRACT',
    'LISTA DE FIGURAS', 'LISTA DE TABELAS', 'LISTA DE ABREVIATURAS E SIGLAS',
    'SUMÁRIO',
}
# Find index of "SUMÁRIO" in heads and start from next
toc_start = 0
for i, (lv, t) in enumerate(heads):
    if lv == 1 and t.upper().startswith('SUMÁRIO'):
        toc_start = i + 1
        break
toc_heads = heads[toc_start:]
# Also drop any remaining pretextual-style entry
toc_heads = [(lv, t) for lv, t in toc_heads if t.upper() not in pretextual_drop]

_body_start = _body_start_hint()
print(f'Body starts at page {_body_start + 1} (0-indexed {_body_start})')

def find_heading_page(text):
    # Search for heading text in PDF pages, SKIPPING TOC pages (which would
    # otherwise match "1. INTRODUÇÃO" etc. before the real body heading).
    needle = text[:40].strip().replace('\u00a0', ' ')
    needle_esc = re.escape(needle)
    pat = re.compile(needle_esc.replace(r'\ ', r'\s+'))
    for i in range(_body_start, len(page_texts)):
        pt_norm = page_texts[i].replace('\u00a0', ' ')
        if pat.search(pt_norm):
            return i + 1
    for i, pt in enumerate(page_texts):
        pt_norm = pt.replace('\u00a0', ' ')
        if pat.search(pt_norm):
            return i + 1
    return None

toc_entries = [(lv, t, find_heading_page(t)) for lv, t in toc_heads]

def build_toc_entry(level, title, page):
    page_str = str(page) if page else '—'
    line_esc = esc(title)
    indent_left = {1: 0, 2: 360, 3: 720}.get(level, 0)
    return (
        '<w:p>'
        '<w:pPr>'
        '<w:tabs><w:tab w:val="right" w:leader="dot" w:pos="8640"/></w:tabs>'
        '<w:spacing w:after="0" w:line="360" w:lineRule="auto"/>'
        f'<w:ind w:left="{indent_left}" w:right="0"/>'
        '</w:pPr>'
        f'<w:r><w:t xml:space="preserve">{line_esc}</w:t></w:r>'
        '<w:r><w:tab/>'
        f'<w:t xml:space="preserve">{page_str}</w:t></w:r>'
        '</w:p>'
    )

toc_xml = wrap('TOC', ''.join(build_toc_entry(lv, t, pg) for lv, t, pg in toc_entries))
print(f'SUMÁRIO entries: {len(toc_entries)}')
for lv, t, pg in toc_entries[:8]: print(f'  [H{lv}] {t[:50]} -> p{pg}')

TOC_PH = 'gerar o sumário'
doc = replace_placeholder(doc, TOC_PH, toc_xml, 'TOC')

# Verify XML
import xml.etree.ElementTree as ET
try:
    ET.fromstring(doc)
    print('XML valid after list injection.')
except ET.ParseError as e:
    print('XML INVALID:', e)
    raise

entries['word/document.xml'] = doc.encode('utf-8')
with zipfile.ZipFile(TMP, 'w', zipfile.ZIP_DEFLATED) as zout:
    for n, d in entries.items():
        zout.writestr(n, d)
shutil.move(TMP, DOCX)
print(f'Lists injected. docx size: {os.path.getsize(DOCX)}')
