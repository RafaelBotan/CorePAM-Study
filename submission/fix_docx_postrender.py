"""Post-render fix: remove first of duplicated <w:pPr> in paragraphs.

When Pandoc's custom-style filter produces a caption paragraph, it injects a second
<w:pPr> containing <w:pStyle>, but the first <w:pPr> (from crossref alignment) remains.
Two pPrs in one paragraph is invalid OOXML and breaks LibreOffice rendering.
Safest fix: remove the FIRST pPr (always just jc=center), keep the SECOND (has pStyle).
"""
import zipfile, re, shutil, os

DOCX = os.environ.get('DOCX', 'TESE_CorePAM_v0.2.docx')
TMP = DOCX + '.fixtmp'

with zipfile.ZipFile(DOCX) as zin:
    entries = {n: zin.read(n) for n in zin.namelist()}

doc_xml = entries['word/document.xml'].decode('utf-8')

# Pattern: duplicate pPr pair. Remove the FIRST one entirely.
# Non-greedy match of first pPr content; then whitespace; then second pPr opening.
pattern = re.compile(r'<w:pPr>[^<]*(?:<w:[a-zA-Z]+(?:\s[^/]*)?/>[^<]*)*</w:pPr>(\s*)<w:pPr>', re.S)

# Better: specifically look for paragraph-level duplicate pPr
# Match: <w:pPr>FIRST</w:pPr>\s*<w:pPr>
dup_pattern = re.compile(r'(<w:pPr>(?:(?!</w:pPr>).)*</w:pPr>)(\s*<w:pPr>)', re.S)

def drop_first(match):
    return match.group(2).lstrip()  # drop the first pPr; keep the opening of second

# Run until no more duplicates (in case of triple)
iterations = 0
while True:
    new = dup_pattern.sub(drop_first, doc_xml)
    if new == doc_xml:
        break
    doc_xml = new
    iterations += 1
    if iterations > 5:
        break

# Verify XML validity
import xml.etree.ElementTree as ET
try:
    ET.fromstring(doc_xml)
    print('XML valid after fix.')
except ET.ParseError as e:
    print('XML STILL INVALID:', e)
    raise

# Count remaining duplicates
remaining = len(re.findall(r'<w:pPr>.*?</w:pPr>\s*<w:pPr>', doc_xml, re.S))
print(f'Iterations: {iterations}, remaining duplicates: {remaining}')

# NOTE: we keep tblLayout as fixed below (after rebalancing) so our computed
# gridCol widths are honored exactly.

# ABNT caption fix: Quarto strips leading space from title-delim, so captions render as
# "Figura 1— Texto". Insert a space before the em-dash in captions.
cap_fix_pat = re.compile(r'(Figura|Tabela)(\s|\u00a0)(\d+)(—)')
n_caps = len(cap_fix_pat.findall(doc_xml))
doc_xml = cap_fix_pat.sub(r'\1\2\3 \4', doc_xml)
print(f'Inserted space before em-dash in {n_caps} captions.')

# Patch SourceCode style: add left alignment so pseudocode blocks don't inherit
# Normal's justify (which produces huge word spacing).
styles_xml = entries['word/styles.xml'].decode('utf-8')
old_sc = ('<w:style w:type="paragraph" w:customStyle="1" w:styleId="SourceCode">'
          '<w:name w:val="Source Code" /><w:basedOn w:val="Normal" />'
          '<w:link w:val="VerbatimChar" /><w:pPr><w:wordWrap w:val="off" />'
          '<w:shd w:val="clear" w:fill="f1f3f5" /></w:pPr></w:style>')
new_sc = ('<w:style w:type="paragraph" w:customStyle="1" w:styleId="SourceCode">'
          '<w:name w:val="Source Code" /><w:basedOn w:val="Normal" />'
          '<w:link w:val="VerbatimChar" /><w:pPr><w:wordWrap w:val="off" />'
          '<w:shd w:val="clear" w:fill="f1f3f5" />'
          '<w:jc w:val="left" /></w:pPr></w:style>')
if old_sc in styles_xml:
    styles_xml = styles_xml.replace(old_sc, new_sc)
    print('Patched SourceCode style to left-align.')
else:
    print('WARN: SourceCode style pattern not found.')

# Shrink Compact (table cell) font to 10pt — ABNT minimum for table content.
# (Do NOT go below 10pt; ABNT NBR 14724 requires >=10pt for tables/captions.)
old_compact = ('<w:style w:customStyle="1" w:styleId="Compact" w:type="paragraph">'
               '<w:name w:val="Compact" /><w:basedOn w:val="BodyText" /><w:qFormat />'
               '<w:pPr><w:spacing w:after="36" w:before="36" w:line="240" '
               'w:lineRule="auto" /><w:ind w:firstLine="0" /></w:pPr></w:style>')
new_compact = ('<w:style w:customStyle="1" w:styleId="Compact" w:type="paragraph">'
               '<w:name w:val="Compact" /><w:basedOn w:val="BodyText" /><w:qFormat />'
               '<w:pPr><w:spacing w:after="36" w:before="36" w:line="240" '
               'w:lineRule="auto" /><w:ind w:firstLine="0" /></w:pPr>'
               '<w:rPr><w:sz w:val="20" /><w:szCs w:val="20" /></w:rPr></w:style>')
if old_compact in styles_xml:
    styles_xml = styles_xml.replace(old_compact, new_compact)
    print('Set Compact (table) style to 10pt (ABNT minimum).')
else:
    print('WARN: Compact style pattern not found.')

entries['word/styles.xml'] = styles_xml.encode('utf-8')

# Rebalance table column widths based on max content width per column.
# Pandoc assigns gridCol widths proportional to markdown delimiter lengths,
# producing very narrow columns (e.g. Gene=833 twips) that break words.
def rebalance_tables(xml):
    """Rebalance column widths of INNERMOST tables based on max content width.
    Pandoc wraps data tables in an outer caption-bearing table; only the inner
    one carries the actual column structure we see in the PDF.
    """
    import html as _html
    total_width_twips = 9070
    char_twips = 140  # base; bold/italic cells get weighted up
    cell_pad = 320
    min_col = 650
    # Find all <w:tbl> positions with matched </w:tbl> via depth count.
    opens = [m.start() for m in re.finditer(r'<w:tbl>', xml)]
    closes = [m.start() for m in re.finditer(r'</w:tbl>', xml)]
    # Match each open to its close using a stack
    events = sorted(
        [(p, 'o') for p in opens] + [(p, 'c') for p in closes]
    )
    stack = []
    spans = []  # list of (start, end_of_close_tag, depth, has_inner)
    inner_flag = {}  # id -> has any inner open
    for pos, kind in events:
        if kind == 'o':
            idx = len(spans)
            spans.append([pos, None, len(stack), False])
            if stack:
                spans[stack[-1]][3] = True
            stack.append(idx)
        else:
            end = pos + len('</w:tbl>')
            if stack:
                spans[stack.pop()][1] = end
    # Keep only innermost (no nested tbl inside)
    innermost = [s for s in spans if not s[3] and s[1] is not None]
    # Sort by start desc so replacements don't shift other spans
    innermost.sort(key=lambda s: s[0], reverse=True)
    new_xml = xml
    count = 0
    for s_start, s_end, _, _ in innermost:
        tbl = new_xml[s_start:s_end]
        if '<w:tblGrid>' not in tbl:
            continue
        rows = re.findall(r'<w:tr\b.*?</w:tr>', tbl, re.S)
        if not rows:
            continue
        col_max = []  # list of (desired_chars, weight_factor)
        for row in rows:
            cells = re.findall(r'<w:tc\b.*?</w:tc>', row, re.S)
            for ci, cell in enumerate(cells):
                texts = re.findall(r'<w:t[^>]*>([^<]*)</w:t>', cell)
                text = _html.unescape(''.join(texts))
                longest_word = max((len(w) for w in text.split()), default=0)
                total_len = len(text)
                desired = max(longest_word, min(total_len, 18))
                # Bold/italic runs render wider; bump width requirement
                weight = 1.0
                if '<w:b ' in cell or '<w:b/>' in cell or '<w:b/>' in cell.replace(' ', ''):
                    weight *= 1.15
                if '<w:i ' in cell or '<w:i/>' in cell or '<w:i/>' in cell.replace(' ', ''):
                    weight *= 1.08
                if ci >= len(col_max):
                    col_max.append((desired, weight))
                else:
                    d_prev, w_prev = col_max[ci]
                    col_max[ci] = (max(d_prev, desired), max(w_prev, weight))
        if not col_max:
            continue
        raw_widths = [max(min_col, int(c * char_twips * w) + cell_pad) for c, w in col_max]
        total_raw = sum(raw_widths)
        if total_raw > total_width_twips:
            scale = total_width_twips / total_raw
            widths = [max(min_col, int(w * scale)) for w in raw_widths]
        else:
            extra = total_width_twips - total_raw
            widths = [int(w + extra * (w / total_raw)) for w in raw_widths]
        diff = total_width_twips - sum(widths)
        widths[-1] += diff
        # Cohort-name columns must not become so narrow that TCGA-BRCA or
        # METABRIC split vertically. Pandoc/LibreOffice otherwise overweights
        # long numeric-result columns and makes the first column unreadable.
        first_row = rows[0] if rows else ''
        first_cells = re.findall(r'<w:tc\b.*?</w:tc>', first_row, re.S)
        first_header = ''
        if first_cells:
            first_header = _html.unescape(''.join(re.findall(r'<w:t[^>]*>([^<]*)</w:t>', first_cells[0])))
        tbl_text = _html.unescape(' '.join(
            ''.join(re.findall(r'<w:t[^>]*>([^<]*)</w:t>', cell))
            for row in rows
            for cell in re.findall(r'<w:tc\b.*?</w:tc>', row, re.S)
        ))
        if widths and 'Coorte' in first_header and any(x in tbl_text for x in ('TCGA-BRCA', 'METABRIC', 'GSE20685')):
            required_first = 1600 if len(widths) <= 7 else 1300
            if widths[0] < required_first:
                need = required_first - widths[0]
                widths[0] = required_first
                donors = sorted(range(1, len(widths)), key=lambda i: widths[i], reverse=True)
                for di in donors:
                    if need <= 0:
                        break
                    spare = max(0, widths[di] - min_col)
                    take = min(need, spare)
                    widths[di] -= take
                    need -= take
                # If all donors were already at minimum, keep total width exact
                # by subtracting the residual from the last column.
                if need > 0 and len(widths) > 1:
                    widths[-1] = max(min_col, widths[-1] - need)
                widths[-1] += total_width_twips - sum(widths)
        new_grid = '<w:tblGrid>' + ''.join(f'<w:gridCol w:w="{w}" />' for w in widths) + '</w:tblGrid>'
        tbl_new = re.sub(r'<w:tblGrid>.*?</w:tblGrid>', new_grid, tbl, count=1, flags=re.S)
        # Force fixed table layout + total width in twips so our gridCol is honored.
        tbl_new = re.sub(
            r'<w:tblLayout w:type="(?:fixed|autofit)" />',
            '<w:tblLayout w:type="fixed" />',
            tbl_new, count=1,
        )
        tbl_new = re.sub(
            r'<w:tblW w:type="pct" w:w="\d+" />',
            f'<w:tblW w:type="dxa" w:w="{total_width_twips}" />',
            tbl_new, count=1,
        )
        # Inject explicit tcW on every <w:tc> in the tbl so LO honors widths.
        tr_blocks = re.findall(r'<w:tr\b.*?</w:tr>', tbl_new, re.S)
        for tr in tr_blocks:
            cells = re.findall(r'<w:tc\b.*?</w:tc>', tr, re.S)
            new_tr = tr
            for ci, cell in enumerate(cells):
                if ci >= len(widths):
                    break
                w = widths[ci]
                # Replace <w:tcPr /> or <w:tcPr>...</w:tcPr> to include tcW.
                if '<w:tcPr />' in cell:
                    new_cell = cell.replace(
                        '<w:tcPr />',
                        f'<w:tcPr><w:tcW w:w="{w}" w:type="dxa" /></w:tcPr>',
                        1,
                    )
                elif '<w:tcPr>' in cell and '<w:tcW ' not in cell:
                    new_cell = cell.replace(
                        '<w:tcPr>',
                        f'<w:tcPr><w:tcW w:w="{w}" w:type="dxa" />',
                        1,
                    )
                else:
                    new_cell = cell
                new_tr = new_tr.replace(cell, new_cell, 1)
            tbl_new = tbl_new.replace(tr, new_tr, 1)
        new_xml = new_xml[:s_start] + tbl_new + new_xml[s_end:]
        count += 1
    return new_xml, count

doc_xml, n_rebalanced = rebalance_tables(doc_xml)
print(f'Rebalanced column widths in {n_rebalanced} innermost tables.')

# Flatten Pandoc's outer caption/figure wrappers. LibreOffice renders these as
# a stray "X" glyph at the wrapper start (bug with tblLook w:val="0000" + orphan
# bookmark anchor inside a single-cell tc).
def flatten_outer_wrappers(xml):
    opens = [m.start() for m in re.finditer(r'<w:tbl>', xml)]
    closes = [m.start() for m in re.finditer(r'</w:tbl>', xml)]
    events = sorted([(p, 'o', p + len('<w:tbl>')) for p in opens] +
                    [(p, 'c', p + len('</w:tbl>')) for p in closes])
    stack = []
    spans = []
    for s, kind, e in events:
        if kind == 'o':
            stack.append(len(spans))
            spans.append([s, None])
        else:
            if stack:
                spans[stack.pop()][1] = e
    flattened = 0
    for s, e in sorted([(s, e) for s, e in spans if e is not None],
                      key=lambda x: x[0], reverse=True):
        tbl = xml[s:e]
        lk = re.search(r'<w:tblLook[^/]*/>', tbl)
        if not lk or 'w:val="0000"' not in lk.group(0):
            continue
        # outer wrapper has single tr/tc before any nested <w:tbl>
        inner_pos = tbl.find('<w:tbl>', 10)
        outer_head = tbl if inner_pos < 0 else tbl[:inner_pos]
        if outer_head.count('<w:tc>') != 1 or outer_head.count('<w:tr>') != 1:
            continue
        # Balanced extraction of the outer tc's content
        tr_start = tbl.find('<w:tr>')
        tc_start = tbl.find('<w:tc>', tr_start) + len('<w:tc>')
        i = tc_start; depth = 1; tc_end = -1
        while i < len(tbl) and depth > 0:
            next_o = tbl.find('<w:tc>', i)
            next_c = tbl.find('</w:tc>', i)
            if next_c < 0:
                break
            if next_o >= 0 and next_o < next_c:
                depth += 1; i = next_o + 6
            else:
                depth -= 1; i = next_c + 7
                if depth == 0:
                    tc_end = next_c
                    break
        if tc_end < 0:
            continue
        inner = tbl[tc_start:tc_end]
        # strip tcPr
        inner = re.sub(r'^<w:tcPr[^/>]*/>|^<w:tcPr>[\s\S]*?</w:tcPr>', '', inner)
        xml = xml[:s] + inner + xml[e:]
        flattened += 1
    return xml, flattened

doc_xml, n_flat = flatten_outer_wrappers(doc_xml)
print(f'Flattened {n_flat} outer caption/figure wrappers.')

# Move bookmark anchors (from flatten) INSIDE the next caption/image paragraph.
# A bare <w:bookmarkStart/> between paragraphs causes LibreOffice to render a
# stray "X" marker. Relocating inside the paragraph fixes this without changing
# the anchor's logical target for hyperlinks.
# Case A: bookmark sits between paragraphs. Move it inside the following <w:p>
# AFTER its <w:pPr>...</w:pPr> block so it attaches to the first run of that
# paragraph (OOXML-valid placement that LibreOffice renders cleanly).
move_pat = re.compile(
    r'(</w:p>\s*)(<w:bookmarkStart w:id="\d+" w:name="[^"]*" />)'
    r'(<w:p>\s*<w:pPr>[\s\S]*?</w:pPr>)'
)
n_moved = 0
def mv(m):
    global n_moved
    n_moved += 1
    return m.group(1) + m.group(3) + m.group(2)
prev_doc = None
while prev_doc != doc_xml:
    prev_doc = doc_xml
    doc_xml = move_pat.sub(mv, doc_xml)
print(f'Moved {n_moved} bookmark anchors into the following paragraph.')

# Convert stand-alone page-break paragraphs into pageBreakBefore on the
# FOLLOWING paragraph — but ONLY when that paragraph is a Heading1/2/3. For
# body paragraphs, keep the explicit page break (LibreOffice honors explicit
# breaks more consistently than pageBreakBefore on non-heading runs).
# The empty page-break paragraph before a heading otherwise renders as a blank
# page (bug: the empty paragraph occupies a full page before the heading).
page_break_pat = re.compile(
    r'<w:p><w:r><w:br w:type="page"\s*/></w:r></w:p>\s*'
    r'(?:<w:bookmark(?:Start|End)[^/]*/>\s*)*'
    r'<w:p>(?:<w:pPr>(?P<ppr>[\s\S]*?)</w:pPr>)?'
)
n_pbb = 0
def convert_pb(m):
    global n_pbb
    ppr_content = m.group('ppr') or ''
    # Only convert when followed by a Heading1/2/3 paragraph (ABNT front matter).
    if not re.search(r'<w:pStyle w:val="Heading[123]"', ppr_content):
        return m.group(0)
    n_pbb += 1
    between = m.group(0)
    bkms = re.findall(r'<w:bookmark(?:Start|End)[^/]*/>', between)
    bk_xml = ''.join(bkms)
    if '<w:pageBreakBefore' in ppr_content:
        new_ppr = ppr_content
    else:
        new_ppr = '<w:pageBreakBefore/>' + ppr_content
    return bk_xml + '<w:p><w:pPr>' + new_ppr + '</w:pPr>'
doc_xml = page_break_pat.sub(convert_pb, doc_xml)
print(f'Converted {n_pbb} page-break paragraphs to pageBreakBefore on next Heading1/2/3.')

entries['word/document.xml'] = doc_xml.encode('utf-8')

with zipfile.ZipFile(TMP, 'w', zipfile.ZIP_DEFLATED) as zout:
    for name, data in entries.items():
        zout.writestr(name, data)

shutil.move(TMP, DOCX)
print('docx repacked:', os.path.getsize(DOCX))
