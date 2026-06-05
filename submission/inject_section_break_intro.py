"""Inject a section break BEFORE the INTRODUÇÃO heading so the front matter uses
a sectPr without headerReference (no page numbers), while the body inherits the
final sectPr with the PAGE-field header. Pagination continues (does not restart),
so INTRODUÇÃO simply shows its absolute page number.

Runs after inject_figtab_lists.py (which also edits document.xml).
"""
import zipfile, re, shutil, os

DOCX = os.environ.get('DOCX', 'TESE_CorePAM_v0.2.docx')
TMP = DOCX + '.secbreak'

with zipfile.ZipFile(DOCX) as z:
    entries = {n: z.read(n) for n in z.namelist()}

doc = entries['word/document.xml'].decode('utf-8')

# Grab pgSz/pgMar from final sectPr to mirror in the front-matter sectPr
final_sec = re.search(r'<w:sectPr\b[^>]*>[\s\S]*?</w:sectPr>', doc)
if not final_sec:
    raise SystemExit('Final sectPr not found.')
pgSz = re.search(r'<w:pgSz\b[^/]*/>', final_sec.group(0))
pgMar = re.search(r'<w:pgMar\b[^/]*/>', final_sec.group(0))
pgSz_xml = pgSz.group(0) if pgSz else '<w:pgSz w:w="11906" w:h="16838"/>'
pgMar_xml = pgMar.group(0) if pgMar else '<w:pgMar w:top="1701" w:right="1134" w:bottom="1134" w:left="1701" w:header="708" w:footer="708" w:gutter="0"/>'

# Front-matter sectPr: NO headerReference → no page numbers shown
front_sectPr_inner = pgSz_xml + pgMar_xml + '<w:cols w:space="708"/><w:docGrid w:linePitch="360"/>'
break_paragraph = (
    '<w:p><w:pPr><w:sectPr>' + front_sectPr_inner + '</w:sectPr></w:pPr></w:p>'
)

# Locate paragraph containing the INTRODUÇÃO heading text
target = re.search(
    r'<w:p[^>]*>(?:(?!</w:p>).)*?w:pStyle w:val="Heading1"(?:(?!</w:p>).)*?INTRODU[ÇC][ÃA]O(?:(?!</w:p>).)*?</w:p>',
    doc, flags=re.S,
)
if not target:
    raise SystemExit('INTRODUÇÃO heading paragraph not found.')

p_start = target.start()
# Insert the section-break paragraph right BEFORE the INTRODUÇÃO paragraph.
doc = doc[:p_start] + break_paragraph + doc[p_start:]
print('Injected front-matter section break before INTRODUÇÃO.')

# Validate
import xml.etree.ElementTree as ET
ET.fromstring(doc)
print('XML valid after section-break injection.')

entries['word/document.xml'] = doc.encode('utf-8')
with zipfile.ZipFile(TMP, 'w', zipfile.ZIP_DEFLATED) as zout:
    for n, d in entries.items():
        zout.writestr(n, d)
shutil.move(TMP, DOCX)
print(f'Section break injected. docx size: {os.path.getsize(DOCX)}')
