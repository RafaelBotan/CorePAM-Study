"""Patch reference_abnt.docx to add:
  - PlaceholderRed character style (red bold) — for [text]{custom-style="PlaceholderRed"}
  - TOC1, TOC2, TOC3 paragraph styles with TIGHT spacing (Coura-style sumario)
  - CapaHeader, CapaAutor, CapaTitulo, CapaLocal, FolhaRostoAutor, BlocoTese, Dedicatoria paragraph styles
"""
import re, zipfile, shutil, os

SRC = r'Y:/Doutorado Botan/CorePAM_Study_submitted/submission/reference_abnt.docx'
BAK = r'Y:/Doutorado Botan/CorePAM_Study_submitted/submission/reference_abnt.docx.bak'
EXTRACT = r'Y:/Doutorado Botan/CorePAM_Study_submitted/submission/_refabnt_extract'

if not os.path.exists(BAK):
    shutil.copy(SRC, BAK)
    print('Backup created:', BAK)

styles_path = os.path.join(EXTRACT, 'word', 'styles.xml')
s = open(styles_path, encoding='utf-8').read()

# Character style: PlaceholderRed (red, bold) — character-level
placeholder_char = '''<w:style w:type="character" w:customStyle="1" w:styleId="PlaceholderRed"><w:name w:val="PlaceholderRed" /><w:basedOn w:val="DefaultParagraphFont" /><w:rPr><w:b /><w:bCs /><w:color w:val="C00000" /></w:rPr></w:style>'''

# TOC1/2/3 paragraph styles — TIGHT: single spacing, zero before/after, no indent for 1, small indent for 2/3
# w:spacing w:before="0" w:after="0" w:line="240" w:lineRule="auto" = single line, no space before/after
# w:ind w:left="xxx" = left indent
toc1 = '''<w:style w:type="paragraph" w:styleId="TOC1"><w:name w:val="toc 1" /><w:basedOn w:val="Normal" /><w:next w:val="Normal" /><w:uiPriority w:val="39" /><w:pPr><w:spacing w:before="0" w:after="0" w:line="240" w:lineRule="auto" /><w:ind w:left="0" /></w:pPr><w:rPr><w:sz w:val="22" /></w:rPr></w:style>'''

toc2 = '''<w:style w:type="paragraph" w:styleId="TOC2"><w:name w:val="toc 2" /><w:basedOn w:val="Normal" /><w:next w:val="Normal" /><w:uiPriority w:val="39" /><w:pPr><w:spacing w:before="0" w:after="0" w:line="240" w:lineRule="auto" /><w:ind w:left="220" /></w:pPr><w:rPr><w:sz w:val="22" /></w:rPr></w:style>'''

toc3 = '''<w:style w:type="paragraph" w:styleId="TOC3"><w:name w:val="toc 3" /><w:basedOn w:val="Normal" /><w:next w:val="Normal" /><w:uiPriority w:val="39" /><w:pPr><w:spacing w:before="0" w:after="0" w:line="240" w:lineRule="auto" /><w:ind w:left="440" /></w:pPr><w:rPr><w:sz w:val="22" /></w:rPr></w:style>'''

# Capa/folha-rosto styles — centered
capa_header = '''<w:style w:type="paragraph" w:customStyle="1" w:styleId="CapaHeader"><w:name w:val="CapaHeader" /><w:basedOn w:val="Normal" /><w:pPr><w:spacing w:before="0" w:after="0" w:line="240" w:lineRule="auto" /><w:jc w:val="center" /></w:pPr><w:rPr><w:b /><w:bCs /><w:sz w:val="28" /></w:rPr></w:style>'''

capa_autor = '''<w:style w:type="paragraph" w:customStyle="1" w:styleId="CapaAutor"><w:name w:val="CapaAutor" /><w:basedOn w:val="Normal" /><w:pPr><w:spacing w:before="0" w:after="0" w:line="360" w:lineRule="auto" /><w:jc w:val="center" /></w:pPr><w:rPr><w:b /><w:bCs /><w:sz w:val="28" /></w:rPr></w:style>'''

capa_titulo = '''<w:style w:type="paragraph" w:customStyle="1" w:styleId="CapaTitulo"><w:name w:val="CapaTitulo" /><w:basedOn w:val="Normal" /><w:pPr><w:spacing w:before="0" w:after="0" w:line="360" w:lineRule="auto" /><w:jc w:val="center" /></w:pPr><w:rPr><w:b /><w:bCs /><w:sz w:val="32" /></w:rPr></w:style>'''

capa_local = '''<w:style w:type="paragraph" w:customStyle="1" w:styleId="CapaLocal"><w:name w:val="CapaLocal" /><w:basedOn w:val="Normal" /><w:pPr><w:spacing w:before="0" w:after="0" w:line="300" w:lineRule="auto" /><w:jc w:val="center" /></w:pPr><w:rPr><w:b /><w:bCs /><w:sz w:val="28" /></w:rPr></w:style>'''

folha_autor = '''<w:style w:type="paragraph" w:customStyle="1" w:styleId="FolhaRostoAutor"><w:name w:val="FolhaRostoAutor" /><w:basedOn w:val="Normal" /><w:pPr><w:spacing w:before="0" w:after="0" w:line="360" w:lineRule="auto" /><w:jc w:val="center" /></w:pPr><w:rPr><w:b /><w:bCs /><w:sz w:val="26" /></w:rPr></w:style>'''

# BlocoTese — ABNT: right-aligned paragraph, slightly indented, size 11
bloco_tese = '''<w:style w:type="paragraph" w:customStyle="1" w:styleId="BlocoTese"><w:name w:val="BlocoTese" /><w:basedOn w:val="Normal" /><w:pPr><w:spacing w:before="0" w:after="0" w:line="276" w:lineRule="auto" /><w:ind w:left="4536" /><w:jc w:val="both" /></w:pPr><w:rPr><w:sz w:val="22" /></w:rPr></w:style>'''

# Dedicatoria — right-aligned, italic
dedicatoria = '''<w:style w:type="paragraph" w:customStyle="1" w:styleId="Dedicatoria"><w:name w:val="Dedicatoria" /><w:basedOn w:val="Normal" /><w:pPr><w:spacing w:before="0" w:after="0" w:line="360" w:lineRule="auto" /><w:ind w:left="4536" /><w:jc w:val="right" /></w:pPr><w:rPr><w:i /><w:iCs /><w:sz w:val="24" /></w:rPr></w:style>'''

new_styles = [placeholder_char, toc1, toc2, toc3, capa_header, capa_autor, capa_titulo, capa_local, folha_autor, bloco_tese, dedicatoria]

# Insert before </w:styles>
insertion = ''.join(new_styles)
new_s = s.replace('</w:styles>', insertion + '</w:styles>')

with open(styles_path, 'w', encoding='utf-8') as f:
    f.write(new_s)

print('styles.xml patched:', len(new_styles), 'styles added')
print('Size before/after:', len(s), '->', len(new_s))

# Repack docx
# Need to preserve structure and compression. Use zipfile.
# Walk the extract dir, create new zip.
out_docx = SRC  # overwrite
tmp_out = SRC + '.tmp'

with zipfile.ZipFile(tmp_out, 'w', zipfile.ZIP_DEFLATED) as zf:
    for root, dirs, files in os.walk(EXTRACT):
        for name in files:
            full = os.path.join(root, name)
            arcname = os.path.relpath(full, EXTRACT).replace('\\', '/')
            zf.write(full, arcname)

shutil.move(tmp_out, out_docx)
print('Reference docx repacked:', out_docx)
print('Size:', os.path.getsize(out_docx))
