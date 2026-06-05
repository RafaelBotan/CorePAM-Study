"""v3 patches on reference_abnt.docx:

1. Caption style: add <w:jc w:val="both"/> (justify) and firstLine=0.
2. Heading3/Heading4: drop themeable font references so Times New Roman is really applied;
   force black color; normalize size; keep italic.
3. Normal: ensure `<w:jc w:val="both"/>` is present.
"""
import zipfile, shutil, os, re

SRC = 'reference_abnt.docx'
TMP = SRC + '.v3tmp'

with zipfile.ZipFile(SRC) as z:
    entries = {n: z.read(n) for n in z.namelist()}

styles = entries['word/styles.xml'].decode('utf-8')

# 1) Caption: inject jc=both + firstLine=0 in pPr
old_caption_pPr = '<w:pPr><w:spacing w:line="240" w:lineRule="auto"/></w:pPr>'
new_caption_pPr = ('<w:pPr><w:spacing w:line="240" w:lineRule="auto"/>'
                   '<w:ind w:firstLine="0"/>'
                   '<w:jc w:val="both"/></w:pPr>')
if old_caption_pPr in styles:
    styles = styles.replace(old_caption_pPr, new_caption_pPr, 1)
    print('Caption: justification + firstLine=0 applied.')
else:
    print('WARN: Caption pPr not matched verbatim.')

# 2) Heading3 — strip themeable refs, force Times New Roman black
old_h3 = re.search(r'<w:style w:type="paragraph" w:styleId="Heading3">.*?</w:style>', styles, re.S).group(0)
new_h3 = (
    '<w:style w:type="paragraph" w:styleId="Heading3"><w:name w:val="heading 3"/>'
    '<w:basedOn w:val="Normal"/><w:next w:val="Normal"/><w:link w:val="Heading3Char"/>'
    '<w:uiPriority w:val="9"/><w:unhideWhenUsed/><w:qFormat/><w:rsid w:val="00FC693F"/>'
    '<w:pPr><w:keepNext/><w:keepLines/><w:spacing w:before="240" w:after="120"/>'
    '<w:ind w:firstLine="0"/><w:jc w:val="left"/><w:outlineLvl w:val="2"/></w:pPr>'
    '<w:rPr><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman" w:cs="Times New Roman"/>'
    '<w:b/><w:bCs/><w:color w:val="000000"/><w:sz w:val="24"/><w:szCs w:val="24"/></w:rPr></w:style>'
)
styles = styles.replace(old_h3, new_h3)
print('Heading3 rewritten (Times New Roman, black, 12pt, no italic).')

# 3) Heading4 — same treatment; drop blue color
old_h4 = re.search(r'<w:style w:type="paragraph" w:styleId="Heading4">.*?</w:style>', styles, re.S).group(0)
new_h4 = (
    '<w:style w:type="paragraph" w:styleId="Heading4"><w:name w:val="heading 4"/>'
    '<w:basedOn w:val="Normal"/><w:next w:val="Normal"/><w:link w:val="Heading4Char"/>'
    '<w:uiPriority w:val="9"/><w:semiHidden/><w:unhideWhenUsed/><w:qFormat/><w:rsid w:val="00FC693F"/>'
    '<w:pPr><w:keepNext/><w:keepLines/><w:spacing w:before="200" w:after="0"/>'
    '<w:ind w:firstLine="0"/><w:jc w:val="left"/><w:outlineLvl w:val="3"/></w:pPr>'
    '<w:rPr><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman" w:cs="Times New Roman"/>'
    '<w:b/><w:bCs/><w:i/><w:iCs/><w:color w:val="000000"/><w:sz w:val="24"/><w:szCs w:val="24"/></w:rPr></w:style>'
)
styles = styles.replace(old_h4, new_h4)
print('Heading4 rewritten (Times New Roman, black, 12pt, italic).')

# 4) Heading2: strip themeable refs too for consistency
old_h2 = re.search(r'<w:style w:type="paragraph" w:styleId="Heading2">.*?</w:style>', styles, re.S).group(0)
new_h2 = (
    '<w:style w:type="paragraph" w:styleId="Heading2"><w:name w:val="heading 2"/>'
    '<w:basedOn w:val="Normal"/><w:next w:val="Normal"/><w:link w:val="Heading2Char"/>'
    '<w:uiPriority w:val="9"/><w:unhideWhenUsed/><w:qFormat/><w:rsid w:val="00FC693F"/>'
    '<w:pPr><w:keepNext/><w:keepLines/><w:spacing w:before="360" w:after="120"/>'
    '<w:ind w:firstLine="0"/><w:jc w:val="left"/><w:outlineLvl w:val="1"/></w:pPr>'
    '<w:rPr><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman" w:cs="Times New Roman"/>'
    '<w:b/><w:bCs/><w:color w:val="000000"/><w:sz w:val="24"/><w:szCs w:val="26"/></w:rPr></w:style>'
)
styles = styles.replace(old_h2, new_h2)
print('Heading2 rewritten (Times New Roman, black, 12pt).')

# 5) Heading1: same treatment
old_h1 = re.search(r'<w:style w:type="paragraph" w:styleId="Heading1">.*?</w:style>', styles, re.S).group(0)
new_h1 = (
    '<w:style w:type="paragraph" w:styleId="Heading1"><w:name w:val="heading 1"/>'
    '<w:basedOn w:val="Normal"/><w:next w:val="Normal"/><w:link w:val="Heading1Char"/>'
    '<w:uiPriority w:val="9"/><w:qFormat/><w:rsid w:val="00FC693F"/>'
    '<w:pPr><w:keepNext/><w:keepLines/><w:spacing w:before="480" w:after="240"/>'
    '<w:ind w:firstLine="0"/><w:jc w:val="left"/><w:outlineLvl w:val="0"/></w:pPr>'
    '<w:rPr><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman" w:cs="Times New Roman"/>'
    '<w:b/><w:bCs/><w:color w:val="000000"/><w:sz w:val="28"/><w:szCs w:val="28"/></w:rPr></w:style>'
)
styles = styles.replace(old_h1, new_h1)
print('Heading1 rewritten (Times New Roman, black, 14pt).')

# 6) Heading{1..4}Char (linked character styles) — same theme-font stripping.
# Some renderers consult the Char style when resolving fonts inside the heading
# paragraph. If the Char style still points at majorHAnsi theme, runs can
# render in whatever theme font LibreOffice falls back to.
import re as _re
for hn, sz in [('Heading1Char', 28), ('Heading2Char', 24), ('Heading3Char', 24), ('Heading4Char', 24)]:
    m = _re.search(rf'<w:style w:type="character" w:(?:custom|default)?[^>]*w:styleId="{hn}"[^>]*>.*?</w:style>', styles, _re.S)
    if not m:
        m = _re.search(rf'<w:style [^>]*w:styleId="{hn}"[^>]*>.*?</w:style>', styles, _re.S)
    if not m:
        print(f'WARN: {hn} not found.')
        continue
    old = m.group(0)
    italic = '<w:i/><w:iCs/>' if hn == 'Heading4Char' else ''
    new = (
        f'<w:style w:type="character" w:customStyle="1" w:styleId="{hn}">'
        f'<w:name w:val="{hn.replace("Char", " Char")}"/>'
        f'<w:basedOn w:val="DefaultParagraphFont"/>'
        f'<w:link w:val="{hn.replace("Char", "")}"/>'
        f'<w:uiPriority w:val="9"/><w:rsid w:val="00FC693F"/>'
        f'<w:rPr><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman" w:cs="Times New Roman"/>'
        f'<w:b/><w:bCs/>{italic}<w:color w:val="000000"/>'
        f'<w:sz w:val="{sz}"/><w:szCs w:val="{sz}"/></w:rPr></w:style>'
    )
    styles = styles.replace(old, new)
    print(f'{hn} rewritten (Times New Roman, black).')

entries['word/styles.xml'] = styles.encode('utf-8')

with zipfile.ZipFile(TMP, 'w', zipfile.ZIP_DEFLATED) as zout:
    for n, d in entries.items():
        zout.writestr(n, d)
shutil.move(TMP, SRC)
print(f'v3 patch done: {os.path.getsize(SRC)} bytes')
