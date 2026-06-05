"""Patch reference_abnt.docx — add AbrevList style and set H1/H2/H3 to blue shades."""
import re, zipfile, shutil, os

SRC = 'reference_abnt.docx'
EXTRACT = '_refabnt_extract'

styles_path = os.path.join(EXTRACT, 'word', 'styles.xml')
s = open(styles_path, encoding='utf-8').read()

# AbrevList: tight single-paragraph style for compact abbreviations listing
abrev_list = '<w:style w:type="paragraph" w:customStyle="1" w:styleId="AbrevList"><w:name w:val="AbrevList" /><w:basedOn w:val="Normal" /><w:pPr><w:spacing w:before="0" w:after="120" w:line="300" w:lineRule="auto" /><w:ind w:firstLine="0" /><w:jc w:val="both" /></w:pPr><w:rPr><w:sz w:val="22" /></w:rPr></w:style>'

if 'w:styleId="AbrevList"' not in s:
    s = s.replace('</w:styles>', abrev_list + '</w:styles>')
    print('AbrevList added')

# Change H1/H2/H3 color to blue shades
def set_color(xml, style_id, color):
    pattern = r'(<w:style[^>]*w:styleId="' + style_id + r'".*?</w:style>)'
    m = re.search(pattern, xml, re.S)
    if not m:
        print(style_id, 'not found')
        return xml
    block = m.group(1)
    new_block = re.sub(r'<w:color w:val="[0-9A-Fa-f]+"\s*/>', '<w:color w:val="' + color + '"/>', block)
    if new_block == block:
        # insert color if no color
        new_block = block.replace('<w:b/>', '<w:b/><w:color w:val="'+color+'"/>', 1)
    return xml.replace(block, new_block)

s = set_color(s, 'Heading1', '1F3A68')
s = set_color(s, 'Heading2', '2E5A9C')
s = set_color(s, 'Heading3', '4C7CB8')

open(styles_path, 'w', encoding='utf-8').write(s)
print('Heading colors updated')

# Repack
tmp = SRC + '.tmp'
with zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as zf:
    for root, dirs, files in os.walk(EXTRACT):
        for name in files:
            full = os.path.join(root, name)
            arc = os.path.relpath(full, EXTRACT).replace('\\', '/')
            zf.write(full, arc)
shutil.move(tmp, SRC)
print('repacked:', os.path.getsize(SRC))
