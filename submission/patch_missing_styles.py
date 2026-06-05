"""Add missing Pandoc-expected styles to reference_abnt.docx.

Root cause of broken table rendering: reference_abnt.docx lacks the Table, Compact,
Caption, ImageCaption, TableCaption, Hyperlink styles that Pandoc's docx output uses.
Without these, LibreOffice falls back to defaults that break table layout.
"""
import zipfile, re, shutil, os

DOCX = 'reference_abnt.docx'
TMP = DOCX + '.tmp'

# Styles to inject (keeping single-line format to be safe)
styles_to_add = {
    'Table': '<w:style w:type="table" w:default="1" w:styleId="Table"><w:name w:val="Table"/><w:basedOn w:val="TableNormal"/><w:qFormat/><w:tblPr><w:tblInd w:w="0" w:type="dxa"/><w:tblBorders><w:top w:val="single" w:sz="4" w:space="0" w:color="auto"/><w:bottom w:val="single" w:sz="4" w:space="0" w:color="auto"/><w:insideH w:val="single" w:sz="4" w:space="0" w:color="BFBFBF"/></w:tblBorders><w:tblCellMar><w:top w:w="40" w:type="dxa"/><w:left w:w="108" w:type="dxa"/><w:bottom w:w="40" w:type="dxa"/><w:right w:w="108" w:type="dxa"/></w:tblCellMar></w:tblPr><w:tblStylePr w:type="firstRow"><w:rPr><w:b/></w:rPr><w:tblPr/><w:tcPr><w:tcBorders><w:bottom w:val="single" w:sz="6" w:space="0" w:color="auto"/></w:tcBorders><w:vAlign w:val="bottom"/></w:tcPr></w:tblStylePr></w:style>',
    'Compact': '<w:style w:type="paragraph" w:customStyle="1" w:styleId="Compact"><w:name w:val="Compact"/><w:basedOn w:val="BodyText"/><w:qFormat/><w:pPr><w:spacing w:before="36" w:after="36" w:line="240" w:lineRule="auto"/><w:ind w:firstLine="0"/></w:pPr></w:style>',
    'Caption': '<w:style w:type="paragraph" w:styleId="Caption"><w:name w:val="Caption"/><w:basedOn w:val="Normal"/><w:qFormat/><w:pPr><w:spacing w:before="120" w:after="120" w:line="276" w:lineRule="auto"/><w:ind w:firstLine="0"/><w:jc w:val="left"/></w:pPr><w:rPr><w:sz w:val="20"/><w:szCs w:val="20"/></w:rPr></w:style>',
    'ImageCaption': '<w:style w:type="paragraph" w:customStyle="1" w:styleId="ImageCaption"><w:name w:val="Image Caption"/><w:basedOn w:val="Caption"/><w:qFormat/></w:style>',
    'TableCaption': '<w:style w:type="paragraph" w:customStyle="1" w:styleId="TableCaption"><w:name w:val="Table Caption"/><w:basedOn w:val="Caption"/><w:qFormat/><w:pPr><w:keepNext/></w:pPr></w:style>',
    'Hyperlink': '<w:style w:type="character" w:styleId="Hyperlink"><w:name w:val="Hyperlink"/><w:basedOn w:val="DefaultParagraphFont"/><w:uiPriority w:val="99"/><w:unhideWhenUsed/><w:rPr><w:color w:val="0563C1"/><w:u w:val="single"/></w:rPr></w:style>',
}

with zipfile.ZipFile(DOCX) as zin:
    entries = {n: zin.read(n) for n in zin.namelist()}

styles_xml = entries['word/styles.xml'].decode('utf-8')

added = 0
to_insert = []
for sid, xml in styles_to_add.items():
    if re.search(r'w:styleId="' + sid + r'"', styles_xml):
        print(f'  [{sid}] already present — skip')
    else:
        to_insert.append(xml)
        added += 1
        print(f'  [{sid}] will add')

if to_insert:
    insertion = ''.join(to_insert)
    styles_xml = styles_xml.replace('</w:styles>', insertion + '</w:styles>')
    entries['word/styles.xml'] = styles_xml.encode('utf-8')

    # Validate XML
    import xml.etree.ElementTree as ET
    try:
        ET.fromstring(styles_xml)
        print('styles.xml valid.')
    except ET.ParseError as e:
        print('XML PARSE ERROR:', e)
        raise

    with zipfile.ZipFile(TMP, 'w', zipfile.ZIP_DEFLATED) as zout:
        for name, data in entries.items():
            zout.writestr(name, data)
    shutil.move(TMP, DOCX)
    print(f'\n{added} styles injected. docx size: {os.path.getsize(DOCX)}')
else:
    print('No changes needed.')
