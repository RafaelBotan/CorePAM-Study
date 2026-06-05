"""Patch reference_abnt.docx for ABNT compliance:

1. Heading colors → black (ABNT disapproves colored headings)
2. Hyperlink color → black, no underline (cross-refs should be plain text)
3. Add header with page number (top-right, Arabic numerals)

Writes to reference_abnt.docx in-place (already has a .bak from an earlier session).
"""
import zipfile, shutil, os, re

SRC = 'reference_abnt.docx'
TMP = SRC + '.v2tmp'

with zipfile.ZipFile(SRC) as z:
    entries = {n: z.read(n) for n in z.namelist()}

styles = entries['word/styles.xml'].decode('utf-8')

# --- Remove blue colors from Heading1/2/3 and Hyperlink ---
# Heading1: color 1F3A68 → 000000
styles = styles.replace('<w:color w:val="1F3A68"/>', '<w:color w:val="000000"/>')
# Heading2: color 2E5A9C → 000000
styles = styles.replace('<w:color w:val="2E5A9C"/>', '<w:color w:val="000000"/>')
# Heading3: color 4C7CB8 → 000000
styles = styles.replace('<w:color w:val="4C7CB8"/>', '<w:color w:val="000000"/>')

# Hyperlink: change color from 0563C1 → auto, remove underline
old_hyper = ('<w:style w:type="character" w:styleId="Hyperlink"><w:name w:val="Hyperlink"/>'
             '<w:basedOn w:val="DefaultParagraphFont"/><w:uiPriority w:val="99"/>'
             '<w:unhideWhenUsed/><w:rPr><w:color w:val="0563C1"/><w:u w:val="single"/></w:rPr></w:style>')
new_hyper = ('<w:style w:type="character" w:styleId="Hyperlink"><w:name w:val="Hyperlink"/>'
             '<w:basedOn w:val="DefaultParagraphFont"/><w:uiPriority w:val="99"/>'
             '<w:unhideWhenUsed/><w:rPr><w:color w:val="000000"/></w:rPr></w:style>')
if old_hyper in styles:
    styles = styles.replace(old_hyper, new_hyper)
    print('Patched Hyperlink → black, no underline.')
else:
    print('WARN: Hyperlink exact pattern not matched; attempting loose replace.')
    styles = re.sub(
        r'<w:style w:type="character" w:styleId="Hyperlink">.*?</w:style>',
        new_hyper, styles, count=1, flags=re.S,
    )

# Ensure a FollowedHyperlink style exists and is also black with no underline
if 'w:styleId="FollowedHyperlink"' not in styles:
    fh_style = (
        '<w:style w:type="character" w:styleId="FollowedHyperlink">'
        '<w:name w:val="FollowedHyperlink"/>'
        '<w:basedOn w:val="DefaultParagraphFont"/><w:uiPriority w:val="99"/>'
        '<w:semiHidden/><w:unhideWhenUsed/>'
        '<w:rPr><w:color w:val="000000"/></w:rPr></w:style>'
    )
    styles = styles.replace('</w:styles>', fh_style + '</w:styles>')
    print('Inserted FollowedHyperlink (black).')

entries['word/styles.xml'] = styles.encode('utf-8')

# --- Add header1.xml with page number (top-right, Times 12pt) ---
HEADER_XML = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:hdr xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">
  <w:p>
    <w:pPr>
      <w:pStyle w:val="Header"/>
      <w:jc w:val="right"/>
    </w:pPr>
    <w:r>
      <w:rPr><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman"/><w:sz w:val="24"/></w:rPr>
      <w:fldChar w:fldCharType="begin"/>
    </w:r>
    <w:r>
      <w:rPr><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman"/><w:sz w:val="24"/></w:rPr>
      <w:instrText xml:space="preserve"> PAGE </w:instrText>
    </w:r>
    <w:r>
      <w:rPr><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman"/><w:sz w:val="24"/></w:rPr>
      <w:fldChar w:fldCharType="end"/>
    </w:r>
  </w:p>
</w:hdr>
'''
entries['word/header1.xml'] = HEADER_XML.encode('utf-8')

# --- Add relationship for header1.xml ---
rels = entries['word/_rels/document.xml.rels'].decode('utf-8')
if 'header1.xml' not in rels:
    # Pick an unused rId
    existing_ids = set(re.findall(r'Id="(rId\d+)"', rels))
    nxt = max(int(x[3:]) for x in existing_ids) + 1
    new_rel = (f'<Relationship Id="rId{nxt}" '
               f'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/header" '
               f'Target="header1.xml"/>')
    rels = rels.replace('</Relationships>', new_rel + '</Relationships>')
    entries['word/_rels/document.xml.rels'] = rels.encode('utf-8')
    hdr_rid = f'rId{nxt}'
    print(f'Added header1.xml relationship: {hdr_rid}')
else:
    m = re.search(r'Id="(rId\d+)"[^/]*Target="header1\.xml"', rels)
    hdr_rid = m.group(1) if m else 'rId99'
    print(f'Header1 rel already exists: {hdr_rid}')

# --- Add Override in [Content_Types].xml ---
ct = entries['[Content_Types].xml'].decode('utf-8')
if '/word/header1.xml' not in ct:
    override = ('<Override PartName="/word/header1.xml" '
                'ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.header+xml"/>')
    ct = ct.replace('</Types>', override + '</Types>')
    entries['[Content_Types].xml'] = ct.encode('utf-8')
    print('Added header1.xml override in Content_Types.')

# --- Insert headerReference in sectPr of document.xml ---
doc = entries['word/document.xml'].decode('utf-8')
hdr_ref = f'<w:headerReference w:type="default" r:id="{hdr_rid}"/>'
if 'headerReference' not in doc:
    # sectPr xmlns declarations may not include r: — we need it. Check namespaces.
    # The root <w:document> usually already declares xmlns:r=...
    # Insert headerReference as first child of <w:sectPr ...>
    doc = re.sub(
        r'(<w:sectPr[^>]*>)',
        lambda m: m.group(1) + hdr_ref,
        doc, count=1,
    )
    entries['word/document.xml'] = doc.encode('utf-8')
    print('Inserted headerReference into sectPr.')
else:
    print('sectPr already has headerReference — skipping.')

# Write
with zipfile.ZipFile(TMP, 'w', zipfile.ZIP_DEFLATED) as zout:
    for n, d in entries.items():
        zout.writestr(n, d)
shutil.move(TMP, SRC)
print(f'Patched {SRC}: {os.path.getsize(SRC)} bytes')
