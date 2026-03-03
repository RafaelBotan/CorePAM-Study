"""Create ABNT/UnB reference DOCX with blue headings, table styles, and caption styles for Quarto."""
from docx import Document
from docx.shared import Pt, Cm, RGBColor, Emu
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.enum.style import WD_STYLE_TYPE
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn, nsdecls
from docx.oxml import parse_xml
import copy

doc = Document()

# Page setup ABNT/UnB: 3cm left/top, 2cm right/bottom
for section in doc.sections:
    section.top_margin = Cm(3)
    section.bottom_margin = Cm(2)
    section.left_margin = Cm(3)
    section.right_margin = Cm(2)
    section.page_height = Cm(29.7)
    section.page_width = Cm(21)

style = doc.styles
BLUE = RGBColor(0, 51, 153)  # dark academic blue

# ---------------------------------------------------------------------------
# Normal: Times New Roman 12pt, 1.5 spacing, 1.25cm indent, justified
# ---------------------------------------------------------------------------
normal = style['Normal']
normal.font.name = 'Times New Roman'
normal.font.size = Pt(12)
normal.font.color.rgb = RGBColor(0, 0, 0)
pf = normal.paragraph_format
pf.space_before = Pt(0)
pf.space_after = Pt(6)
pf.line_spacing_rule = WD_LINE_SPACING.ONE_POINT_FIVE
pf.first_line_indent = Cm(1.25)
pf.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY

# ---------------------------------------------------------------------------
# Heading 1: 14pt bold BLUE
# ---------------------------------------------------------------------------
h1 = style['Heading 1']
h1.font.name = 'Times New Roman'
h1.font.size = Pt(14)
h1.font.bold = True
h1.font.color.rgb = BLUE
h1pf = h1.paragraph_format
h1pf.space_before = Pt(24)
h1pf.space_after = Pt(12)
h1pf.first_line_indent = Cm(0)
h1pf.alignment = WD_ALIGN_PARAGRAPH.LEFT
h1pf.line_spacing_rule = WD_LINE_SPACING.ONE_POINT_FIVE

# ---------------------------------------------------------------------------
# Heading 2: 12pt bold BLUE
# ---------------------------------------------------------------------------
h2 = style['Heading 2']
h2.font.name = 'Times New Roman'
h2.font.size = Pt(12)
h2.font.bold = True
h2.font.color.rgb = BLUE
h2pf = h2.paragraph_format
h2pf.space_before = Pt(18)
h2pf.space_after = Pt(6)
h2pf.first_line_indent = Cm(0)
h2pf.alignment = WD_ALIGN_PARAGRAPH.LEFT
h2pf.line_spacing_rule = WD_LINE_SPACING.ONE_POINT_FIVE

# ---------------------------------------------------------------------------
# Heading 3: 12pt bold italic BLUE
# ---------------------------------------------------------------------------
h3 = style['Heading 3']
h3.font.name = 'Times New Roman'
h3.font.size = Pt(12)
h3.font.bold = True
h3.font.italic = True
h3.font.color.rgb = BLUE
h3pf = h3.paragraph_format
h3pf.space_before = Pt(12)
h3pf.space_after = Pt(6)
h3pf.first_line_indent = Cm(0)
h3pf.line_spacing_rule = WD_LINE_SPACING.ONE_POINT_FIVE

# ---------------------------------------------------------------------------
# Title: 14pt bold BLUE centered
# ---------------------------------------------------------------------------
title = style['Title']
title.font.name = 'Times New Roman'
title.font.size = Pt(14)
title.font.bold = True
title.font.color.rgb = BLUE
tpf = title.paragraph_format
tpf.alignment = WD_ALIGN_PARAGRAPH.CENTER
tpf.space_after = Pt(24)
tpf.first_line_indent = Cm(0)

# ---------------------------------------------------------------------------
# Caption: 10pt centered, no first-line indent (ABNT captions)
# Pandoc/Quarto uses the "Caption" style for figure/table captions.
# We also set Table Caption and Image Caption if they exist.
# ---------------------------------------------------------------------------
for sname in ['Caption', 'Table Caption', 'Image Caption']:
    if sname in [s.name for s in style]:
        cap = style[sname]
        cap.font.name = 'Times New Roman'
        cap.font.size = Pt(10)
        cap.font.color.rgb = RGBColor(0, 0, 0)
        cap_pf = cap.paragraph_format
        cap_pf.alignment = WD_ALIGN_PARAGRAPH.CENTER
        cap_pf.first_line_indent = Cm(0)
        cap_pf.space_before = Pt(6)
        cap_pf.space_after = Pt(6)
        cap_pf.line_spacing_rule = WD_LINE_SPACING.SINGLE

# Also create a custom "Table Caption" style if it does not exist, since
# Quarto may use it depending on configuration.
existing_names = [s.name for s in style]
if 'Table Caption' not in existing_names:
    tc_style = style.add_style('Table Caption', WD_STYLE_TYPE.PARAGRAPH)
    tc_style.base_style = style['Caption']
    tc_style.font.name = 'Times New Roman'
    tc_style.font.size = Pt(10)
    tc_style.font.color.rgb = RGBColor(0, 0, 0)
    tc_pf = tc_style.paragraph_format
    tc_pf.alignment = WD_ALIGN_PARAGRAPH.CENTER
    tc_pf.first_line_indent = Cm(0)
    tc_pf.space_before = Pt(6)
    tc_pf.space_after = Pt(6)
    tc_pf.line_spacing_rule = WD_LINE_SPACING.SINGLE

if 'Image Caption' not in existing_names:
    ic_style = style.add_style('Image Caption', WD_STYLE_TYPE.PARAGRAPH)
    ic_style.base_style = style['Caption']
    ic_style.font.name = 'Times New Roman'
    ic_style.font.size = Pt(10)
    ic_style.font.color.rgb = RGBColor(0, 0, 0)
    ic_pf = ic_style.paragraph_format
    ic_pf.alignment = WD_ALIGN_PARAGRAPH.CENTER
    ic_pf.first_line_indent = Cm(0)
    ic_pf.space_before = Pt(6)
    ic_pf.space_after = Pt(6)
    ic_pf.line_spacing_rule = WD_LINE_SPACING.SINGLE

# ---------------------------------------------------------------------------
# Figure paragraph style: centered, no indent, spacing 6pt before/after
# Quarto renders images inside a paragraph. If we create a custom style,
# we can assign it; but more practically, Quarto uses "Body Text" or
# "First Paragraph" for image paragraphs. We create a dedicated style
# and also ensure "Body Text" has no indent for compatibility.
# ---------------------------------------------------------------------------
if 'Body Text' in existing_names:
    bt = style['Body Text']
    bt.font.name = 'Times New Roman'
    bt.font.size = Pt(12)
    bt_pf = bt.paragraph_format
    bt_pf.alignment = WD_ALIGN_PARAGRAPH.CENTER
    bt_pf.first_line_indent = Cm(0)
    bt_pf.space_before = Pt(6)
    bt_pf.space_after = Pt(6)
    bt_pf.line_spacing_rule = WD_LINE_SPACING.SINGLE

# "Figure" paragraph style (custom) -- centered, no indent
if 'Figure' not in existing_names:
    fig_style = style.add_style('Figure', WD_STYLE_TYPE.PARAGRAPH)
    fig_style.base_style = style['Normal']
    fig_style.font.name = 'Times New Roman'
    fig_style.font.size = Pt(12)
    fig_pf = fig_style.paragraph_format
    fig_pf.alignment = WD_ALIGN_PARAGRAPH.CENTER
    fig_pf.first_line_indent = Cm(0)
    fig_pf.space_before = Pt(6)
    fig_pf.space_after = Pt(6)
    fig_pf.line_spacing_rule = WD_LINE_SPACING.SINGLE

# ---------------------------------------------------------------------------
# "First Paragraph" style (Pandoc uses this for paragraphs after headings)
# -- no indent, same as Normal otherwise
# ---------------------------------------------------------------------------
if 'First Paragraph' not in existing_names:
    fp_style = style.add_style('First Paragraph', WD_STYLE_TYPE.PARAGRAPH)
    fp_style.base_style = style['Normal']
    fp_style.font.name = 'Times New Roman'
    fp_style.font.size = Pt(12)
    fp_pf = fp_style.paragraph_format
    fp_pf.first_line_indent = Cm(0)
    fp_pf.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    fp_pf.line_spacing_rule = WD_LINE_SPACING.ONE_POINT_FIVE

# ---------------------------------------------------------------------------
# Block Text (for long citations): 10pt, single spacing, 4cm left indent
# ---------------------------------------------------------------------------
if 'Block Text' in existing_names:
    bt2 = style['Block Text']
    bt2.font.name = 'Times New Roman'
    bt2.font.size = Pt(10)
    bt2pf = bt2.paragraph_format
    bt2pf.left_indent = Cm(4)
    bt2pf.first_line_indent = Cm(0)
    bt2pf.line_spacing_rule = WD_LINE_SPACING.SINGLE

# ---------------------------------------------------------------------------
# TOC Heading
# ---------------------------------------------------------------------------
if 'TOC Heading' in existing_names:
    toch = style['TOC Heading']
    toch.font.name = 'Times New Roman'
    toch.font.size = Pt(14)
    toch.font.bold = True
    toch.font.color.rgb = BLUE

# ---------------------------------------------------------------------------
# Table Grid style: Times New Roman 10pt, single spacing, no first-line
# indent, with proper borders (academic three-line table style).
# ---------------------------------------------------------------------------
# python-docx has a built-in "Table Grid" style. We modify it via XML to get
# the exact border configuration we want (horizontal lines only: top of
# header, bottom of header, bottom of table -- "three-line" / "bookman" style).
# ---------------------------------------------------------------------------
tbl_grid = None
for s in style:
    if s.name == 'Table Grid':
        tbl_grid = s
        break

if tbl_grid is not None:
    # Set font properties on the table style's linked paragraph properties
    tbl_grid.font.name = 'Times New Roman'
    tbl_grid.font.size = Pt(10)
    tbl_grid.font.color.rgb = RGBColor(0, 0, 0)

    # Access the underlying XML element for the table style
    tbl_style_elem = tbl_grid.element

    # Set paragraph properties (single spacing, no indent, centered)
    # We need to add pPr to the style
    existing_pPr = tbl_style_elem.find(qn('w:pPr'))
    if existing_pPr is not None:
        tbl_style_elem.remove(existing_pPr)

    pPr_xml = (
        '<w:pPr %s>'
        '  <w:spacing w:after="0" w:before="0" w:line="240" w:lineRule="auto"/>'
        '  <w:ind w:firstLine="0"/>'
        '  <w:jc w:val="center"/>'
        '</w:pPr>' % nsdecls('w')
    )
    tbl_style_elem.append(parse_xml(pPr_xml))

    # Set table-level properties: borders (three-line academic style)
    # Remove existing tblPr if any and replace
    existing_tblPr = tbl_style_elem.find(qn('w:tblPr'))
    if existing_tblPr is not None:
        tbl_style_elem.remove(existing_tblPr)

    # Border size: 4 = 0.5pt, 12 = 1.5pt, 8 = 1pt
    tblPr_xml = (
        '<w:tblPr %s>'
        '  <w:tblBorders>'
        '    <w:top w:val="single" w:sz="12" w:space="0" w:color="000000"/>'
        '    <w:bottom w:val="single" w:sz="12" w:space="0" w:color="000000"/>'
        '    <w:insideH w:val="single" w:sz="4" w:space="0" w:color="000000"/>'
        '  </w:tblBorders>'
        '  <w:jc w:val="center"/>'
        '  <w:tblCellMar>'
        '    <w:top w:w="57" w:type="dxa"/>'
        '    <w:start w:w="108" w:type="dxa"/>'
        '    <w:bottom w:w="57" w:type="dxa"/>'
        '    <w:end w:w="108" w:type="dxa"/>'
        '  </w:tblCellMar>'
        '</w:tblPr>' % nsdecls('w')
    )
    tbl_style_elem.append(parse_xml(tblPr_xml))

    # Conditional formatting for the header row (first row):
    # thicker bottom border to separate header from body
    # Remove existing tblStylePr elements
    for existing_tsp in tbl_style_elem.findall(qn('w:tblStylePr')):
        tbl_style_elem.remove(existing_tsp)

    # Header row (firstRow): bold text + thicker bottom border
    firstRow_xml = (
        '<w:tblStylePr w:type="firstRow" %s>'
        '  <w:rPr>'
        '    <w:b/>'
        '    <w:bCs/>'
        '  </w:rPr>'
        '  <w:tblPr/>'
        '  <w:tcPr>'
        '    <w:tcBorders>'
        '      <w:top w:val="single" w:sz="12" w:space="0" w:color="000000"/>'
        '      <w:bottom w:val="single" w:sz="8" w:space="0" w:color="000000"/>'
        '    </w:tcBorders>'
        '  </w:tcPr>'
        '</w:tblStylePr>' % nsdecls('w')
    )
    tbl_style_elem.append(parse_xml(firstRow_xml))

    # Last row: thicker bottom border
    lastRow_xml = (
        '<w:tblStylePr w:type="lastRow" %s>'
        '  <w:tblPr/>'
        '  <w:tcPr>'
        '    <w:tcBorders>'
        '      <w:bottom w:val="single" w:sz="12" w:space="0" w:color="000000"/>'
        '    </w:tcBorders>'
        '  </w:tcPr>'
        '</w:tblStylePr>' % nsdecls('w')
    )
    tbl_style_elem.append(parse_xml(lastRow_xml))

    # Set the table style to use firstRow and lastRow conditional formatting
    # by setting tblLook (this tells Word to apply the conditional styles)
    # We need to set this in tblPr
    tblPr_elem = tbl_style_elem.find(qn('w:tblPr'))
    if tblPr_elem is not None:
        tblLook = parse_xml(
            '<w:tblLook %s w:val="04A0" w:firstRow="1" w:lastRow="1" '
            'w:firstColumn="0" w:lastColumn="0" w:noHBand="0" w:noVBand="1"/>'
            % nsdecls('w')
        )
        tblPr_elem.append(tblLook)

    # Table row properties: keep rows together (prevent page break inside row)
    trPr_xml = (
        '<w:trPr %s>'
        '  <w:cantSplit/>'
        '</w:trPr>' % nsdecls('w')
    )
    tbl_style_elem.append(parse_xml(trPr_xml))

# ---------------------------------------------------------------------------
# Also create a "Plain Table" style (no borders) for cases where we want
# a clean layout table (e.g., figure placement).
# ---------------------------------------------------------------------------
existing_table_styles = [s.name for s in style if s.type == WD_STYLE_TYPE.TABLE]
if 'Plain Table' not in existing_table_styles:
    # We cannot easily add table styles with python-docx's high-level API,
    # so we'll skip this if it does not exist. Table Grid is the main one.
    pass

# ---------------------------------------------------------------------------
# Compact paragraph style (for table body text -- 10pt, single, no indent)
# ---------------------------------------------------------------------------
if 'Compact' not in existing_names:
    compact = style.add_style('Compact', WD_STYLE_TYPE.PARAGRAPH)
    compact.base_style = style['Normal']
    compact.font.name = 'Times New Roman'
    compact.font.size = Pt(10)
    compact.font.color.rgb = RGBColor(0, 0, 0)
    compact_pf = compact.paragraph_format
    compact_pf.space_before = Pt(0)
    compact_pf.space_after = Pt(0)
    compact_pf.first_line_indent = Cm(0)
    compact_pf.line_spacing_rule = WD_LINE_SPACING.SINGLE
    compact_pf.alignment = WD_ALIGN_PARAGRAPH.LEFT

# ---------------------------------------------------------------------------
# "Source Note" style -- for table source/notes below table (ABNT)
# 10pt, single spacing, no indent
# ---------------------------------------------------------------------------
if 'Source Note' not in existing_names:
    sn = style.add_style('Source Note', WD_STYLE_TYPE.PARAGRAPH)
    sn.base_style = style['Normal']
    sn.font.name = 'Times New Roman'
    sn.font.size = Pt(10)
    sn.font.color.rgb = RGBColor(0, 0, 0)
    sn_pf = sn.paragraph_format
    sn_pf.space_before = Pt(2)
    sn_pf.space_after = Pt(6)
    sn_pf.first_line_indent = Cm(0)
    sn_pf.line_spacing_rule = WD_LINE_SPACING.SINGLE
    sn_pf.alignment = WD_ALIGN_PARAGRAPH.LEFT

# ---------------------------------------------------------------------------
# Header / Footer styles
# ---------------------------------------------------------------------------
for hf_name in ['Header', 'Footer']:
    if hf_name in existing_names:
        hf = style[hf_name]
        hf.font.name = 'Times New Roman'
        hf.font.size = Pt(10)
        hf_pf = hf.paragraph_format
        hf_pf.first_line_indent = Cm(0)

# ---------------------------------------------------------------------------
# Save the reference document
# ---------------------------------------------------------------------------
out_path = 'Y:/Phd-Genomic-claude/tese/reference_abnt_v2.docx'
doc.save(out_path)
print(f"reference_abnt_v2.docx created at {out_path}")
print("Styles included:")
print("  - Normal (12pt, 1.5 spacing, 1.25cm indent, justified)")
print("  - Heading 1/2/3 (blue, bold)")
print("  - Title (14pt, blue, centered)")
print("  - Caption / Table Caption / Image Caption (10pt, centered, no indent)")
print("  - Figure (centered, no indent, 6pt spacing)")
print("  - First Paragraph (no indent)")
print("  - Body Text (centered, for images)")
print("  - Table Grid (10pt, single spacing, three-line borders, cell padding, cantSplit)")
print("  - Compact (10pt, single, for table body)")
print("  - Source Note (10pt, for table notes)")
print("  - Block Text (10pt, 4cm indent, for citations)")
print("  - Header / Footer (10pt)")
