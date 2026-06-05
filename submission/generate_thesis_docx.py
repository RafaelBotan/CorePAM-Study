#!/usr/bin/env python3
"""
Generate TESE_CorePAM_v0.1.docx from the markdown thesis.
Features:
  - ABNT/UnB formatting (A4, Times New Roman 12pt, 1.5 spacing)
  - Heading styles (Heading 1-3) for auto-generated TOC
  - Auto-numbered figure and table captions via Word SEQ fields
  - Numbered bibliography references [1]-[46] in text
  - Proper page breaks between major sections
"""

import re, os, sys
from docx import Document
from docx.shared import Pt, Cm, Inches, RGBColor, Emu
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.enum.section import WD_ORIENT
from docx.enum.style import WD_STYLE_TYPE
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn, nsdecls
from docx.oxml import parse_xml
from lxml import etree
import copy

# ─── paths ───────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MD_PATH    = os.path.join(SCRIPT_DIR, "TESE_CorePAM_v0.1.md")
OUT_PATH   = os.path.join(SCRIPT_DIR, "TESE_CorePAM_v0.1.docx")
DESKTOP    = os.path.join(os.path.expanduser("~"), "Desktop", "TESE_CorePAM_v0.1.docx")

# ─── style config ────────────────────────────────────────────────────
FONT_NAME   = "Times New Roman"
FONT_SIZE   = Pt(12)
LINE_SPACE  = WD_LINE_SPACING.ONE_POINT_FIVE
MARGIN_TOP  = Cm(3)
MARGIN_BOT  = Cm(2)
MARGIN_LEFT = Cm(3)
MARGIN_RIGHT= Cm(2)


def setup_styles(doc):
    """Configure document styles for ABNT thesis."""
    # Default paragraph style
    style = doc.styles['Normal']
    font = style.font
    font.name = FONT_NAME
    font.size = FONT_SIZE
    font.color.rgb = RGBColor(0, 0, 0)
    pf = style.paragraph_format
    pf.line_spacing_rule = LINE_SPACE
    pf.space_after = Pt(6)
    pf.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY

    # Heading 1 — chapter titles
    h1 = doc.styles['Heading 1']
    h1.font.name = FONT_NAME
    h1.font.size = Pt(14)
    h1.font.bold = True
    h1.font.color.rgb = RGBColor(0, 0, 0)
    h1.paragraph_format.space_before = Pt(24)
    h1.paragraph_format.space_after = Pt(12)
    h1.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.LEFT
    h1.paragraph_format.page_break_before = False  # we control breaks manually

    # Heading 2 — sections
    h2 = doc.styles['Heading 2']
    h2.font.name = FONT_NAME
    h2.font.size = Pt(13)
    h2.font.bold = True
    h2.font.color.rgb = RGBColor(0, 0, 0)
    h2.paragraph_format.space_before = Pt(18)
    h2.paragraph_format.space_after = Pt(6)
    h2.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.LEFT

    # Heading 3 — subsections
    h3 = doc.styles['Heading 3']
    h3.font.name = FONT_NAME
    h3.font.size = Pt(12)
    h3.font.bold = True
    h3.font.color.rgb = RGBColor(0, 0, 0)
    h3.paragraph_format.space_before = Pt(12)
    h3.paragraph_format.space_after = Pt(6)
    h3.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.LEFT

    # Caption style for figures and tables
    try:
        cap = doc.styles.add_style('FigCaption', WD_STYLE_TYPE.PARAGRAPH)
    except ValueError:
        cap = doc.styles['FigCaption']
    cap.font.name = FONT_NAME
    cap.font.size = Pt(10)
    cap.font.italic = False
    cap.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.LEFT
    cap.paragraph_format.space_before = Pt(6)
    cap.paragraph_format.space_after = Pt(12)

    # Table cell style
    try:
        tc = doc.styles.add_style('TableCell', WD_STYLE_TYPE.PARAGRAPH)
    except ValueError:
        tc = doc.styles['TableCell']
    tc.font.name = FONT_NAME
    tc.font.size = Pt(10)
    tc.paragraph_format.space_before = Pt(2)
    tc.paragraph_format.space_after = Pt(2)
    tc.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE

    # Quote / epigraph style
    try:
        eq = doc.styles.add_style('Epigraph', WD_STYLE_TYPE.PARAGRAPH)
    except ValueError:
        eq = doc.styles['Epigraph']
    eq.font.name = FONT_NAME
    eq.font.size = Pt(12)
    eq.font.italic = True
    eq.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    eq.paragraph_format.left_indent = Cm(4)
    eq.paragraph_format.space_before = Pt(12)
    eq.paragraph_format.space_after = Pt(6)


def setup_page(doc):
    """Configure page layout (A4, ABNT margins)."""
    section = doc.sections[0]
    section.page_width  = Cm(21)
    section.page_height = Cm(29.7)
    section.top_margin    = MARGIN_TOP
    section.bottom_margin = MARGIN_BOT
    section.left_margin   = MARGIN_LEFT
    section.right_margin  = MARGIN_RIGHT


def add_page_break(doc):
    doc.add_page_break()


def add_toc(doc):
    """Insert a TOC field that Word will update on open."""
    p = doc.add_paragraph()
    p.style = doc.styles['Heading 1']
    run = p.add_run("SUMARIO")

    # Instruction paragraph
    pi = doc.add_paragraph()
    pi.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = pi.add_run("[Clique com botao direito e selecione 'Atualizar campo' para gerar o sumario automatico]")
    r.font.size = Pt(10)
    r.font.italic = True
    r.font.color.rgb = RGBColor(128, 128, 128)

    # Add the actual TOC field
    paragraph = doc.add_paragraph()
    run = paragraph.add_run()
    fldChar = parse_xml(f'<w:fldChar {nsdecls("w")} w:fldCharType="begin"/>')
    run._r.append(fldChar)

    run2 = paragraph.add_run()
    instrText = parse_xml(f'<w:instrText {nsdecls("w")} xml:space="preserve"> TOC \\o "1-3" \\h \\z \\u </w:instrText>')
    run2._r.append(instrText)

    run3 = paragraph.add_run()
    fldChar2 = parse_xml(f'<w:fldChar {nsdecls("w")} w:fldCharType="separate"/>')
    run3._r.append(fldChar2)

    run4 = paragraph.add_run("[Sumario sera gerado automaticamente ao atualizar campos no Word]")
    run4.font.color.rgb = RGBColor(128, 128, 128)

    run5 = paragraph.add_run()
    fldChar3 = parse_xml(f'<w:fldChar {nsdecls("w")} w:fldCharType="end"/>')
    run5._r.append(fldChar3)


def add_seq_field(paragraph, seq_name, prefix=""):
    """Add an auto-numbering SEQ field (e.g., FIGURA 1, TABELA 1)."""
    if prefix:
        run_pre = paragraph.add_run(prefix)
        run_pre.bold = True
        run_pre.font.size = Pt(10)
        run_pre.font.name = FONT_NAME

    run = paragraph.add_run()
    fldChar = parse_xml(f'<w:fldChar {nsdecls("w")} w:fldCharType="begin"/>')
    run._r.append(fldChar)

    run2 = paragraph.add_run()
    instrText = parse_xml(f'<w:instrText {nsdecls("w")} xml:space="preserve"> SEQ {seq_name} \\* ARABIC </w:instrText>')
    run2._r.append(instrText)

    run3 = paragraph.add_run()
    fldChar2 = parse_xml(f'<w:fldChar {nsdecls("w")} w:fldCharType="separate"/>')
    run3._r.append(fldChar2)

    # Placeholder number (Word will update)
    run4 = paragraph.add_run("?")
    run4.bold = True
    run4.font.size = Pt(10)
    run4.font.name = FONT_NAME

    run5 = paragraph.add_run()
    fldChar3 = parse_xml(f'<w:fldChar {nsdecls("w")} w:fldCharType="end"/>')
    run5._r.append(fldChar3)


def add_figure_caption(doc, caption_text, fig_counter):
    """Add a figure caption with auto-numbering."""
    p = doc.add_paragraph(style='FigCaption')
    add_seq_field(p, "Figura", "FIGURA ")
    run = p.add_run(f" — {caption_text}")
    run.font.size = Pt(10)
    run.font.name = FONT_NAME
    return p


def add_table_caption(doc, caption_text, tab_counter):
    """Add a table caption with auto-numbering (ABNT: caption above table)."""
    p = doc.add_paragraph(style='FigCaption')
    add_seq_field(p, "Tabela", "TABELA ")
    run = p.add_run(f" — {caption_text}")
    run.font.size = Pt(10)
    run.font.name = FONT_NAME
    return p


def add_cover_page(doc):
    """Create the cover page."""
    for _ in range(3):
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER

    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run("UNIVERSIDADE DE BRASILIA")
    r.bold = True; r.font.size = Pt(14); r.font.name = FONT_NAME

    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run("FACULDADE DE MEDICINA")
    r.bold = True; r.font.size = Pt(14); r.font.name = FONT_NAME

    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run("PROGRAMA DE POS-GRADUACAO EM CIENCIAS MEDICAS")
    r.bold = True; r.font.size = Pt(13); r.font.name = FONT_NAME

    for _ in range(4):
        doc.add_paragraph()

    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run("CorePAM: Escore prognostico de 24 genes derivado do PAM50\ncom validacao externa cross-plataforma para cancer de mama")
    r.bold = True; r.font.size = Pt(16); r.font.name = FONT_NAME

    for _ in range(2):
        doc.add_paragraph()

    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run("Versao 0.1 — Rascunho para revisao")
    r.font.size = Pt(12); r.font.name = FONT_NAME; r.italic = True

    for _ in range(3):
        doc.add_paragraph()

    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run("RAFAEL DE NEGREIROS BOTAN")
    r.bold = True; r.font.size = Pt(14); r.font.name = FONT_NAME

    doc.add_paragraph()

    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run("Orientador: Prof. Dr. Joao Batista de Sousa")
    r.font.size = Pt(12); r.font.name = FONT_NAME

    for _ in range(3):
        doc.add_paragraph()

    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run("Brasilia — DF\n2026")
    r.font.size = Pt(14); r.font.name = FONT_NAME


def add_folha_rosto(doc):
    """Create the title page (folha de rosto)."""
    for _ in range(5):
        doc.add_paragraph()

    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run("RAFAEL DE NEGREIROS BOTAN")
    r.bold = True; r.font.size = Pt(14); r.font.name = FONT_NAME

    for _ in range(4):
        doc.add_paragraph()

    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run("CorePAM: Escore prognostico de 24 genes derivado do PAM50\ncom validacao externa cross-plataforma para cancer de mama")
    r.bold = True; r.font.size = Pt(16); r.font.name = FONT_NAME

    for _ in range(3):
        doc.add_paragraph()

    # The "natureza" box (right-aligned block)
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    p.paragraph_format.left_indent = Cm(7)
    r = p.add_run("Tese apresentada ao Programa de Pos-Graduacao em Ciencias Medicas da Faculdade de Medicina da Universidade de Brasilia como requisito parcial para a obtencao do titulo de Doutor em Ciencias Medicas.")
    r.font.size = Pt(11); r.font.name = FONT_NAME

    doc.add_paragraph()

    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    p.paragraph_format.left_indent = Cm(7)
    r = p.add_run("Orientador: Prof. Dr. Joao Batista de Sousa\nFaculdade de Medicina — Departamento de Cirurgia Colorretal")
    r.font.size = Pt(11); r.font.name = FONT_NAME

    for _ in range(5):
        doc.add_paragraph()

    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run("Brasilia, 2026")
    r.font.size = Pt(14); r.font.name = FONT_NAME


def parse_markdown_table(lines):
    """Parse a markdown table into a list of rows, each a list of cells."""
    rows = []
    for line in lines:
        line = line.strip()
        if not line.startswith('|'):
            continue
        # Skip separator lines (|---|---|)
        if re.match(r'^\|[\s\-:]+\|', line):
            continue
        cells = [c.strip() for c in line.split('|')[1:-1]]
        if cells:
            rows.append(cells)
    return rows


def add_word_table(doc, rows, has_header=True):
    """Add a formatted table to the document."""
    if not rows:
        return
    n_cols = max(len(r) for r in rows)
    table = doc.add_table(rows=len(rows), cols=n_cols)
    table.style = 'Table Grid'
    table.alignment = WD_TABLE_ALIGNMENT.CENTER

    for i, row_data in enumerate(rows):
        row = table.rows[i]
        for j, cell_text in enumerate(row_data):
            if j < n_cols:
                cell = row.cells[j]
                cell.text = ""
                p = cell.paragraphs[0]
                p.style = doc.styles['TableCell']
                r = p.add_run(cell_text)
                r.font.size = Pt(9)
                r.font.name = FONT_NAME
                if i == 0 and has_header:
                    r.bold = True
                    # Gray header background
                    shading = parse_xml(f'<w:shd {nsdecls("w")} w:fill="D9E2F3"/>')
                    cell._tc.get_or_add_tcPr().append(shading)

    doc.add_paragraph()  # spacing after table


def process_inline_formatting(paragraph, text, base_bold=False, base_italic=False, base_size=None):
    """Process bold and italic markdown formatting within a paragraph."""
    if base_size is None:
        base_size = FONT_SIZE

    # Split by bold markers (**text**)
    parts = re.split(r'(\*\*[^*]+\*\*)', text)
    for part in parts:
        if part.startswith('**') and part.endswith('**'):
            r = paragraph.add_run(part[2:-2])
            r.bold = True
            r.font.name = FONT_NAME
            r.font.size = base_size
        else:
            # Check for italic within non-bold parts
            iparts = re.split(r'(\*[^*]+\*)', part)
            for ip in iparts:
                if ip.startswith('*') and ip.endswith('*') and not ip.startswith('**'):
                    r = paragraph.add_run(ip[1:-1])
                    r.italic = True
                    r.font.name = FONT_NAME
                    r.font.size = base_size
                else:
                    if ip:
                        r = paragraph.add_run(ip)
                        r.bold = base_bold
                        r.italic = base_italic
                        r.font.name = FONT_NAME
                        r.font.size = base_size


def build_docx():
    """Main function: read MD and build the DOCX."""
    print("Reading markdown thesis...")
    with open(MD_PATH, 'r', encoding='utf-8') as f:
        md_text = f.read()

    lines = md_text.split('\n')

    doc = Document()
    setup_styles(doc)
    setup_page(doc)

    # Remove default empty paragraph
    if doc.paragraphs:
        p = doc.paragraphs[0]
        p_element = p._element
        p_element.getparent().remove(p_element)

    # ─── Cover page ──────────────────────────────────────
    print("Building cover page...")
    add_cover_page(doc)
    add_page_break(doc)

    # ─── Folha de rosto ──────────────────────────────────
    print("Building folha de rosto...")
    add_folha_rosto(doc)
    add_page_break(doc)

    # ─── Process markdown content ────────────────────────
    print("Processing markdown content...")

    fig_counter = [0]
    tab_counter = [0]

    # State machine
    i = 0
    in_table = False
    table_lines = []
    current_table_caption = ""
    skip_front_matter = True  # skip until we hit EPIGRAFE
    found_epigrafe = False

    # Track sections that need page breaks before them
    page_break_sections = {
        "EPIGRAFE", "ABSTRACT", "RESUMO",
        "LISTA DE FIGURAS", "LISTA DE TABELAS",
        "LISTA DE ABREVIATURAS E SIGLAS", "SUMARIO",
        "1. INTRODUCAO", "2. OBJETIVOS", "3. METODO",
        "4. RESULTADOS", "5. DISCUSSAO", "6. CONCLUSAO",
        "7. REFERENCIAS BIBLIOGRAFICAS",
        "APENDICE A", "APENDICE B"
    }

    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        # Skip everything until EPIGRAFE (we already did cover + folha de rosto)
        if skip_front_matter:
            if stripped == "# EPIGRAFE":
                skip_front_matter = False
                found_epigrafe = True
            else:
                i += 1
                continue

        # ─── Horizontal rule (---) → skip ────────────
        if re.match(r'^---+$', stripped):
            i += 1
            continue

        # ─── Empty line ──────────────────────────────
        if not stripped:
            if in_table and table_lines:
                # End of table: flush it
                rows = parse_markdown_table(table_lines)
                if rows:
                    if current_table_caption:
                        tab_counter[0] += 1
                        add_table_caption(doc, current_table_caption, tab_counter[0])
                        current_table_caption = ""
                    add_word_table(doc, rows)
                table_lines = []
                in_table = False
            i += 1
            continue

        # ─── Table lines ────────────────────────────
        if stripped.startswith('|'):
            in_table = True
            table_lines.append(stripped)
            i += 1
            continue

        # ─── Heading 1: # Title ──────────────────────
        m1 = re.match(r'^# (.+)$', stripped)
        if m1:
            title = m1.group(1).strip()

            # Check for page break
            for pb_key in page_break_sections:
                if title.startswith(pb_key) or title == pb_key:
                    add_page_break(doc)
                    break

            # Special handling for EPIGRAFE
            if title == "EPIGRAFE":
                p = doc.add_paragraph(style='Heading 1')
                p.add_run("EPIGRAFE")
                # Read the quote
                i += 1
                quote_lines = []
                while i < len(lines):
                    ql = lines[i].strip()
                    if ql == "---":
                        break
                    if ql.startswith('>'):
                        ql = ql.lstrip('> ').strip()
                        if ql:
                            quote_lines.append(ql)
                    elif ql:
                        quote_lines.append(ql)
                    i += 1

                for ql in quote_lines:
                    p = doc.add_paragraph(style='Epigraph')
                    if ql.startswith('—') or ql.startswith('-'):
                        r = p.add_run(ql)
                        r.font.size = Pt(11)
                    else:
                        r = p.add_run(f'"{ql}"')
                        r.font.size = Pt(12)
                i += 1
                continue

            # Special: SUMARIO → insert TOC
            if title == "SUMARIO":
                add_toc(doc)
                # Skip the manual TOC lines in the MD
                i += 1
                while i < len(lines) and not lines[i].strip().startswith('#'):
                    if lines[i].strip() == '---':
                        i += 1
                        break
                    i += 1
                continue

            # Special: skip LISTA DE sections (they'll be auto in Word)
            if title.startswith("LISTA DE"):
                p = doc.add_paragraph(style='Heading 1')
                p.add_run(title)

                # Read list items
                i += 1
                while i < len(lines):
                    ll = lines[i].strip()
                    if ll.startswith('#') or ll == '---':
                        break
                    if ll.startswith('- '):
                        item_text = ll[2:]
                        p = doc.add_paragraph(style='List Bullet')
                        process_inline_formatting(p, item_text, base_size=Pt(11))
                    elif ll:
                        p = doc.add_paragraph()
                        process_inline_formatting(p, ll)
                    i += 1
                continue

            # Regular H1
            p = doc.add_paragraph(style='Heading 1')
            p.add_run(title)
            i += 1
            continue

        # ─── Heading 2: ## Title ─────────────────────
        m2 = re.match(r'^## (.+)$', stripped)
        if m2:
            title = m2.group(1).strip()
            p = doc.add_paragraph(style='Heading 2')
            p.add_run(title)
            i += 1
            continue

        # ─── Heading 3: ### Title ────────────────────
        m3 = re.match(r'^### (.+)$', stripped)
        if m3:
            title = m3.group(1).strip()
            p = doc.add_paragraph(style='Heading 3')
            p.add_run(title)
            i += 1
            continue

        # ─── Table caption detection ─────────────────
        # "**TABELA N — description**"
        tm = re.match(r'^\*\*TABELA \d+\s*[—–-]\s*(.+?)\*\*', stripped)
        if tm:
            current_table_caption = tm.group(1).strip()
            if current_table_caption.endswith('**'):
                current_table_caption = current_table_caption[:-2]
            i += 1
            continue

        # ─── Figure reference (FIGURA N — ...) ───────
        fm = re.match(r'^- FIGURA \d+\s*[—–-]\s*(.+)$', stripped)
        if fm:
            p = doc.add_paragraph(style='List Bullet')
            process_inline_formatting(p, stripped[2:], base_size=Pt(11))
            i += 1
            continue

        # ─── Bullet list ────────────────────────────
        if stripped.startswith('- '):
            text = stripped[2:]
            p = doc.add_paragraph(style='List Bullet')
            process_inline_formatting(p, text)
            i += 1
            continue

        # ─── Numbered list ──────────────────────────
        nm = re.match(r'^(\d+)\.\s+(.+)$', stripped)
        if nm and not re.match(r'^(\d+)\.\s+(INTRODUCAO|OBJETIVOS|METODO|RESULTADOS|DISCUSSAO|CONCLUSAO|REFERENCIAS|APENDICE)', stripped):
            text = nm.group(2)
            p = doc.add_paragraph(style='List Number')
            process_inline_formatting(p, text)
            i += 1
            continue

        # ─── Block quote (> text) ────────────────────
        if stripped.startswith('> '):
            text = stripped[2:].strip()
            p = doc.add_paragraph(style='Epigraph')
            r = p.add_run(text)
            r.font.name = FONT_NAME
            i += 1
            continue

        # ─── Code block (indented or ```) ────────────
        if stripped.startswith('```') or stripped.startswith('    '):
            if stripped.startswith('```'):
                i += 1
                code_lines = []
                while i < len(lines) and not lines[i].strip().startswith('```'):
                    code_lines.append(lines[i])
                    i += 1
                i += 1  # skip closing ```
                for cl in code_lines:
                    p = doc.add_paragraph()
                    p.paragraph_format.left_indent = Cm(1)
                    r = p.add_run(cl)
                    r.font.name = "Consolas"
                    r.font.size = Pt(10)
                continue
            else:
                p = doc.add_paragraph()
                p.paragraph_format.left_indent = Cm(1)
                r = p.add_run(stripped)
                r.font.name = "Consolas"
                r.font.size = Pt(10)
                i += 1
                continue

        # ─── Footnote-like lines (* or **) ───────────
        if stripped.startswith('*') and not stripped.startswith('**'):
            if stripped.startswith('*Coorte') or stripped.startswith('*Versao') or stripped.startswith('*Orientador') or stripped.startswith('*Autor') or stripped.startswith('*Programa') or stripped.startswith('*Documento'):
                p = doc.add_paragraph()
                p.paragraph_format.space_before = Pt(2)
                text = stripped.strip('*').strip()
                r = p.add_run(text)
                r.font.size = Pt(10)
                r.font.name = FONT_NAME
                r.italic = True
                i += 1
                continue

        # ─── Regular paragraph ───────────────────────
        p = doc.add_paragraph()
        process_inline_formatting(p, stripped)
        i += 1

    # Flush any remaining table
    if table_lines:
        rows = parse_markdown_table(table_lines)
        if rows:
            if current_table_caption:
                tab_counter[0] += 1
                add_table_caption(doc, current_table_caption, tab_counter[0])
            add_word_table(doc, rows)

    # ─── Add page numbers ────────────────────────────
    print("Adding page numbers...")
    for section in doc.sections:
        footer = section.footer
        footer.is_linked_to_previous = False
        p = footer.paragraphs[0] if footer.paragraphs else footer.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER

        run = p.add_run()
        fldChar = parse_xml(f'<w:fldChar {nsdecls("w")} w:fldCharType="begin"/>')
        run._r.append(fldChar)

        run2 = p.add_run()
        instrText = parse_xml(f'<w:instrText {nsdecls("w")} xml:space="preserve"> PAGE </w:instrText>')
        run2._r.append(instrText)

        run3 = p.add_run()
        fldChar2 = parse_xml(f'<w:fldChar {nsdecls("w")} w:fldCharType="separate"/>')
        run3._r.append(fldChar2)

        run4 = p.add_run("1")
        run4.font.size = Pt(10)
        run4.font.name = FONT_NAME

        run5 = p.add_run()
        fldChar3 = parse_xml(f'<w:fldChar {nsdecls("w")} w:fldCharType="end"/>')
        run5._r.append(fldChar3)

    # ─── Set document to update fields on open ────────
    # This makes Word update TOC and SEQ fields when opened
    settings = doc.settings.element
    update_fields = parse_xml(f'<w:updateFields {nsdecls("w")} w:val="true"/>')
    settings.append(update_fields)

    # ─── Save ─────────────────────────────────────────
    print(f"Saving to {OUT_PATH}...")
    doc.save(OUT_PATH)

    # Copy to desktop
    import shutil
    print(f"Copying to {DESKTOP}...")
    shutil.copy2(OUT_PATH, DESKTOP)

    print(f"\nDone! Files saved:")
    print(f"  1. {OUT_PATH}")
    print(f"  2. {DESKTOP}")
    print(f"\nIMPORTANTE: Ao abrir no Word, clique 'Sim' quando perguntar")
    print(f"  se deseja atualizar campos. Isso ativara o sumario automatico")
    print(f"  e a numeracao de figuras/tabelas.")


if __name__ == "__main__":
    build_docx()
