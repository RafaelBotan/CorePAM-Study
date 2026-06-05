"""Force visible blue/underlined hyperlinks in the final thesis DOCX.

Pandoc creates real DOCX hyperlinks for DOI/URL references, but the reference
DOCX style can render them as black in LibreOffice PDF export. This post-process
adds direct WordprocessingML formatting to every hyperlink run so the exported
PDF clearly signals clickable links while remaining readable in print.
"""
import os
import re
import shutil
import zipfile
import xml.etree.ElementTree as ET


DOCX = os.environ.get("DOCX", "TESE_CorePAM_FINAL_envio_Dr_Joao_impressao.docx")
TMP = DOCX + ".hyperlinktmp"

with zipfile.ZipFile(DOCX) as zin:
    entries = {name: zin.read(name) for name in zin.namelist()}

doc = entries["word/document.xml"].decode("utf-8")


def ensure_run_format(run_xml: str) -> str:
    """Add blue color and single underline to one <w:r> run."""
    if "<w:rPr" in run_xml:
        start = run_xml.find("<w:rPr")
        open_end = run_xml.find(">", start) + 1
        close = run_xml.find("</w:rPr>", open_end)
        rpr = run_xml[open_end:close]
        rpr = re.sub(r"<w:color\b[^>]*/>", "", rpr)
        rpr = re.sub(r"<w:u\b[^>]*/>", "", rpr)
        rpr += '<w:color w:val="0563C1"/><w:u w:val="single"/>'
        return run_xml[:open_end] + rpr + run_xml[close:]
    insert = run_xml.find(">") + 1
    rpr = '<w:rPr><w:rStyle w:val="Hyperlink"/><w:color w:val="0563C1"/><w:u w:val="single"/></w:rPr>'
    return run_xml[:insert] + rpr + run_xml[insert:]


def patch_hyperlink(match: re.Match) -> str:
    block = match.group(0)
    return re.sub(r"<w:r\b[^>]*>.*?</w:r>", lambda m: ensure_run_format(m.group(0)), block, flags=re.S)


doc = re.sub(r"<w:hyperlink\b[^>]*>.*?</w:hyperlink>", patch_hyperlink, doc, flags=re.S)

# Also make the style definition explicit when present.
if "word/styles.xml" in entries:
    styles = entries["word/styles.xml"].decode("utf-8")
    style_match = re.search(
        r'<w:style\b(?=[^>]*w:styleId="Hyperlink")[\s\S]*?</w:style>',
        styles,
    )
    if style_match:
        style = style_match.group(0)
        if "<w:rPr>" in style:
            style = re.sub(r"<w:color\b[^>]*/>", "", style)
            style = re.sub(r"<w:u\b[^>]*/>", "", style)
            style = style.replace("</w:rPr>", '<w:color w:val="0563C1"/><w:u w:val="single"/></w:rPr>')
        else:
            style = style.replace("</w:style>", '<w:rPr><w:color w:val="0563C1"/><w:u w:val="single"/></w:rPr></w:style>')
        styles = styles[: style_match.start()] + style + styles[style_match.end() :]
        ET.fromstring(styles)
        entries["word/styles.xml"] = styles.encode("utf-8")

ET.fromstring(doc)
entries["word/document.xml"] = doc.encode("utf-8")

with zipfile.ZipFile(TMP, "w", zipfile.ZIP_DEFLATED) as zout:
    for name, data in entries.items():
        zout.writestr(name, data)

shutil.move(TMP, DOCX)
print(f"Hyperlinks forced blue/underlined in {DOCX}")
