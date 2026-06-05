"""Fix approval-page paragraph alignment after DOCX render.

The thesis reference style justifies Normal paragraphs. In the approval sheet,
the banca entries use manual line breaks to keep the page compact; justified
paragraphs make those short lines expand with large spaces between words.
This patch sets the approval-sheet paragraphs to explicit left/center alignment.
"""
import html
import os
import re
import shutil
import zipfile
import xml.etree.ElementTree as ET

DOCX = os.environ.get("DOCX", "TESE_CorePAM_FINAL_envio_Dr_Joao_impressao.docx")
TMP = DOCX + ".approvaltmp"


def paragraph_blocks(xml):
    pos = 0
    while True:
        starts = [s for s in (xml.find("<w:p>", pos), xml.find("<w:p ", pos)) if s >= 0]
        if not starts:
            break
        start = min(starts)
        end = xml.find("</w:p>", start)
        if end < 0:
            break
        end += len("</w:p>")
        yield start, end, xml[start:end]
        pos = end


def text_of(block):
    parts = re.findall(r"<w:t[^>]*>(.*?)</w:t>", block)
    return html.unescape(" ".join(parts)).replace("\u00a0", " ").strip()


def set_alignment(block, val, clear_indent=False):
    ppr_match = re.search(r"<w:pPr\b[^>]*>[\s\S]*?</w:pPr>", block)
    jc = f'<w:jc w:val="{val}"/>'
    if ppr_match:
        ppr = ppr_match.group(0)
        ppr = re.sub(r"<w:jc\b[^/]*/>", "", ppr)
        if clear_indent:
            ppr = re.sub(r"<w:ind\b[^/]*/>", "", ppr)
            ppr = ppr.replace("</w:pPr>", '<w:ind w:left="0" w:right="0" w:firstLine="0"/>' + "</w:pPr>")
        ppr = ppr.replace("</w:pPr>", jc + "</w:pPr>")
        return block[: ppr_match.start()] + ppr + block[ppr_match.end() :]

    open_end = block.find(">")
    return block[: open_end + 1] + f"<w:pPr>{jc}</w:pPr>" + block[open_end + 1 :]


with zipfile.ZipFile(DOCX) as zin:
    entries = {name: zin.read(name) for name in zin.namelist()}

doc = entries["word/document.xml"].decode("utf-8")
blocks = list(paragraph_blocks(doc))
texts = [text_of(block) for _, _, block in blocks]

start_idx = next((i for i, txt in enumerate(texts) if "FOLHA DE APROVAÇÃO" in txt), None)
end_idx = next((i for i, txt in enumerate(texts) if "DEDICATÓRIA" in txt), None)
if start_idx is None or end_idx is None or end_idx <= start_idx:
    raise SystemExit("Approval-page range not found.")

rebuilt = []
last = 0
patched = 0
for i, (start, end, block) in enumerate(blocks):
    rebuilt.append(doc[last:start])
    txt = texts[i]
    if start_idx <= i < end_idx:
        if txt == "BANCA EXAMINADORA":
            block = set_alignment(block, "center", clear_indent=True)
        elif txt == "RAFAEL DE NEGREIROS BOTAN" or txt.startswith("CorePAM:"):
            block = set_alignment(block, "center")
        else:
            block = set_alignment(block, "left")
        patched += 1
    rebuilt.append(block)
    last = end
rebuilt.append(doc[last:])
doc = "".join(rebuilt)

ET.fromstring(doc)
entries["word/document.xml"] = doc.encode("utf-8")

with zipfile.ZipFile(TMP, "w", zipfile.ZIP_DEFLATED) as zout:
    for name, data in entries.items():
        zout.writestr(name, data)
shutil.move(TMP, DOCX)
print(f"Approval-page alignment patched in {DOCX}; paragraphs={patched}")
