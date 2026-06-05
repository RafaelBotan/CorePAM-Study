"""Hide visible page numbers before the first INTRODUÇÃO page in the final PDF.

ABNT-style pagination counts front matter but displays page numbers only from the
textual section. LibreOffice may keep PAGE fields visible in the pre-textual
section even after the DOCX section break. This PDF-level patch removes only
standalone numeric blocks in the top-right header area before INTRODUÇÃO.
"""
import os
import pathlib
import re
import shutil

import fitz

PDF = pathlib.Path(os.environ.get("PDF", "TESE_CorePAM_FINAL_envio_Dr_Joao_impressao.pdf"))
TMP = PDF.with_suffix(".frontmatter-redact.tmp.pdf")

doc = fitz.open(PDF)

intro_page = None
intro_pat = re.compile(r"^\s*1\.\s+INTRODU[ÇC][ÃA]O\s*$", re.I)
for i, page in enumerate(doc):
    text = page.get_text("text") or ""
    for line in text.splitlines():
        if intro_pat.match(line.strip()):
            intro_page = i
            break
    if intro_page is not None:
        break

if intro_page is None:
    raise SystemExit("INTRODUÇÃO page not found; front-matter page numbers not redacted.")

count = 0
for i in range(intro_page):
    page = doc[i]
    blocks = page.get_text("blocks")
    for block in blocks:
        x0, y0, x1, y1, text = block[:5]
        if x0 >= 500 and y0 <= 75 and re.fullmatch(r"\s*\d+\s*", text or ""):
            rect = fitz.Rect(x0 - 2, y0 - 2, x1 + 2, y1 + 2)
            page.add_redact_annot(rect, fill=(1, 1, 1))
            count += 1
    if count:
        page.apply_redactions()

doc.save(TMP, garbage=4, deflate=True)
doc.close()
shutil.move(TMP, PDF)
print(f"Front-matter page numbers redacted before page {intro_page + 1}: {count}")
