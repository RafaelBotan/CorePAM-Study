"""Insert online publication status into reference 76.

The Vancouver CSL used by Pandoc does not render the BibTeX `note` field for
this article entry. This post-process inserts the CorePAM online publication
date in the `ref-botan2026` paragraph before the DOI hyperlink.
"""
import os
import re
import shutil
import zipfile
import xml.etree.ElementTree as ET


DOCX = os.environ.get("DOCX", "TESE_CorePAM_FINAL_envio_Dr_Joao_impressao.docx")
TMP = DOCX + ".ref76tmp"
OLD_STATUS = " Artigo aceito após revisão por pares; ainda não publicado."
STATUS = " Publicado online em 21 de maio de 2026."
STATUS_XML = (
    '<w:r><w:t xml:space="preserve"> Publicado online em 21 de maio de '
    "2026.</w:t></w:r>"
)

with zipfile.ZipFile(DOCX) as zin:
    entries = {name: zin.read(name) for name in zin.namelist()}

doc = entries["word/document.xml"].decode("utf-8")

start = doc.find('w:name="ref-botan2026"')
if start < 0:
    raise SystemExit("ref-botan2026 bookmark not found")
p_start = doc.rfind("<w:p", 0, start)
p_end = doc.find("</w:p>", start) + len("</w:p>")
para = doc[p_start:p_end]

if OLD_STATUS in re.sub(r"<[^>]+>", "", para):
    para = para.replace(
        '<w:r><w:t xml:space="preserve"> Artigo aceito após revisão por pares; ainda não publicado.</w:t></w:r>',
        "",
    )

if STATUS not in re.sub(r"<[^>]+>", "", para):
    target = '<w:t xml:space="preserve">-derived expression score with cross-platform external validation for breast cancer prognosis. Breast Cancer Res. 2026.</w:t></w:r>'
    if target not in para:
        raise SystemExit("Reference 76 journal/year run not found")
    para = para.replace(target, target + STATUS_XML, 1)
    doc = doc[:p_start] + para + doc[p_end:]
    print("Inserted online publication status into reference 76.")
else:
    print("Reference 76 already contains online publication status.")

ET.fromstring(doc)
entries["word/document.xml"] = doc.encode("utf-8")

with zipfile.ZipFile(TMP, "w", zipfile.ZIP_DEFLATED) as zout:
    for name, data in entries.items():
        zout.writestr(name, data)

shutil.move(TMP, DOCX)
print(f"Reference 76 publication status patched in {DOCX}")
