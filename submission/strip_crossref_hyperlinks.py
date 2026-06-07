# -*- coding: utf-8 -*-
# POS-RENDER OBRIGATORIO para a tese CorePAM.
#
# Problema: o Quarto renderiza cada cross-ref (@fig-x, @tbl-y) como um
# <w:hyperlink w:anchor="..."> interno no docx. O LibreOffice, ao exportar
# esse docx para PDF, imprime um "X" visivel grudado no fim do link
# (ex.: "...desenho do estudo.X"). Isso NAO esta no qmd nem no texto do docx —
# e artefato de renderizacao do LibreOffice sobre o hyperlink interno.
#
# Solucao: remover os <w:hyperlink w:anchor=...> (mantendo o texto interno) e
# o estilo de caractere Hyperlink, ANTES de exportar para PDF.
#
# Uso:
#   python strip_crossref_hyperlinks.py <entrada.docx> <saida.docx>
# Depois exportar saida.docx -> PDF com LibreOffice (perfil isolado, sem
# instancia aberta) conforme HANDOFF.
#
# IMPORTANTE: rodar este passo TODA VEZ que o qmd for re-renderizado para docx.

import sys, re, zipfile, shutil, os

src, dst = sys.argv[1], sys.argv[2]
shutil.copyfile(src, dst)

# Le document.xml de dentro do docx (zip)
with zipfile.ZipFile(dst, "r") as z:
    names = z.namelist()
    xml = z.read("word/document.xml").decode("utf-8")

# 1) Remove TODOS os <w:hyperlink> (internos w:anchor= e externos r:id=),
#    preservando o texto interno. O LibreOffice imprime um "X" visivel no fim
#    de QUALQUER hyperlink, inclusive os DOIs externos das referencias. Como as
#    referencias ABNT sao texto puro (DOI nao precisa ser clicavel), desembrulhar
#    o link e a correcao certa. O texto do DOI (inclusive os que terminam em "-X"
#    de verdade, como 10.1016/S1470-2045(17)30904-X) permanece intacto.
xml = re.sub(
    r'<w:hyperlink\b[^>]*>(.*?)</w:hyperlink>',
    r'\1', xml, flags=re.S)
# 2) Remove o estilo de caractere Hyperlink (azul/sublinhado residual)
xml = re.sub(r'<w:rStyle w:val="Hyperlink"\s*/>', '', xml)

# Reescreve o zip com o document.xml modificado
tmp = dst + ".tmp"
with zipfile.ZipFile(dst, "r") as zin, zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
    for item in zin.infolist():
        data = zin.read(item.filename)
        if item.filename == "word/document.xml":
            data = xml.encode("utf-8")
        zout.writestr(item, data)
os.replace(tmp, dst)

n = len(re.findall(r'w:anchor=', xml))
print(f"OK: {dst} | hyperlinks internos restantes (deve ser 0 de cross-ref): {n}")
