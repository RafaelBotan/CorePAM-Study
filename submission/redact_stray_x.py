"""Redact stray standalone 'X' glyphs LibreOffice emits near bookmark anchors.

Run as the final pipeline step after the second LibreOffice PDF pass.
The 'X' spans originate from `<w:bookmarkStart/End>` elements that Pandoc
places around captioned tables/figures; fix_docx_postrender.py removes most,
but LibreOffice still renders a few single-glyph 'X' spans in the PDF.
"""
import fitz, os
PDF = os.environ.get('PDF', 'TESE_CorePAM_v0.2.pdf')
doc = fitz.open(PDF)
n = 0
for page in doc:
    tp = page.get_text('dict')
    for block in tp['blocks']:
        if block.get('type') != 0:
            continue
        for line in block['lines']:
            for span in line['spans']:
                t = span['text']
                x0, y0, x1, y1 = span['bbox']
                if t.strip() == 'X':
                    page.add_redact_annot(fitz.Rect(x0, y0, x1 + 1, y1), fill=(1, 1, 1))
                    n += 1
                elif len(t) > 2 and t.endswith('X') and t[-2] in '.,;:)0123456789 ':
                    char_w = (x1 - x0) / max(1, len(t))
                    page.add_redact_annot(
                        fitz.Rect(x1 - char_w * 1.1, y0, x1 + 1, y1), fill=(1, 1, 1)
                    )
                    n += 1
    page.apply_redactions()
doc.saveIncr()
print(f'Stray X glyphs redacted: {n}; pages: {len(doc)}')
