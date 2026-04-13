# Breast Cancer Research (BCR) - Formatting Guide

## Journal: Breast Cancer Research (BMC/Springer Nature)
- ISSN: 1465-542X
- Publisher: BioMed Central (BMC), Springer Nature

---

## Manuscript Format

### General
- **Double-line spacing** (use `referee` option in sn-jnl.cls)
- **Line numbering** (use `lineno` option)
- **Page numbering** required
- **SI units** required; embed all special characters
- **No page breaks** in manuscript
- **Single .tex file** (do not use `\input{...}` to include other tex files)
- All figures/files attached separately, not embedded in TeX source

### Accepted file formats
- Microsoft Word (DOC, DOCX)
- Rich Text Format (RTF)
- **TeX/LaTeX** (preferred for this project)
  - Compiled with **pdfLaTeX** and **TexLive 2021**
  - All editable source files must be uploaded

### LaTeX template
- Springer Nature template: `sn-jnl.cls` (Version 3.1, December 2024)
- Reference style for BCR: `sn-vancouver-num` (Vancouver numbered)
- Use `\unnumbered` for BCR (unnumbered section heads)
- documentclass: `\documentclass[referee,lineno,pdflatex,sn-vancouver-num]{sn-jnl}`

---

## Article Structure (Research Article)

1. **Title** (concise, informative)
2. **Abstract** (structured: Background, Methods, Results, Conclusions)
3. **Keywords** (3-10)
4. **Background** (= Introduction for BCR)
5. **Methods**
6. **Results**
7. **Discussion**
8. **Conclusions**
9. **List of abbreviations**
10. **Declarations** (Ethics, Consent, Data availability, Competing interests, Funding, Authors' contributions)
11. **Additional files** (supplementary)
12. **References**

---

## Figures

### In manuscript
- Use `\includegraphics[width=\textwidth]{filename.png}`
- Avoid subfigures; use multi-panel figures with A/B/C labels instead
- Each figure from a single input file
- For pdfLaTeX: use **PDF, JPG, or PNG** format

### Submission requirements
- **Resolution**: minimum 300 dpi (600 dpi preferred for line art)
- **Format**: TIFF, PNG, or PDF preferred; EPS also accepted
- **Maximum width**: ~170mm (single column) or ~210mm (full page)
- **Font**: consistent across all figures; minimum 8pt after scaling
- **Background**: white (no dark/screenshot backgrounds)
- **Color**: use colorblind-friendly palettes when possible

### Supplementary figures
- Numbered as "Additional file N: Figure SN"
- Same quality requirements as main figures
- Can be larger since they don't occupy main text space

---

## Tables

### Formatting
- Use `\begin{table}` with `booktabs` package (`\toprule`, `\midrule`, `\botrule`)
- Use `\begin{tabular*}{\textwidth}{@{\extracolsep\fill}...}` for full-width tables
- Footnotes inside tables via `\footnotetext[N]{...}`
- Keep tables concise; avoid redundancy
- For wide tables: use `\small` or `\footnotesize` font

### Supplementary tables
- Numbered as "Additional file N: Table SN"
- Can be provided as CSV for machine-readability (BCR encourages this)

---

## References

### Style
- **Vancouver numbered** (`sn-vancouver-num.bst`)
- Numbered in order of appearance in text
- Web links/URLs get reference numbers (not inline)
- Include DOI when available

### What can be cited
- Published articles, in-press articles, preprints
- NOT: unpublished data, personal communications (mention in text only)

### Format for data repositories
> "The dataset(s) supporting the conclusions of this article is(are) available in the [repository name] repository, [unique persistent identifier and hyperlink in https:// format]."

### Software citation format
- Project name, home page URL, archived version (DOI/Zenodo), OS, language, license

---

## Data and Code Availability

### Required section format
For this project:
```
All raw data are publicly available (GEO, GDC, cBioPortal).
Analysis code is available at https://github.com/RafaelBotan/CorePAM-Study.
All results can be reproduced from raw public data using the provided pipeline.
```

### Software details to include
- Project name: CorePAM-Study
- Project home page: https://github.com/RafaelBotan/CorePAM-Study
- Operating system: Platform independent (R)
- Programming language: R (version 4.5.2)
- License: (specify)
- Restrictions: None

---

## Additional Files (Supplementary)

- Referenced as "Additional file N" in main text
- Each with a brief title and description
- Figures: PNG/TIFF at 300+ dpi
- Tables: CSV preferred for machine-readability
- Data: deposit in public repositories when possible

---

## Checklist Before Submission

- [ ] Double-line spacing enabled
- [ ] Line and page numbering present
- [ ] Structured abstract (Background/Methods/Results/Conclusions)
- [ ] All figures at 300+ dpi, white background, consistent fonts
- [ ] All tables concise with booktabs formatting
- [ ] References in Vancouver numbered style
- [ ] Data availability statement present
- [ ] Ethics/declarations section complete
- [ ] All additional files numbered and described
- [ ] No page breaks in manuscript
- [ ] Single .tex file (no \input)
- [ ] .bib file included
- [ ] GitHub repo URL included and accessible
