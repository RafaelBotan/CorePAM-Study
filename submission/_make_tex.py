"""Modify CorePAM_manuscript.tex:
- switch to clean publication mode (no referee/lineno)
- force figures/tables to appear inline (where mentioned), not as floats
- add a SEPARATE acceptance cover page before the title
"""
import re
from pathlib import Path

HERE = Path(__file__).parent
SRC = HERE / "CorePAM_manuscript.tex"
OUT = HERE / "CorePAM_FINAL_accepted.tex"

src = SRC.read_text(encoding="utf-8")

OLD = r"\documentclass[referee,lineno,pdflatex,sn-vancouver-num]{sn-jnl}"
NEW = r"\documentclass[pdflatex,sn-vancouver-num]{sn-jnl}"
assert OLD in src, "documentclass line not found"
src = src.replace(OLD, NEW)

# Force figures and tables to be placed near mention.
# sn-jnl's class wraps tables in threeparttable, which conflicts with float pkg [H].
# Use [!ht] (forced here-or-top) plus `placeins[section]` and tight float fractions
# to keep them very close to the text reference.
def _force_inline(match):
    return match.group(1) + "[!ht]"

src = re.sub(r"(\\begin\{(?:figure|table)\})\[[^\]]*\]", _force_inline, src)
src = re.sub(r"(\\begin\{(?:figure|table)\})(?!\[)", lambda m: m.group(1) + "[!ht]", src)

preamble_extra = r"""
%% ---- keep figures/tables near their text mention ----
\usepackage[section]{placeins}      % \FloatBarrier at every \section
\renewcommand{\topfraction}{0.95}
\renewcommand{\bottomfraction}{0.95}
\renewcommand{\textfraction}{0.05}
\renewcommand{\floatpagefraction}{0.85}
\setcounter{topnumber}{4}
\setcounter{bottomnumber}{4}
\setcounter{totalnumber}{8}

%% ---- COVER PAGE PACKAGES ----
\usepackage{xcolor}
\usepackage{etoolbox}
\definecolor{bcrblue}{HTML}{0B5394}
\definecolor{bcrlight}{HTML}{E8F0FA}
\newcommand{\acceptancecover}{%
  \thispagestyle{empty}%
  \begingroup
  \null\vspace*{2cm}
  \begin{center}
  {\color{bcrblue}\rule{\textwidth}{1.5pt}}\\[1.2em]
  {\color{bcrblue}\Huge\bfseries FINAL ACCEPTED MANUSCRIPT}\\[1em]
  {\color{bcrblue}\rule{\textwidth}{0.6pt}}\\[2em]

  {\Large Accepted for publication in}\\[0.6em]
  {\Huge\itshape Breast Cancer Research}\\[0.4em]
  {\large (BMC, Springer Nature)}\\[3em]

  \fcolorbox{bcrblue}{bcrlight}{%
    \begin{minipage}{0.85\textwidth}
    \vspace{0.6em}
    \begin{center}
    \textbf{Manuscript title}\\[0.3em]
    {\itshape CorePAM: a 24-gene PAM50-derived expression score with cross-platform external validation for breast cancer prognosis}\\[1.2em]

    \begin{tabular}{rl}
    \textbf{Acceptance date:} & 27 April 2026 \\[0.3em]
    \textbf{Submission ID:} & c2980ad5-9b3d-41a2-bf47-f093c5dc8c5f \\[0.3em]
    \textbf{Version:} & v3.0 (R2 revision) \\[0.3em]
    \textbf{Editorial status:} & Publishing and rights \\[0.3em]
    \textbf{DOI:} & to be assigned during typesetting
    \end{tabular}
    \vspace{0.6em}
    \end{center}
    \end{minipage}}\\[3em]

  \begin{minipage}{0.85\textwidth}
  \itshape\small This document is the author's final accepted manuscript, immediately following peer-review and editorial acceptance, prior to journal copy-editing and typesetting. The scientific content is identical to the version that will appear in the journal. Figures and tables are embedded near their first mention in the text, as in the published article.
  \end{minipage}

  \vfill
  {\color{bcrblue}\rule{\textwidth}{0.6pt}}\\[0.4em]
  {\small Rafael de Negreiros Botan \quad\textbullet\quad Universidade de Bras\'{i}lia \quad\textbullet\quad oncologista@gmail.com}
  \end{center}
  \endgroup
  \clearpage}
"""

anchor = r"\usepackage{url}"
assert anchor in src
src = src.replace(anchor, anchor + "\n" + preamble_extra)

# Insert cover page right after \begin{document}, before \maketitle
begin_doc = r"\begin{document}"
assert begin_doc in src
src = src.replace(begin_doc, begin_doc + "\n\n\\acceptancecover\n", 1)

OUT.write_text(src, encoding="utf-8")
print("wrote", OUT)
print("figures/tables forced to [H] (inline placement):",
      sum(1 for _ in re.finditer(r"\\begin\{(?:figure|table)\}\[H\]", src)))
