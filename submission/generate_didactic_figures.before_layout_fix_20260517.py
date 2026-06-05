"""
Gera figuras didáticas para TESE_CorePAM (Português, ABNT).

Estilo: editorial sofisticado (NYT / Economist / The Atlantic).
- Paleta warm-tinted: grafite, prussian muted, cream, terracotta muted
- Tipografia DejaVu Sans, respiração generosa
- Sem ornamentos "infográficos"; uso restrito de cor

Fase 1:
  Fig1 (A1) - Desenho geral do estudo
  B2  - Os 4 subtipos moleculares
  B3  - Timeline das assinaturas gênicas
  B6  - Fluxograma das coortes
  B9  - Pipeline de derivação do CorePAM
"""
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import (FancyBboxPatch, FancyArrowPatch,
                                Circle, Rectangle)

ROOT = Path(r"Y:/Doutorado Botan/CorePAM_Study_submitted")
OUT_PNG = ROOT / "figures" / "didaticas" / "pt" / "png"
OUT_PDF = ROOT / "figures" / "didaticas" / "pt" / "pdf"
MAIN_PNG = ROOT / "figures" / "main" / "pt" / "png"
MAIN_PDF = ROOT / "figures" / "main" / "pt" / "pdf"
for d in (OUT_PNG, OUT_PDF, MAIN_PNG, MAIN_PDF):
    d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Paleta warm editorial (inspirada em NYT / Economist / The Atlantic)
# Princípio: quase-preto grafite, prussian muted, cinzas warm,
# cream em vez de branco puro. Único acento = terracotta-burgundy muted
# reservado para destacar o objeto desta tese (CorePAM / pior prognóstico).
# ---------------------------------------------------------------------------
C = dict(
    ink="#1C1F24",           # grafite deep (texto primário, headers escuros)
    graphite="#2E3544",      # grafite escuro
    navy="#3E5573",          # prussian muted (primário)
    steel="#5C7795",          # editorial blue-grey (secundário)
    fog="#98A8BC",           # dusty blue light (terciário)
    sand="#D8D0C0",          # warm grey border
    cream="#F5F1EA",         # warm off-white (bg agrupamento)
    paper="#FBF8F3",         # warm near-white (bg muito sutil)
    stone="#78726A",         # warm grey (captions, subtítulos)
    grey_mid="#5A5F6A",      # medium grey
    grey_soft="#A5A09A",     # light warm grey
    white="#FFFFFF",
    accent="#7D3E3A",        # terracotta-burgundy muted
    accent_light="#A56657",
)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.edgecolor": C["grey_mid"],
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})


def save(fig, name, main=False):
    png_dir = MAIN_PNG if main else OUT_PNG
    pdf_dir = MAIN_PDF if main else OUT_PDF
    fig.savefig(png_dir / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(pdf_dir / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  OK  {name}")


def rbox(ax, x, y, w, h, text, fc, ec=None, fontsize=10, fontcolor="white",
         weight="bold", zorder=3, ha="center", va="center", pad=0.03,
         rounding=0.05, lw=1.0):
    ec = ec or fc
    box = FancyBboxPatch((x, y), w, h,
                         boxstyle=f"round,pad={pad},rounding_size={rounding}",
                         fc=fc, ec=ec, lw=lw, zorder=zorder)
    ax.add_patch(box)
    ax.text(x + w/2, y + h/2, text, ha=ha, va=va,
            fontsize=fontsize, color=fontcolor, weight=weight, zorder=zorder+1)


def arrow(ax, x1, y1, x2, y2, color=None, lw=1.4, style="->"):
    color = color or C["grey_mid"]
    ar = FancyArrowPatch((x1, y1), (x2, y2),
                         arrowstyle=style, mutation_scale=12,
                         color=color, lw=lw, zorder=2)
    ax.add_patch(ar)


# ============================================================================
# Fig1 (A1) — Desenho geral do estudo
# ============================================================================
def fig_A1_desenho_estudo():
    fig, ax = plt.subplots(figsize=(13, 9.5))
    ax.set_xlim(0, 13); ax.set_ylim(0, 10)
    ax.axis("off")

    ax.text(6.5, 9.60, "Desenho do estudo",
            ha="center", fontsize=17, weight="bold", color=C["ink"])
    ax.text(6.5, 9.22,
            "Derivação em coorte populacional, validação externa e "
            "análise secundária de resposta à quimioterapia neoadjuvante",
            ha="center", fontsize=10, style="italic", color=C["stone"])

    # ===== Bloco 1 — Derivação =====
    rbox(ax, 0.5, 5.7, 4.5, 2.8, "", fc=C["paper"],
         ec=C["sand"], pad=0.02, rounding=0.06, zorder=1, lw=1.0)
    ax.text(2.75, 8.22, "DERIVAÇÃO", ha="center",
            fontsize=11, weight="bold", color=C["navy"])
    ax.plot([0.85, 4.65], [8.03, 8.03], color=C["navy"], lw=1.0, zorder=2)

    rbox(ax, 0.9, 7.05, 3.7, 0.80,
         "SCAN-B (GSE96058)\nN = 3.069 · RNA-seq · 322 eventos (OS)",
         fc=C["navy"], fontsize=9.5, pad=0.03, lw=0)
    rbox(ax, 0.9, 5.90, 1.7, 0.95, "PAM50\n(50 genes)",
         fc=C["white"], ec=C["steel"], fontcolor=C["ink"],
         fontsize=10, weight="bold", pad=0.03, lw=1.3)
    arrow(ax, 2.6, 6.38, 2.95, 6.38, color=C["steel"], lw=1.5)
    rbox(ax, 2.95, 5.90, 1.7, 0.95, "CorePAM\n(24 genes)",
         fc=C["accent"], fontsize=10, pad=0.03, lw=0)

    # ===== Bloco 2 — Validação externa =====
    rbox(ax, 5.5, 5.7, 7.0, 2.8, "", fc=C["paper"],
         ec=C["sand"], pad=0.02, rounding=0.06, zorder=1, lw=1.0)
    ax.text(9.0, 8.22, "VALIDAÇÃO EXTERNA · SOBREVIDA", ha="center",
            fontsize=11, weight="bold", color=C["navy"])
    ax.plot([5.85, 12.15], [8.03, 8.03], color=C["navy"], lw=1.0, zorder=2)

    val = [("TCGA-BRCA", "RNA-seq", "N = 1.072", "150 eventos · OS"),
           ("METABRIC",  "Microarray", "N = 1.978", "646 eventos · DSS"),
           ("GSE20685",  "Microarray", "N = 327",   "83 eventos · OS")]
    cx = [5.8, 8.0, 10.2]; cw = 2.2
    for (name, plat, n, ev), x in zip(val, cx):
        rbox(ax, x, 7.05, cw, 0.80, f"{name}\n{plat}",
             fc=C["steel"], fontsize=9.5, pad=0.03, lw=0)
        rbox(ax, x, 5.90, cw, 0.95, f"{n}\n{ev}",
             fc=C["white"], ec=C["steel"], fontcolor=C["ink"],
             fontsize=9.0, weight="normal", pad=0.03, lw=1.0)
    ax.text(9.0, 5.75,
            "GSE1456 · N = 159 · 40 eventos (OS) — sensibilidade",
            ha="center", fontsize=8.8, style="italic", color=C["stone"])

    arrow(ax, 5.0, 6.65, 5.5, 6.65, color=C["graphite"], lw=1.8)
    ax.text(5.25, 6.88, "modelo\ncongelado", ha="center",
            fontsize=7.5, style="italic", color=C["stone"])

    # ===== Parâmetros analíticos (2 linhas para caber) =====
    rbox(ax, 0.5, 4.10, 12.0, 1.10,
         "Parâmetros congelados: elastic-net (α = 0,5)  ·  K = 10 folds (SHA-256 do ID)\n"
         "grade λ log, 100 valores  ·  ΔC ≤ 0,010  ·  z-score intra-coorte por gene  ·  tempo em meses",
         fc=C["cream"], ec=C["sand"],
         fontcolor=C["ink"], fontsize=9.5, weight="normal", pad=0.03, lw=1.0)

    # ===== Bloco 3 — Análise secundária (pCR) =====
    rbox(ax, 0.5, 1.15, 12.0, 2.95, "", fc=C["paper"],
         ec=C["sand"], pad=0.02, rounding=0.06, zorder=1, lw=1.0)
    ax.text(6.5, 3.83,
            "ANÁLISE SECUNDÁRIA · RESPOSTA À QUIMIOTERAPIA (pCR)",
            ha="center", fontsize=11, weight="bold", color=C["accent"])
    ax.plot([0.85, 12.15], [3.64, 3.64], color=C["accent"], lw=1.0, zorder=2)

    pcr = [("GSE25066", "N = 182", "pCR 23,1 %"),
           ("GSE20194", "N = 278", "pCR 20,1 %"),
           ("GSE32646", "N = 115", "pCR 23,5 %"),
           ("I-SPY1",   "N = 122", "pCR 26,2 %"),
           ("I-SPY2",   "N = 986", "pCR 32,4 %")]
    px = [0.8, 3.2, 5.6, 8.0, 10.4]; pw = 2.2
    for (name, n, p), x in zip(pcr, px):
        rbox(ax, x, 2.60, pw, 0.75, name,
             fc=C["graphite"], fontsize=9.5, pad=0.03, lw=0)
        rbox(ax, x, 1.45, pw, 1.00, f"{n}\n{p}",
             fc=C["white"], ec=C["graphite"], fontcolor=C["ink"],
             fontsize=9.0, weight="normal", pad=0.03, lw=1.0)

    arrow(ax, 3.8, 5.8, 3.8, 4.15, color=C["accent"], lw=1.8)
    ax.text(4.0, 5.0, "aplicação do\nescore congelado",
            fontsize=7.5, style="italic", color=C["accent"], ha="left")

    ax.text(6.5, 0.70,
            "OS: overall survival; DSS: disease-specific survival; "
            "pCR: resposta patológica completa; ΔC: diferença de C-index.",
            ha="center", fontsize=8.5, style="italic", color=C["stone"])

    save(fig, "Fig1_StudyDesign_PT", main=True)


# ============================================================================
# B2 — Os 4 subtipos moleculares (footer FORA dos cards)
# ============================================================================
def fig_B2_subtipos_moleculares():
    # Figura alta com espaço reservado para rodapé abaixo dos cards
    fig, ax = plt.subplots(figsize=(12, 10.5))
    ax.set_xlim(0, 12); ax.set_ylim(0, 10.5)
    ax.axis("off")

    ax.text(6, 10.05, "Os quatro subtipos moleculares do câncer de mama",
            ha="center", fontsize=16, weight="bold", color=C["ink"])
    ax.text(6, 9.65,
            "Classificação intrínseca (PAM50) — prevalência, marcadores e prognóstico",
            ha="center", fontsize=10, style="italic", color=C["stone"])

    # Cards: topo cy=7.0 (spans 5.4–8.6), base cy=3.0 (spans 1.4–4.6)
    # Rodapé em y=0.6 (bem abaixo dos cards) para não cortar texto
    subtypes = [
        dict(cx=3.0, cy=7.0, name="LUMINAL A", pct=40, color=C["navy"],
             markers=["RE +", "RP +", "HER2 −", "Ki67 baixo"],
             prog="Excelente prognóstico",
             chemo="Baixa dependência de quimioterapia"),
        dict(cx=9.0, cy=7.0, name="LUMINAL B", pct=20, color=C["steel"],
             markers=["RE +", "RP +/−", "HER2 −/+", "Ki67 alto"],
             prog="Prognóstico intermediário",
             chemo="Benefício variável de quimioterapia"),
        dict(cx=3.0, cy=3.0, name="HER2-ENRIQUECIDO", pct=15, color=C["graphite"],
             markers=["RE −", "RP −", "HER2 +++", "Ki67 alto"],
             prog="Ruim sem terapia anti-HER2",
             chemo="Resposta a trastuzumabe"),
        dict(cx=9.0, cy=3.0, name="TRIPLO-NEGATIVO", pct=15, color=C["accent"],
             markers=["RE −", "RP −", "HER2 −", "CK5/6 +"],
             prog="Pior prognóstico",
             chemo="Alta resposta à quimioterapia"),
    ]

    for s in subtypes:
        cx, cy = s["cx"], s["cy"]
        w, h = 5.2, 3.2
        # Card base (fundo paper, borda sutil)
        card = FancyBboxPatch((cx - w/2, cy - h/2), w, h,
                              boxstyle="round,pad=0.03,rounding_size=0.08",
                              fc=C["paper"], ec=C["sand"], lw=1.2, zorder=2)
        ax.add_patch(card)

        # Fita lateral colorida (editorial accent bar)
        ribbon = Rectangle((cx - w/2 + 0.05, cy - h/2 + 0.05),
                            0.22, h - 0.10, fc=s["color"], ec="none",
                            zorder=3)
        ax.add_patch(ribbon)

        # Círculo de prevalência
        circ_cx = cx - w/2 + 1.30
        ax.add_patch(Circle((circ_cx, cy), 0.82,
                            fc=s["color"], ec="none", zorder=4))
        ax.text(circ_cx, cy, f"{s['pct']} %",
                ha="center", va="center", fontsize=19,
                weight="bold", color=C["white"], zorder=5)

        # Nome do subtipo (topo direita do card)
        text_x = cx - w/2 + 2.55
        ax.text(text_x, cy + h/2 - 0.48, s["name"],
                ha="left", va="center", fontsize=12,
                weight="bold", color=s["color"], zorder=5)
        # Fina linha sob o nome
        ax.plot([text_x, text_x + 2.0], [cy + h/2 - 0.72, cy + h/2 - 0.72],
                color=s["color"], lw=0.8, zorder=5)

        # Marcadores (2 colunas compactas)
        mk_y = cy + 0.32
        ax.text(text_x, mk_y,
                f"{s['markers'][0]}    {s['markers'][1]}",
                ha="left", fontsize=10.5, color=C["ink"])
        ax.text(text_x, mk_y - 0.40,
                f"{s['markers'][2]}    {s['markers'][3]}",
                ha="left", fontsize=10.5, color=C["ink"])

        # Prognóstico e QT (base do card)
        ax.text(text_x, cy - h/2 + 0.62, s["prog"],
                ha="left", fontsize=10, weight="bold", color=s["color"])
        ax.text(text_x, cy - h/2 + 0.32, s["chemo"],
                ha="left", fontsize=9, color=C["stone"], style="italic")

    # Rodapé BEM ABAIXO dos cards (y=0.6, cards terminam em y=1.4)
    ax.text(6, 0.80,
            "RE: receptor de estrogênio  ·  RP: receptor de progesterona  ·  "
            "HER2: receptor tipo 2 do fator de crescimento epidérmico humano",
            ha="center", fontsize=8.5, style="italic", color=C["stone"])
    ax.text(6, 0.50,
            "Ki67: índice proliferativo  ·  CK5/6: citoqueratina 5/6  ·  "
            "QT: quimioterapia",
            ha="center", fontsize=8.5, style="italic", color=C["stone"])

    save(fig, "B2_Subtipos_Moleculares_PT")


# ============================================================================
# B3 — Timeline (header EMPILHADO: nome em cima, N genes embaixo)
# ============================================================================
def fig_B3_timeline_assinaturas():
    # Cada assinatura é um card unificado com 3 zonas internas:
    #   Header colorido (2 linhas empilhadas): nome + "N genes"
    #   Meio: plataforma em italic
    #   Base: propósito em cinza
    fig, ax = plt.subplots(figsize=(15, 10))
    ax.set_xlim(1999, 2029); ax.set_ylim(-6.2, 6.2)
    ax.axis("off")

    ax.text(2014, 5.75,
            "Evolução histórica das assinaturas gênicas prognósticas em câncer de mama",
            ha="center", fontsize=15, weight="bold", color=C["ink"])
    ax.text(2014, 5.25,
            "Redução progressiva no número de genes e migração para plataformas de leitura clínica",
            ha="center", fontsize=10, style="italic", color=C["stone"])

    # Linha do tempo
    ax.plot([2001, 2027], [0.4, 0.4], color=C["graphite"], lw=1.8,
            solid_capstyle="round", zorder=1)
    for yr in [2002, 2006, 2010, 2014, 2018, 2022, 2026]:
        ax.plot([yr, yr], [0.28, 0.52], color=C["graphite"], lw=1.0, zorder=1)
        ax.text(yr, -0.15, str(yr), ha="center", fontsize=10,
                color=C["stone"], weight="bold")

    sigs = [
        (2002, "MammaPrint",  70, "Microarray Agilent",
         "Recidiva em 5 anos (RE+/RE−)",         C["fog"],    "top"),
        (2004, "OncotypeDX",  21, "RT-PCR (FFPE)",
         "Benefício de QT em RE+",                 C["steel"],  "bot"),
        (2009, "PAM50",       50, "Microarray / RT-PCR",
         "Subtipos intrínsecos + prognóstico",     C["navy"],   "top"),
        (2013, "EndoPredict", 12, "RT-PCR (FFPE)",
         "Recidiva tardia em RE+",                 C["steel"],  "bot"),
        (2026, "CorePAM",     24, "RNA-seq + microarray",
         "Escore derivado do PAM50",               C["accent"], "top"),
    ]

    card_w = 4.6
    card_h = 2.3
    stem_len = 1.3
    hdr_h = 1.0  # maior, pra acomodar 2 linhas empilhadas

    for yr, name, n, plat, purpose, col, side in sigs:
        if side == "top":
            stem_start_y = 0.52
            stem_end_y = stem_start_y + stem_len
            card_y = stem_end_y
        else:
            stem_start_y = 0.28
            stem_end_y = stem_start_y - stem_len
            card_y = stem_end_y - card_h

        # Haste
        ax.plot([yr, yr], [stem_start_y, stem_end_y],
                color=col, lw=1.8, zorder=2)
        # Marcador na linha
        ax.add_patch(Circle((yr, 0.40), 0.20, fc=col, ec=C["white"],
                             lw=1.6, zorder=5))

        # Card base
        card_x = yr - card_w / 2
        card = FancyBboxPatch((card_x, card_y), card_w, card_h,
                              boxstyle="round,pad=0.02,rounding_size=0.08",
                              fc=C["white"], ec=col, lw=1.6, zorder=3)
        ax.add_patch(card)

        # Header colorido (com cantos arredondados no TOPO via rounding)
        hdr = FancyBboxPatch((card_x, card_y + card_h - hdr_h),
                             card_w, hdr_h,
                             boxstyle="round,pad=0.0,rounding_size=0.08",
                             fc=col, ec=col, zorder=4)
        ax.add_patch(hdr)
        # EMPILHADO: linha 1 = nome, linha 2 = N genes
        ax.text(card_x + card_w/2, card_y + card_h - hdr_h/2 + 0.15,
                name, ha="center", va="center",
                fontsize=13, weight="bold", color=C["white"], zorder=5)
        ax.text(card_x + card_w/2, card_y + card_h - hdr_h/2 - 0.25,
                f"{n} genes", ha="center", va="center",
                fontsize=11, weight="normal", color=C["white"], zorder=5)

        # Plataforma (italic no meio)
        ax.text(card_x + card_w/2, card_y + card_h - hdr_h - 0.32,
                plat, ha="center", va="center", fontsize=10.5,
                color=C["ink"], style="italic", zorder=5)

        # Divisor sutil
        ax.plot([card_x + 0.30, card_x + card_w - 0.30],
                [card_y + 0.38, card_y + 0.38],
                color=C["sand"], lw=0.8, zorder=5)

        # Propósito (base)
        ax.text(card_x + card_w/2, card_y + 0.22,
                purpose, ha="center", va="center", fontsize=9.5,
                color=C["grey_mid"], zorder=5)

    # Anel de destaque CorePAM
    ax.add_patch(Circle((2026, 0.40), 0.38, fc="none", ec=C["accent"],
                        lw=1.8, ls="--", zorder=4))

    ax.text(2014, -5.95,
            "Fonte: elaborada a partir de revisão narrativa da literatura. "
            "A assinatura CorePAM é resultado desta tese.",
            ha="center", fontsize=8.8, style="italic", color=C["stone"])

    save(fig, "B3_Timeline_Assinaturas_PT")


# ============================================================================
# B6 — Fluxograma das coortes
# ============================================================================
def fig_B6_consort_coortes():
    fig, ax = plt.subplots(figsize=(13, 11))
    ax.set_xlim(0, 13); ax.set_ylim(0, 11)
    ax.axis("off")

    ax.text(6.5, 10.65, "Fluxograma das coortes do estudo",
            ha="center", fontsize=16, weight="bold", color=C["ink"])
    ax.text(6.5, 10.25,
            "9 coortes públicas — 5 para validação prognóstica (sobrevida) "
            "e 4 + 1 para análise secundária (pCR)",
            ha="center", fontsize=10, style="italic", color=C["stone"])

    # Derivação (topo, central)
    rbox(ax, 5.2, 8.95, 2.6, 0.60, "COORTE DE DERIVAÇÃO",
         fc=C["navy"], fontsize=10.5, lw=0)
    rbox(ax, 5.0, 8.10, 3.0, 0.75,
         "SCAN-B (GSE96058)\nN = 3.069 · RNA-seq",
         fc=C["white"], ec=C["navy"], fontcolor=C["ink"],
         fontsize=10, weight="normal", lw=1.3)

    # Nó central Score CorePAM
    ax.add_patch(Circle((6.5, 7.15), 0.55, fc=C["accent"], ec=C["white"],
                        lw=1.8, zorder=4))
    ax.text(6.5, 7.15, "Score\nCorePAM", ha="center", va="center",
            fontsize=9.5, weight="bold", color=C["white"], zorder=5)

    arrow(ax, 6.5, 8.05, 6.5, 7.72, color=C["navy"], lw=1.8)

    # Bifurcação
    arrow(ax, 6.2, 7.00, 3.0, 6.20, color=C["navy"], lw=1.6)
    arrow(ax, 6.8, 7.00, 10.0, 6.20, color=C["accent"], lw=1.6)

    # === Sobrevida (esquerda) ===
    panel_surv = FancyBboxPatch((0.3, 1.2), 5.9, 4.9,
                                boxstyle="round,pad=0.03,rounding_size=0.08",
                                fc=C["paper"], ec=C["sand"],
                                lw=1.0, zorder=1)
    ax.add_patch(panel_surv)
    ax.text(3.25, 5.80, "VALIDAÇÃO DE SOBREVIDA", ha="center",
            fontsize=11.5, weight="bold", color=C["navy"])
    ax.text(3.25, 5.50, "5 coortes — N = 6.605 pacientes",
            ha="center", fontsize=9.5, style="italic", color=C["stone"])

    surv = [
        ("SCAN-B *",  "N = 3.069", "OS · 322 eventos", 0.5, 4.30),
        ("TCGA-BRCA", "N = 1.072", "OS · 150 eventos", 3.20, 4.30),
        ("METABRIC",  "N = 1.978", "DSS · 646 eventos", 0.5, 2.75),
        ("GSE20685",  "N = 327",   "OS · 83 eventos",   3.20, 2.75),
        ("GSE1456",   "N = 159",   "OS · 40 eventos",   1.85, 1.35),
    ]
    for name, n, ev, x, y in surv:
        w, h = 2.6, 1.05
        rbox(ax, x, y, w, h, "", fc=C["white"], ec=C["steel"], pad=0.02,
             rounding=0.06, lw=1.1)
        ax.add_patch(FancyBboxPatch((x + 0.08, y + h - 0.40),
                                    w - 0.16, 0.35,
                                    boxstyle="round,pad=0.0,rounding_size=0.12",
                                    fc=C["steel"], ec=C["steel"], zorder=4))
        ax.text(x + w/2, y + h - 0.225, name, ha="center", va="center",
                fontsize=10, weight="bold", color=C["white"], zorder=5)
        ax.text(x + w/2, y + 0.45, n, ha="center", fontsize=10,
                color=C["ink"], weight="bold", zorder=5)
        ax.text(x + w/2, y + 0.18, ev, ha="center", fontsize=8.8,
                color=C["stone"], style="italic", zorder=5)

    # === pCR (direita) ===
    panel_pcr = FancyBboxPatch((6.8, 1.2), 5.9, 4.9,
                               boxstyle="round,pad=0.03,rounding_size=0.08",
                               fc=C["paper"], ec=C["sand"],
                               lw=1.0, zorder=1)
    ax.add_patch(panel_pcr)
    ax.text(9.75, 5.80, "ANÁLISE DE pCR", ha="center",
            fontsize=11.5, weight="bold", color=C["accent"])
    ax.text(9.75, 5.50, "4 + 1 coortes — N = 1.683 pacientes",
            ha="center", fontsize=9.5, style="italic", color=C["stone"])

    pcr = [
        ("GSE25066",     "N = 182", "pCR = 23,1 %", 7.00, 4.30),
        ("GSE20194",     "N = 278", "pCR = 20,1 %", 9.70, 4.30),
        ("GSE32646",     "N = 115", "pCR = 23,5 %", 7.00, 2.75),
        ("I-SPY1",       "N = 122", "pCR = 26,2 %", 9.70, 2.75),
        ("I-SPY2 (exp.)","N = 986", "pCR = 32,4 %", 8.35, 1.35),
    ]
    for name, n, p, x, y in pcr:
        w, h = 2.6, 1.05
        rbox(ax, x, y, w, h, "", fc=C["white"], ec=C["graphite"], pad=0.02,
             rounding=0.06, lw=1.1)
        ax.add_patch(FancyBboxPatch((x + 0.08, y + h - 0.40),
                                    w - 0.16, 0.35,
                                    boxstyle="round,pad=0.0,rounding_size=0.12",
                                    fc=C["graphite"], ec=C["graphite"], zorder=4))
        ax.text(x + w/2, y + h - 0.225, name, ha="center", va="center",
                fontsize=10, weight="bold", color=C["white"], zorder=5)
        ax.text(x + w/2, y + 0.45, n, ha="center", fontsize=10,
                color=C["ink"], weight="bold", zorder=5)
        ax.text(x + w/2, y + 0.18, p, ha="center", fontsize=8.8,
                color=C["stone"], style="italic", zorder=5)

    ax.text(6.5, 0.65,
            "TOTAL GLOBAL: 10 coortes públicas · N = 8.288 pacientes",
            ha="center", fontsize=12, weight="bold", color=C["ink"])
    ax.text(6.5, 0.25,
            "* Coorte de treinamento (não usada como validação independente). "
            "OS: overall survival; DSS: disease-specific survival; "
            "pCR: resposta patológica completa.",
            ha="center", fontsize=8.5, style="italic", color=C["stone"])

    save(fig, "B6_Fluxograma_Coortes_PT")


# ============================================================================
# B9 — Pipeline de derivação (SEM tags)
# ============================================================================
def fig_B9_derivacao_corepam():
    fig, ax = plt.subplots(figsize=(12, 11))
    ax.set_xlim(0, 12); ax.set_ylim(0, 12)
    ax.axis("off")

    ax.text(6, 11.55, "Derivação do CorePAM — pipeline reprodutível",
            ha="center", fontsize=16, weight="bold", color=C["ink"])
    ax.text(6, 11.15,
            "Do painel PAM50 (50 genes candidatos) ao CorePAM (24 genes selecionados)",
            ha="center", fontsize=10, style="italic", color=C["stone"])

    # 6 steps — limpo, sem tags laterais
    steps = [
        ("ENTRADA · SCAN-B (N = 3.069)\nExpressão RNA-seq  ·  50 genes PAM50 candidatos",
         C["navy"],     9.80, 1.0),
        ("PRÉ-PROCESSAMENTO\nLog2(FPKM + 1)  →  z-score intra-coorte por gene",
         C["steel"],    8.20, 1.0),
        ("ATRIBUIÇÃO DETERMINÍSTICA DE FOLDS\nSHA-256 do ID do paciente  →  10 folds · sem random seed",
         C["steel"],    6.60, 1.0),
        ("VALIDAÇÃO CRUZADA 10-FOLD\nCox elastic-net (α = 0,5)  ·  grade λ log (100 valores)",
         C["navy"],     5.00, 1.0),
        ("FRONTEIRA DE PARETO\nC-index out-of-fold × número de genes  ·  critério ΔC ≤ 0,010",
         C["graphite"], 3.40, 1.0),
        ("SAÍDA · CorePAM — 24 genes selecionados\nC-index OOF = 0,670 (máximo 50 = 0,679)  ·  gap = 0,009",
         C["accent"],   1.80, 1.0),
    ]

    box_x, box_w = 1.0, 10.0
    for txt, col, y, h in steps:
        rbox(ax, box_x, y, box_w, h, txt, fc=col, fontsize=11, pad=0.04,
             rounding=0.06, lw=0)

    # Setas entre steps — centradas no gap
    for i in range(len(steps) - 1):
        y_top = steps[i][2]
        y_bot = steps[i+1][2] + steps[i+1][3]
        mid_x = box_x + box_w / 2
        arrow(ax, mid_x, y_top - 0.02,
              mid_x, y_bot + 0.02,
              color=C["grey_mid"], lw=2)

    # Rodapé
    ax.text(6, 0.90,
            "OOF: out-of-fold (predição fora da dobra)  ·  "
            "elastic-net: combinação de regularização L1 (lasso) e L2 (ridge)  ·  "
            "ΔC: diferença de C-index.",
            ha="center", fontsize=8.8, style="italic", color=C["stone"])
    ax.text(6, 0.55,
            "Fonte: elaborada pelo autor.",
            ha="center", fontsize=8.5, style="italic", color=C["stone"])

    save(fig, "B9_Derivacao_CorePAM_PT")


if __name__ == "__main__":
    print("Gerando figuras didáticas (paleta warm editorial)")
    fig_A1_desenho_estudo()
    fig_B2_subtipos_moleculares()
    fig_B3_timeline_assinaturas()
    fig_B6_consort_coortes()
    fig_B9_derivacao_corepam()
    print("\nConcluído.")
    print(f"  PNG didáticas: {OUT_PNG}")
    print(f"  PDF didáticas: {OUT_PDF}")
    print(f"  PNG principal: {MAIN_PNG}")
