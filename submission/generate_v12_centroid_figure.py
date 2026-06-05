"""Generate a didactic v1.2 figure explaining PAM50 centroids and CorePAM.

This script writes a new v1.2-only asset. It does not overwrite the thesis
v1.1 figures.
"""
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle


ROOT = Path(__file__).resolve().parents[1]
OUT_PNG = ROOT / "figures" / "didaticas" / "pt" / "png"
OUT_PDF = ROOT / "figures" / "didaticas" / "pt" / "pdf"
OUT_PNG.mkdir(parents=True, exist_ok=True)
OUT_PDF.mkdir(parents=True, exist_ok=True)

C = {
    "ink": "#1C1F24",
    "graphite": "#2E3544",
    "navy": "#3E5573",
    "steel": "#5C7795",
    "fog": "#98A8BC",
    "sand": "#D8D0C0",
    "cream": "#F5F1EA",
    "paper": "#FBF8F3",
    "stone": "#78726A",
    "white": "#FFFFFF",
    "accent": "#7D3E3A",
    "green": "#557A5D",
}

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})


def rounded(ax, xy, wh, fc, ec=None, lw=1.0, radius=0.04, z=2):
    ec = ec or fc
    patch = FancyBboxPatch(
        xy, wh[0], wh[1],
        boxstyle=f"round,pad=0.035,rounding_size={radius}",
        fc=fc, ec=ec, lw=lw, zorder=z
    )
    ax.add_patch(patch)
    return patch


def arrow(ax, x1, y1, x2, y2, color=None, lw=1.4):
    ax.add_patch(FancyArrowPatch(
        (x1, y1), (x2, y2), arrowstyle="->", mutation_scale=15,
        color=color or C["stone"], lw=lw, zorder=4
    ))


def gene_bar(ax, x, y, w, label, value, color):
    ax.text(x, y + 0.06, label, ha="left", va="bottom",
            fontsize=8.0, color=C["ink"], zorder=6)
    ax.add_patch(Rectangle((x + 0.55, y), 1.2, 0.10,
                           fc=C["cream"], ec=C["sand"], lw=0.4, zorder=5))
    ax.add_patch(Rectangle((x + 0.55, y), 1.2 * value, 0.10,
                           fc=color, ec="none", zorder=6))


def centroid_card(ax, x, y, title, corr, color, selected=False):
    fc = C["paper"] if not selected else "#F9EFEA"
    ec = color if selected else C["sand"]
    lw = 1.6 if selected else 0.8
    rounded(ax, (x, y), (2.15, 0.72), fc=fc, ec=ec, lw=lw, radius=0.035, z=3)
    ax.text(x + 0.18, y + 0.48, title, ha="left", va="center",
            fontsize=9.2, weight="bold", color=color, zorder=5)
    ax.text(x + 0.18, y + 0.23, f"correlação: {corr}",
            ha="left", va="center", fontsize=7.7, color=C["stone"], zorder=5)
    if selected:
        ax.text(x + 1.98, y + 0.49, "maior", ha="right", va="center",
                fontsize=7.5, color=color, weight="bold", zorder=6)


def main():
    fig, ax = plt.subplots(figsize=(14.2, 8.4))
    ax.set_xlim(0, 14.2)
    ax.set_ylim(0, 8.4)
    ax.axis("off")

    ax.text(7.1, 8.05, "PAM50: centroides moleculares e a pergunta do CorePAM",
            ha="center", fontsize=16.5, weight="bold", color=C["ink"])
    ax.text(7.1, 7.68,
            "O PAM50 compara a amostra com perfis de referência; o CorePAM pergunta se o sinal prognóstico pode ser comprimido.",
            ha="center", fontsize=9.8, style="italic", color=C["stone"])

    # Left panel: patient expression profile
    rounded(ax, (0.45, 1.0), (3.55, 6.1), fc=C["paper"], ec=C["sand"], lw=1.0, radius=0.06)
    ax.text(2.22, 6.72, "1. Tumor da paciente", ha="center",
            fontsize=11.2, weight="bold", color=C["navy"])
    ax.text(2.22, 6.33, "expressão dos genes PAM50", ha="center",
            fontsize=8.7, style="italic", color=C["stone"])

    genes = [
        ("ESR1", 0.88, C["navy"]),
        ("PGR", 0.75, C["navy"]),
        ("ERBB2", 0.32, C["steel"]),
        ("MKI67", 0.58, C["accent"]),
        ("KRT5", 0.18, C["graphite"]),
        ("FOXC1", 0.22, C["graphite"]),
        ("BCL2", 0.70, C["green"]),
        ("UBE2C", 0.52, C["accent"]),
    ]
    yy = 5.65
    for lab, val, col in genes:
        gene_bar(ax, 1.0, yy, 1.2, lab, val, col)
        yy -= 0.45

    rounded(ax, (0.9, 1.45), (2.65, 1.02), fc=C["white"], ec=C["sand"], lw=0.8, radius=0.03)
    ax.text(2.22, 2.13, "Vetor de expressão", ha="center",
            fontsize=9.0, weight="bold", color=C["ink"])
    ax.text(2.22, 1.82, "um padrão multigênico,\nnão uma única lâmina de IHQ",
            ha="center", fontsize=8.1, color=C["stone"], linespacing=1.2)

    # Middle panel: centroids
    rounded(ax, (5.05, 1.0), (4.15, 6.1), fc=C["paper"], ec=C["sand"], lw=1.0, radius=0.06)
    ax.text(7.13, 6.72, "2. Comparação com centroides", ha="center",
            fontsize=11.2, weight="bold", color=C["navy"])
    ax.text(7.13, 6.34, "centroide = perfil de referência do subtipo", ha="center",
            fontsize=8.7, style="italic", color=C["stone"])

    centroid_card(ax, 6.02, 5.40, "Luminal A", "0,82", C["navy"], selected=True)
    centroid_card(ax, 6.02, 4.48, "Luminal B", "0,63", C["steel"])
    centroid_card(ax, 6.02, 3.56, "HER2-enriched", "0,21", C["graphite"])
    centroid_card(ax, 6.02, 2.64, "Basal-like", "0,08", C["accent"])
    centroid_card(ax, 6.02, 1.72, "Normal-like", "0,12", C["fog"])

    rounded(ax, (5.42, 0.78), (3.38, 0.55), fc=C["white"], ec=C["navy"], lw=1.1, radius=0.03)
    ax.text(7.11, 1.05, "subtipo = centroide mais parecido", ha="center",
            fontsize=8.7, weight="bold", color=C["navy"])

    # Right panel: CorePAM compression
    rounded(ax, (10.25, 1.0), (3.5, 6.1), fc=C["paper"], ec=C["sand"], lw=1.0, radius=0.06)
    ax.text(12.0, 6.72, "3. Pergunta do CorePAM", ha="center",
            fontsize=11.2, weight="bold", color=C["accent"])
    ax.text(12.0, 6.34, "compressão prognóstica, não novo PAM50", ha="center",
            fontsize=8.7, style="italic", color=C["stone"])

    rounded(ax, (10.72, 5.15), (2.56, 0.78), fc=C["white"], ec=C["navy"], lw=1.1, radius=0.035)
    ax.text(12.0, 5.54, "PAM50 canônico", ha="center",
            fontsize=9.1, weight="bold", color=C["navy"])
    ax.text(12.0, 5.27, "50 genes + centroides", ha="center",
            fontsize=8.2, color=C["stone"])

    rounded(ax, (11.06, 4.08), (1.88, 0.58), fc=C["cream"], ec=C["sand"], lw=0.8, radius=0.03)
    ax.text(12.0, 4.37, "seleção penalizada", ha="center",
            fontsize=8.2, color=C["ink"])
    arrow(ax, 12.0, 5.12, 12.0, 4.69, color=C["stone"])
    arrow(ax, 12.0, 4.04, 12.0, 3.62, color=C["stone"])

    rounded(ax, (10.72, 2.78), (2.56, 0.78), fc="#F9EFEA", ec=C["accent"], lw=1.3, radius=0.035)
    ax.text(12.0, 3.17, "CorePAM", ha="center",
            fontsize=9.4, weight="bold", color=C["accent"])
    ax.text(12.0, 2.90, "24 genes prognósticos", ha="center",
            fontsize=8.2, color=C["stone"])

    rounded(ax, (10.58, 1.35), (2.84, 0.90), fc=C["white"], ec=C["sand"], lw=0.8, radius=0.03)
    ax.text(12.0, 1.93, "mantém o sinal de risco\nem coortes e plataformas externas",
            ha="center", va="center", fontsize=8.3, color=C["ink"], linespacing=1.18)

    arrow(ax, 4.05, 4.05, 5.00, 4.05, color=C["navy"], lw=1.8)
    arrow(ax, 9.25, 4.05, 10.18, 4.05, color=C["accent"], lw=1.8)

    ax.text(7.1, 0.42,
            "Valores de correlação são ilustrativos. IHQ: imuno-histoquímica. "
            "O CorePAM não redefine os centroides originais do PAM50; ele reduz o painel para prognóstico.",
            ha="center", fontsize=8.2, style="italic", color=C["stone"])

    name = "B4_PAM50_Centroides_CorePAM_PT"
    fig.savefig(OUT_PNG / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_PDF / f"{name}.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"OK {name}")


if __name__ == "__main__":
    main()
