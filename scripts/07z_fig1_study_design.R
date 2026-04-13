# =============================================================================
# SCRIPT: 07z_fig1_study_design.R
# PURPOSE: Generate Fig1 — Study design + cohort overview + frozen parameters
#          Shows: PAM50 → CorePAM derivation pipeline, cohort roles, N/events,
#          frozen parameters (ΔC=0.010, α=0.5, K=10, SHA-256 folds)
#
# OUTPUTS:
#   figures/main/en/pdf/Fig1_StudyDesign_EN.pdf
#   figures/main/en/png/Fig1_StudyDesign_EN.png
#   figures/main/pt/pdf/Fig1_StudyDesign_PT.pdf
#   figures/main/pt/png/Fig1_StudyDesign_PT.png
# PROJECT: Core-PAM (Memorial v6.1 / Freeze Core-PAM)
# =============================================================================

source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(patchwork)
})

SCRIPT_NAME <- "07z_fig1_study_design.R"

FORCE <- as.logical(Sys.getenv("FORCE_RERUN", "FALSE"))
out_en <- "figures/main/Fig1_StudyDesign_EN.pdf"
if (!FORCE && file.exists(out_en)) {
  message(sprintf("[%s] Output exists — skipping. Set FORCE_RERUN=TRUE.", SCRIPT_NAME))
  quit(save = "no", status = 0)
}
dir.create("figures/main", showWarnings = FALSE, recursive = TRUE)

# --------------------------------------------------------------------------
# Colors
# --------------------------------------------------------------------------
COL_TRAIN  <- "#2166AC"   # blue — training
COL_VAL    <- "#D6604D"   # red — validation
COL_PCR    <- "#4DAC26"   # green — pCR secondary
COL_CORE   <- "#762A83"   # purple — CorePAM
COL_FREEZE <- "#F4A582"   # orange — frozen params box
COL_ARROW  <- "#525252"   # dark grey arrows
COL_BG     <- "white"

# --------------------------------------------------------------------------
# Build plot via ggplot2 geom_rect + geom_text + geom_segment (arrows)
# Coordinate system: x 0-100, y 0-100
# --------------------------------------------------------------------------

# Helper: box
box <- function(xmin, xmax, ymin, ymax, fill, color = "white",
                alpha = 1, linetype = "solid") {
  annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
           fill = fill, color = color, alpha = alpha, linewidth = 0.6,
           linetype = linetype)
}
# Helper: label
lbl <- function(x, y, label, size = 3.2, fontface = "plain", color = "white",
                hjust = 0.5, vjust = 0.5, ...) {
  annotate("text", x = x, y = y, label = label, size = size,
           fontface = fontface, color = color, hjust = hjust, vjust = vjust, ...)
}
# Helper: arrow
arr <- function(x1, y1, x2, y2) {
  annotate("segment", x = x1, xend = x2, y = y1, yend = y2,
           arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
           color = COL_ARROW, linewidth = 0.55)
}

make_fig1 <- function(lang = "EN") {
  tr   <- if (lang == "EN") list(
    title       = "Study Design Overview",
    training    = "Training",
    validation  = "External Validation",
    pcr_block   = "Secondary Analysis\n(pCR — chemotherapy response)",
    scanb       = "SCAN-B\n(RNA-seq)\nN=3,069 | Events=322",
    tcga        = "TCGA-BRCA\n(RNA-seq)\nN=1,072 | Events=150",
    metabric    = "METABRIC\n(Microarray)\nN=1,978 | DSS Events=646",
    gse20685    = "GSE20685\n(Microarray)\nN=327 | Events=83",
    pam50       = "PAM50\n(50 genes)",
    corepam     = "CorePAM\n(24 genes)",
    method      = "Cox Elastic-Net\nOOF C-index Non-Inferiority\nΔC ≤ 0.010",
    freeze      = "Frozen Parameters\nα=0.5 | K=10 | seed=SHA-256(ID)\nδC=0.010 | t in months | Z-score intra-cohort",
    endpoint_os = "Primary: OS/DSS",
    endpoint_pcr= "Endpoint: pCR",
    gse25       = "GSE25066 | N=182",
    gse20194    = "GSE20194 | N=278",
    gse32646    = "GSE32646 | N=115",
    ispy1       = "I-SPY1 | N=122",
    gse1456     = "GSE1456\n(Microarray)\nN=159 | Events=40",
    ispy2       = "I-SPY2*\n(Microarray)\nN=986 | Exploratory",
    score       = "Score = Σ(wᵢ·zᵢ) / Σ|wᵢ|"
  ) else list(
    title       = "Desenho do Estudo",
    training    = "Derivação",
    validation  = "Validação Externa",
    pcr_block   = "Análise Secundária\n(pCR — resposta à quimioterapia)",
    scanb       = "SCAN-B\n(RNA-seq)\nN=3.069 | Eventos=322",
    tcga        = "TCGA-BRCA\n(RNA-seq)\nN=1.072 | Eventos=150",
    metabric    = "METABRIC\n(Microarray)\nN=1.978 | Eventos DSS=646",
    gse20685    = "GSE20685\n(Microarray)\nN=327 | Eventos=83",
    pam50       = "PAM50\n(50 genes)",
    corepam     = "Core-PAM\n(24 genes)",
    method      = "Cox Elastic-Net\nC-index OOF Não-Inferioridade\nΔC ≤ 0,010",
    freeze      = "Parâmetros Congelados\nα=0,5 | K=10 | semente=SHA-256(ID)\nδC=0,010 | t em meses | Z-score intra-coorte",
    endpoint_os = "Primário: OS/DSS",
    endpoint_pcr= "Desfecho: pCR",
    gse25       = "GSE25066 | N=182",
    gse20194    = "GSE20194 | N=278",
    gse32646    = "GSE32646 | N=115",
    gse1456     = "GSE1456\n(Microarray)\nN=159 | Eventos=40",
    ispy1       = "I-SPY1 | N=122",
    ispy2       = "I-SPY2*\n(Microarray)\nN=986 | Exploratória",
    score       = "Score = Σ(wᵢ·zᵢ) / Σ|wᵢ|"
  )

  p <- ggplot() +
    theme_void() +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100), expand = FALSE) +
    theme(plot.background  = element_rect(fill = COL_BG, color = NA),
          panel.background = element_rect(fill = COL_BG, color = NA),
          plot.margin      = margin(4, 4, 4, 4, "mm"))

  # ---- TITLE row (y 94-100) ----
  p <- p + lbl(50, 97, tr$title, size = 5, fontface = "bold", color = "#252525")

  # ---- LEFT: PAM50 → CorePAM derivation (x 2-38, y 52-90) ----
  # Derivation section header
  p <- p +
    box(2, 38, 80, 90, fill = COL_TRAIN, alpha = 0.92) +
    lbl(20, 85.5, tr$training, size = 3.6, fontface = "bold") +
    lbl(20, 82, tr$scanb, size = 2.7)

  # PAM50 box
  p <- p +
    box(2, 17, 62, 78, fill = COL_CORE, alpha = 0.20, color = COL_CORE) +
    lbl(9.5, 74, tr$pam50, size = 3.3, fontface = "bold", color = COL_CORE) +
    lbl(9.5, 67, tr$method, size = 2.4, color = "#333333")

  # Arrow PAM50 → CorePAM
  p <- p + arr(17, 70, 23, 70)

  # CorePAM box
  p <- p +
    box(23, 37, 62, 78, fill = COL_CORE, alpha = 0.85) +
    lbl(30, 74, tr$corepam, size = 3.5, fontface = "bold") +
    lbl(30, 68.5, tr$score, size = 2.1)

  # Freeze params (bottom left)
  p <- p +
    box(2, 38, 51, 60, fill = COL_FREEZE, alpha = 0.9, color = "#D8A060") +
    lbl(20, 55.5, tr$freeze, size = 2.3, color = "#4D2600")

  # ---- CENTER: arrow from CorePAM to validation (x 38-52) ----
  p <- p + arr(37, 70, 43, 70)
  p <- p + lbl(40, 72.5, tr$endpoint_os, size = 2.2, color = "#333333", fontface = "italic")

  # ---- RIGHT: Validation cohorts (x 43-98, y 52-90) ----
  # Validation section header
  p <- p +
    box(43, 98, 80, 90, fill = COL_VAL, alpha = 0.88) +
    lbl(70.5, 85.5, tr$validation, size = 3.6, fontface = "bold") +
    lbl(70.5, 82, tr$endpoint_os, size = 2.7)

  # Validation cohort boxes (4 cohorts)
  for (i in seq_along(c("tcga", "metabric", "gse20685", "gse1456"))) {
    cohort_lbl <- c(tr$tcga, tr$metabric, tr$gse20685, tr$gse1456)[i]
    x1 <- c(43, 58, 73, 86)[i]; x2 <- c(57, 72, 85, 98)[i]
    p <- p +
      box(x1, x2, 62, 77.5, fill = COL_VAL, alpha = 0.65) +
      lbl((x1 + x2) / 2, 69.5, cohort_lbl, size = 2.15)
  }
  # Rebuild cleanly
  p2 <- ggplot() +
    theme_void() +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100), expand = FALSE) +
    theme(plot.background  = element_rect(fill = COL_BG, color = NA),
          panel.background = element_rect(fill = COL_BG, color = NA),
          plot.margin      = margin(3, 3, 3, 3, "mm")) +

    # Title
    annotate("text", x = 50, y = 97.5, label = tr$title,
             size = 5.5, fontface = "bold", color = "#1a1a1a", hjust = 0.5) +

    # ---- TRAINING block (left, y 48-93) ----
    annotate("rect", xmin = 2, xmax = 37, ymin = 48, ymax = 93,
             fill = COL_TRAIN, color = "white", alpha = 0.1, linewidth = 1) +
    annotate("text", x = 19.5, y = 91, label = tr$training,
             size = 3.8, fontface = "bold", color = COL_TRAIN, hjust = 0.5) +

    # SCAN-B box
    annotate("rect", xmin = 3, xmax = 36, ymin = 74, ymax = 88,
             fill = COL_TRAIN, color = NA, alpha = 0.85) +
    annotate("text", x = 19.5, y = 82, label = tr$scanb,
             size = 3, color = "white", hjust = 0.5, lineheight = 1.3) +

    # PAM50 → method → CorePAM
    annotate("rect", xmin = 3, xmax = 18, ymin = 55, ymax = 72,
             fill = COL_CORE, color = NA, alpha = 0.18) +
    annotate("text", x = 10.5, y = 69, label = tr$pam50,
             size = 3.2, fontface = "bold", color = COL_CORE, hjust = 0.5) +
    annotate("text", x = 10.5, y = 62, label = tr$method,
             size = 2.25, color = "#444444", hjust = 0.5, lineheight = 1.25) +

    # Arrow PAM50→CorePAM
    annotate("segment", x = 18, xend = 22, y = 63.5, yend = 63.5,
             arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
             color = COL_CORE, linewidth = 0.8) +

    annotate("rect", xmin = 22, xmax = 36, ymin = 55, ymax = 72,
             fill = COL_CORE, color = NA, alpha = 0.88) +
    annotate("text", x = 29, y = 68, label = tr$corepam,
             size = 3.5, fontface = "bold", color = "white", hjust = 0.5) +
    annotate("text", x = 29, y = 61.5, label = tr$score,
             size = 2.1, color = "white", hjust = 0.5) +

    # Frozen params box
    annotate("rect", xmin = 3, xmax = 36, ymin = 48.5, ymax = 54,
             fill = COL_FREEZE, color = "#B07830", alpha = 0.9, linewidth = 0.5) +
    annotate("text", x = 19.5, y = 51.2, label = tr$freeze,
             size = 2.15, color = "#3D1A00", hjust = 0.5, lineheight = 1.3) +

    # Arrow from CorePAM to validation (through center)
    annotate("segment", x = 36, xend = 40, y = 63.5, yend = 63.5,
             arrow = arrow(length = unit(0.22, "cm"), type = "closed"),
             color = COL_ARROW, linewidth = 0.7) +

    # ---- VALIDATION block (right, y 48-93) ----
    annotate("rect", xmin = 40, xmax = 98, ymin = 48, ymax = 93,
             fill = COL_VAL, color = "white", alpha = 0.08, linewidth = 1) +
    annotate("text", x = 69, y = 91, label = tr$validation,
             size = 3.8, fontface = "bold", color = COL_VAL, hjust = 0.5) +

    # Validation cohort boxes (4 cohorts side by side)
    annotate("rect", xmin = 41, xmax = 55, ymin = 55, ymax = 88,
             fill = COL_VAL, color = NA, alpha = 0.82) +
    annotate("text", x = 48, y = 71.5, label = tr$tcga,
             size = 2.5, color = "white", hjust = 0.5, lineheight = 1.3) +

    annotate("rect", xmin = 56, xmax = 72, ymin = 55, ymax = 88,
             fill = COL_VAL, color = NA, alpha = 0.82) +
    annotate("text", x = 64, y = 71.5, label = tr$metabric,
             size = 2.5, color = "white", hjust = 0.5, lineheight = 1.3) +

    annotate("rect", xmin = 73, xmax = 85, ymin = 55, ymax = 88,
             fill = COL_VAL, color = NA, alpha = 0.82) +
    annotate("text", x = 79, y = 71.5, label = tr$gse20685,
             size = 2.5, color = "white", hjust = 0.5, lineheight = 1.3) +

    annotate("rect", xmin = 86, xmax = 97, ymin = 55, ymax = 88,
             fill = COL_VAL, color = NA, alpha = 0.82) +
    annotate("text", x = 91.5, y = 71.5, label = tr$gse1456,
             size = 2.5, color = "white", hjust = 0.5, lineheight = 1.3) +

    annotate("text", x = 69, y = 50.5, label = tr$endpoint_os,
             size = 2.5, color = COL_VAL, hjust = 0.5, fontface = "italic") +

    # ---- pCR block (bottom, y 4-44) ----
    annotate("rect", xmin = 2, xmax = 98, ymin = 4, ymax = 44,
             fill = COL_PCR, color = "white", alpha = 0.08, linewidth = 0.8,
             linetype = "dashed") +
    annotate("text", x = 50, y = 42, label = tr$pcr_block,
             size = 3.4, fontface = "bold.italic", color = COL_PCR, hjust = 0.5) +

    # 5 pCR cohort boxes (4 primary + 1 exploratory)
    annotate("rect", xmin = 3,  xmax = 21, ymin = 10, ymax = 38,
             fill = COL_PCR, color = NA, alpha = 0.8) +
    annotate("text", x = 12, y = 24, label = tr$gse25,
             size = 2.3, color = "white", hjust = 0.5) +

    annotate("rect", xmin = 22, xmax = 40, ymin = 10, ymax = 38,
             fill = COL_PCR, color = NA, alpha = 0.8) +
    annotate("text", x = 31, y = 24, label = tr$gse20194,
             size = 2.3, color = "white", hjust = 0.5) +

    annotate("rect", xmin = 41, xmax = 57, ymin = 10, ymax = 38,
             fill = COL_PCR, color = NA, alpha = 0.8) +
    annotate("text", x = 49, y = 24, label = tr$gse32646,
             size = 2.3, color = "white", hjust = 0.5) +

    annotate("rect", xmin = 58, xmax = 74, ymin = 10, ymax = 38,
             fill = COL_PCR, color = NA, alpha = 0.8) +
    annotate("text", x = 66, y = 24, label = tr$ispy1,
             size = 2.3, color = "white", hjust = 0.5) +

    annotate("rect", xmin = 75, xmax = 97, ymin = 10, ymax = 38,
             fill = COL_PCR, color = NA, alpha = 0.5, linetype = "dashed") +
    annotate("text", x = 86, y = 24, label = tr$ispy2,
             size = 2.1, color = "white", hjust = 0.5, lineheight = 1.3) +

    annotate("text", x = 50, y = 7, label = tr$endpoint_pcr,
             size = 2.5, color = COL_PCR, hjust = 0.5, fontface = "italic") +

    # Arrow from CorePAM down to pCR block
    annotate("segment", x = 29, xend = 29, y = 55, yend = 44,
             arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
             color = COL_PCR, linewidth = 0.7, linetype = "dashed")

  p2
}

# --------------------------------------------------------------------------
# Save EN and PT
# --------------------------------------------------------------------------
message(sprintf("[%s] Generating Fig1 EN...", SCRIPT_NAME))
p_en <- make_fig1("EN")
ggsave("figures/main/Fig1_StudyDesign_EN.pdf", p_en,
       width = 170, height = 140, units = "mm", device = cairo_pdf)
ggsave("figures/main/Fig1_StudyDesign_EN.png", p_en,
       width = 170, height = 140, units = "mm", dpi = 300, bg = "white")
message(sprintf("[%s] EN saved", SCRIPT_NAME))

message(sprintf("[%s] Generating Fig1 PT...", SCRIPT_NAME))
p_pt <- make_fig1("PT")
ggsave("figures/main/Fig1_StudyDesign_PT.pdf", p_pt,
       width = 170, height = 140, units = "mm", device = cairo_pdf)
ggsave("figures/main/Fig1_StudyDesign_PT.png", p_pt,
       width = 170, height = 140, units = "mm", dpi = 300, bg = "white")
message(sprintf("[%s] PT saved", SCRIPT_NAME))

# Registry
for (f in c("figures/main/Fig1_StudyDesign_EN.pdf",
            "figures/main/Fig1_StudyDesign_PT.pdf")) {
  h  <- sha256_file(f)
  sz <- file.info(f)$size / 1e6
  registry_append("STUDY_DESIGN", "fig1_study_design", f, h, "ok", SCRIPT_NAME, sz, list())
}

message(sprintf("[%s] DONE — Fig1_StudyDesign EN+PT generated", SCRIPT_NAME))
