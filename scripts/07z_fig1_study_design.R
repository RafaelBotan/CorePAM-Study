# =============================================================================
# SCRIPT: 07z_fig1_study_design.R
# PURPOSE: Generate Fig1 — Study design + cohort overview + frozen parameters
#          Full-page vertical flowchart: PAM50 → CorePAM derivation, training,
#          validation cohorts, pCR secondary analysis, frozen parameters.
#
# OUTPUTS:
#   figures/main/en/pdf/Fig1_StudyDesign_EN.pdf
#   figures/main/en/png/Fig1_StudyDesign_EN.png
#   figures/main/pt/pdf/Fig1_StudyDesign_PT.pdf
#   figures/main/pt/png/Fig1_StudyDesign_PT.png
# PROJETO: Core-PAM (Memorial v6.1 / Freeze Core-PAM)
# =============================================================================

source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

SCRIPT_NAME <- "07z_fig1_study_design.R"

FORCE <- as.logical(Sys.getenv("FORCE_RERUN", "FALSE"))
out_en <- "figures/main/en/pdf/Fig1_StudyDesign_EN.pdf"
if (!FORCE && file.exists(out_en)) {
  message(sprintf("[%s] Output exists — skipping. Set FORCE_RERUN=TRUE.", SCRIPT_NAME))
  quit(save = "no", status = 0)
}

# --------------------------------------------------------------------------
# Colors
# --------------------------------------------------------------------------
COL_TRAIN  <- "#2166AC"   # blue — training
COL_VAL    <- "#D6604D"   # red — validation
COL_PCR    <- "#4DAC26"   # green — pCR secondary
COL_CORE   <- "#762A83"   # purple — CorePAM
COL_FREEZE <- "#FFF3E0"   # light orange — frozen params box
COL_ARROW  <- "#424242"   # dark grey arrows
COL_BG     <- "white"
COL_BORDER <- "#BDBDBD"   # light grey section borders

# --------------------------------------------------------------------------
# Build full-page vertical flowchart
# Coordinate system: x 0-100, y 0-100 (saved at 200×280mm portrait)
# --------------------------------------------------------------------------

make_fig1 <- function(lang = "EN") {

  # ---- Translations ----
  tr <- if (lang == "EN") list(
    title       = "Study Design Overview",
    subtitle    = "CorePAM: derivation, validation, and secondary pCR analysis",
    # Derivation row
    pam50_title = "PAM50",
    pam50_sub   = "50 genes",
    method_l1   = "Cox Elastic-Net",
    method_l2   = "\u03b1 = 0.5  |  K = 10 folds",
    method_l3   = "Folds: SHA-256(patient ID)",
    method_l4   = "Non-inferiority: \u0394C \u2264 0.010",
    corepam_title = "CorePAM",
    corepam_sub = "24 genes",
    score       = "Score = \u03a3(w\u1d62\u00b7z\u1d62) / \u03a3|w\u1d62|",
    score_note  = "Z-score intra-cohort | genes present \u2229 cohort",
    # Training
    train_title = "TRAINING (Derivation)",
    scanb_name  = "SCAN-B",
    scanb_plat  = "RNA-seq (GEO: GSE96058)",
    scanb_n     = "N = 3,069  |  Events = 322 (10.5%)",
    scanb_fu    = "Median follow-up: 54.9 months",
    scanb_ep    = "Endpoint: Overall Survival (OS)",
    # Validation
    val_title   = "EXTERNAL VALIDATION",
    val_ep      = "Primary endpoints: OS / DSS",
    tcga_name   = "TCGA-BRCA",
    tcga_plat   = "RNA-seq",
    tcga_n      = "N = 1,072",
    tcga_ev     = "Events = 150 (14.0%)",
    tcga_fu     = "FU: 25.6 mo",
    tcga_ep     = "OS",
    metab_name  = "METABRIC",
    metab_plat  = "Microarray (Illumina)",
    metab_n     = "N = 1,980",
    metab_ev    = "Events = 1,144 (57.8%)",
    metab_fu    = "FU: 157.9 mo",
    metab_ep    = "DSS",
    gse_name    = "GSE20685",
    gse_plat    = "Microarray (Affymetrix)",
    gse_n       = "N = 327",
    gse_ev      = "Events = 83 (25.4%)",
    gse_fu      = "FU: 110.4 mo",
    gse_ep      = "OS",
    # pCR
    pcr_title   = "SECONDARY ANALYSIS \u2014 pCR (Chemotherapy Response)",
    pcr_ep      = "Endpoint: pathological complete response (pCR)",
    gse25_name  = "GSE25066",
    gse25_plat  = "HGU133A",
    gse25_n     = "N = 508",
    gse20194_name = "GSE20194",
    gse20194_plat = "HGU133Plus2",
    gse20194_n  = "N = 278",
    gse32_name  = "GSE32646",
    gse32_plat  = "HGU133Plus2",
    gse32_n     = "N = 154",
    ispy_name   = "I-SPY1",
    ispy_plat   = "Agilent 44K",
    ispy_n      = "N = 149",
    # Freeze
    freeze_title = "FROZEN ANALYSIS PARAMETERS",
    freeze_l1   = "\u03b1 = 0.5 (elastic-net mixing)   |   K = 10 (CV folds)   |   seed = SHA-256(patient ID)",
    freeze_l2   = "\u0394C = 0.010 (non-inferiority margin)   |   time in months (days / 30.4375)   |   Z-score intra-cohort"
  ) else list(
    title       = "Desenho do Estudo",
    subtitle    = "CorePAM: deriva\u00e7\u00e3o, valida\u00e7\u00e3o e an\u00e1lise secund\u00e1ria de pCR",
    pam50_title = "PAM50",
    pam50_sub   = "50 genes",
    method_l1   = "Cox Elastic-Net",
    method_l2   = "\u03b1 = 0,5  |  K = 10 folds",
    method_l3   = "Folds: SHA-256(ID do paciente)",
    method_l4   = "N\u00e3o-inferioridade: \u0394C \u2264 0,010",
    corepam_title = "CorePAM",
    corepam_sub = "24 genes",
    score       = "Score = \u03a3(w\u1d62\u00b7z\u1d62) / \u03a3|w\u1d62|",
    score_note  = "Z-score intra-coorte | genes presentes \u2229 coorte",
    train_title = "DERIVA\u00c7\u00c3O (Treinamento)",
    scanb_name  = "SCAN-B",
    scanb_plat  = "RNA-seq (GEO: GSE96058)",
    scanb_n     = "N = 3.069  |  Eventos = 322 (10,5%)",
    scanb_fu    = "Mediana follow-up: 54,9 meses",
    scanb_ep    = "Desfecho: Sobrevida Global (OS)",
    val_title   = "VALIDA\u00c7\u00c3O EXTERNA",
    val_ep      = "Desfechos prim\u00e1rios: OS / SLE",
    tcga_name   = "TCGA-BRCA",
    tcga_plat   = "RNA-seq",
    tcga_n      = "N = 1.072",
    tcga_ev     = "Eventos = 150 (14,0%)",
    tcga_fu     = "FU: 25,6 m",
    tcga_ep     = "OS",
    metab_name  = "METABRIC",
    metab_plat  = "Microarray (Illumina)",
    metab_n     = "N = 1.980",
    metab_ev    = "Eventos = 1.144 (57,8%)",
    metab_fu    = "FU: 157,9 m",
    metab_ep    = "SLE",
    gse_name    = "GSE20685",
    gse_plat    = "Microarray (Affymetrix)",
    gse_n       = "N = 327",
    gse_ev      = "Eventos = 83 (25,4%)",
    gse_fu      = "FU: 110,4 m",
    gse_ep      = "OS",
    pcr_title   = "AN\u00c1LISE SECUND\u00c1RIA \u2014 pCR (Resposta \u00e0 Quimioterapia)",
    pcr_ep      = "Desfecho: resposta patol\u00f3gica completa (pCR)",
    gse25_name  = "GSE25066",
    gse25_plat  = "HGU133A",
    gse25_n     = "N = 508",
    gse20194_name = "GSE20194",
    gse20194_plat = "HGU133Plus2",
    gse20194_n  = "N = 278",
    gse32_name  = "GSE32646",
    gse32_plat  = "HGU133Plus2",
    gse32_n     = "N = 154",
    ispy_name   = "I-SPY1",
    ispy_plat   = "Agilent 44K",
    ispy_n      = "N = 149",
    freeze_title = "PAR\u00c2METROS CONGELADOS DA AN\u00c1LISE",
    freeze_l1   = "\u03b1 = 0,5 (elastic-net)   |   K = 10 (folds CV)   |   semente = SHA-256(ID do paciente)",
    freeze_l2   = "\u0394C = 0,010 (margem de n\u00e3o-inferioridade)   |   tempo em meses (dias / 30,4375)   |   Z-score intra-coorte"
  )

  # ---- Build plot ----
  p <- ggplot() +
    theme_void() +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100), expand = FALSE) +
    theme(plot.background  = element_rect(fill = COL_BG, color = NA),
          panel.background = element_rect(fill = COL_BG, color = NA),
          plot.margin      = margin(5, 5, 5, 5, "mm"))

  # ============================================================
  # TITLE (y 95-99)
  # ============================================================
  p <- p +
    annotate("text", x = 50, y = 98, label = tr$title,
             size = 7, fontface = "bold", color = "#1a1a1a") +
    annotate("text", x = 50, y = 95.5, label = tr$subtitle,
             size = 4, color = "#555555", fontface = "italic")

  # ============================================================
  # ROW 1: DERIVATION FLOW (y 82-93) — PAM50 → Method → CorePAM
  # ============================================================
  # PAM50 box (x 3-22)
  p <- p +
    annotate("rect", xmin = 3, xmax = 22, ymin = 83, ymax = 93,
             fill = COL_CORE, alpha = 0.15, color = COL_CORE, linewidth = 0.8) +
    annotate("text", x = 12.5, y = 90.5, label = tr$pam50_title,
             size = 5, fontface = "bold", color = COL_CORE) +
    annotate("text", x = 12.5, y = 87.5, label = tr$pam50_sub,
             size = 3.5, color = COL_CORE) +
    annotate("text", x = 12.5, y = 84.8, label = tr$method_l1,
             size = 2.8, color = "#555555")

  # Arrow PAM50 → Method box
  p <- p +
    annotate("segment", x = 22, xend = 28, y = 88, yend = 88,
             arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
             color = COL_CORE, linewidth = 1)

  # Method box (x 28-62)
  p <- p +
    annotate("rect", xmin = 28, xmax = 62, ymin = 83, ymax = 93,
             fill = "#F3E5F5", color = COL_CORE, linewidth = 0.6,
             linetype = "dashed") +
    annotate("text", x = 45, y = 91, label = tr$method_l2,
             size = 3, color = "#333333") +
    annotate("text", x = 45, y = 88.5, label = tr$method_l3,
             size = 3, color = "#333333") +
    annotate("text", x = 45, y = 86, label = tr$method_l4,
             size = 3, fontface = "bold", color = COL_CORE)

  # Arrow Method → CorePAM
  p <- p +
    annotate("segment", x = 62, xend = 68, y = 88, yend = 88,
             arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
             color = COL_CORE, linewidth = 1)

  # CorePAM box (x 68-97)
  p <- p +
    annotate("rect", xmin = 68, xmax = 97, ymin = 83, ymax = 93,
             fill = COL_CORE, color = NA, alpha = 0.9) +
    annotate("text", x = 82.5, y = 90.5, label = tr$corepam_title,
             size = 5.5, fontface = "bold", color = "white") +
    annotate("text", x = 82.5, y = 88, label = tr$corepam_sub,
             size = 3.5, color = "white") +
    annotate("text", x = 82.5, y = 85.5, label = tr$score,
             size = 2.8, color = "#E1BEE7") +
    annotate("text", x = 82.5, y = 83.8, label = tr$score_note,
             size = 2.2, color = "#CE93D8")

  # ============================================================
  # ARROW from derivation to training
  # ============================================================
  p <- p +
    annotate("segment", x = 50, xend = 50, y = 83, yend = 80,
             arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
             color = COL_TRAIN, linewidth = 1)

  # ============================================================
  # ROW 2: TRAINING (y 67-79)
  # ============================================================
  p <- p +
    # Section border
    annotate("rect", xmin = 3, xmax = 97, ymin = 67, ymax = 79,
             fill = COL_TRAIN, alpha = 0.06, color = COL_TRAIN,
             linewidth = 1) +
    # Header
    annotate("text", x = 50, y = 77.5, label = tr$train_title,
             size = 4.5, fontface = "bold", color = COL_TRAIN) +
    # SCAN-B box
    annotate("rect", xmin = 10, xmax = 90, ymin = 68.5, ymax = 76,
             fill = COL_TRAIN, color = NA, alpha = 0.85) +
    annotate("text", x = 50, y = 74.2, label = paste(tr$scanb_name, " \u2014 ", tr$scanb_plat),
             size = 3.5, fontface = "bold", color = "white") +
    annotate("text", x = 50, y = 72, label = tr$scanb_n,
             size = 3.2, color = "white") +
    annotate("text", x = 50, y = 70, label = paste(tr$scanb_fu, " | ", tr$scanb_ep),
             size = 3, color = "#BBDEFB")

  # ============================================================
  # ARROW from training to validation
  # ============================================================
  p <- p +
    annotate("segment", x = 50, xend = 50, y = 67, yend = 64,
             arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
             color = COL_VAL, linewidth = 1)

  # ============================================================
  # ROW 3: EXTERNAL VALIDATION (y 42-63)
  # ============================================================
  p <- p +
    annotate("rect", xmin = 3, xmax = 97, ymin = 42, ymax = 63,
             fill = COL_VAL, alpha = 0.05, color = COL_VAL,
             linewidth = 1) +
    annotate("text", x = 50, y = 61.5, label = tr$val_title,
             size = 4.5, fontface = "bold", color = COL_VAL) +
    annotate("text", x = 50, y = 59.2, label = tr$val_ep,
             size = 3, fontface = "italic", color = COL_VAL)

  # TCGA-BRCA (x 5-34)
  p <- p +
    annotate("rect", xmin = 5, xmax = 34, ymin = 43.5, ymax = 57.5,
             fill = COL_VAL, color = NA, alpha = 0.82) +
    annotate("text", x = 19.5, y = 56, label = tr$tcga_name,
             size = 3.8, fontface = "bold", color = "white") +
    annotate("text", x = 19.5, y = 53.5, label = tr$tcga_plat,
             size = 2.8, color = "#FFCDD2") +
    annotate("text", x = 19.5, y = 51.5, label = tr$tcga_n,
             size = 3, color = "white") +
    annotate("text", x = 19.5, y = 49.5, label = tr$tcga_ev,
             size = 2.8, color = "white") +
    annotate("text", x = 19.5, y = 47.5, label = paste(tr$tcga_fu, " | ", tr$tcga_ep),
             size = 2.6, color = "#FFCDD2")

  # METABRIC (x 36-65)
  p <- p +
    annotate("rect", xmin = 36, xmax = 65, ymin = 43.5, ymax = 57.5,
             fill = COL_VAL, color = NA, alpha = 0.82) +
    annotate("text", x = 50.5, y = 56, label = tr$metab_name,
             size = 3.8, fontface = "bold", color = "white") +
    annotate("text", x = 50.5, y = 53.5, label = tr$metab_plat,
             size = 2.8, color = "#FFCDD2") +
    annotate("text", x = 50.5, y = 51.5, label = tr$metab_n,
             size = 3, color = "white") +
    annotate("text", x = 50.5, y = 49.5, label = tr$metab_ev,
             size = 2.8, color = "white") +
    annotate("text", x = 50.5, y = 47.5, label = paste(tr$metab_fu, " | ", tr$metab_ep),
             size = 2.6, color = "#FFCDD2")

  # GSE20685 (x 67-96)
  p <- p +
    annotate("rect", xmin = 67, xmax = 96, ymin = 43.5, ymax = 57.5,
             fill = COL_VAL, color = NA, alpha = 0.82) +
    annotate("text", x = 81.5, y = 56, label = tr$gse_name,
             size = 3.8, fontface = "bold", color = "white") +
    annotate("text", x = 81.5, y = 53.5, label = tr$gse_plat,
             size = 2.8, color = "#FFCDD2") +
    annotate("text", x = 81.5, y = 51.5, label = tr$gse_n,
             size = 3, color = "white") +
    annotate("text", x = 81.5, y = 49.5, label = tr$gse_ev,
             size = 2.8, color = "white") +
    annotate("text", x = 81.5, y = 47.5, label = paste(tr$gse_fu, " | ", tr$gse_ep),
             size = 2.6, color = "#FFCDD2")

  # ============================================================
  # ARROW from validation to pCR
  # ============================================================
  p <- p +
    annotate("segment", x = 50, xend = 50, y = 42, yend = 39,
             arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
             color = COL_PCR, linewidth = 1, linetype = "dashed")

  # ============================================================
  # ROW 4: pCR BLOCK (y 18-38)
  # ============================================================
  p <- p +
    annotate("rect", xmin = 3, xmax = 97, ymin = 18, ymax = 38,
             fill = COL_PCR, alpha = 0.06, color = COL_PCR,
             linewidth = 0.8, linetype = "dashed") +
    annotate("text", x = 50, y = 36.5, label = tr$pcr_title,
             size = 3.8, fontface = "bold.italic", color = COL_PCR) +
    annotate("text", x = 50, y = 34.5, label = tr$pcr_ep,
             size = 2.8, fontface = "italic", color = COL_PCR)

  # 4 pCR cohort boxes
  pcr_y1 <- 19.5; pcr_y2 <- 33
  # GSE25066
  p <- p +
    annotate("rect", xmin = 4, xmax = 27, ymin = pcr_y1, ymax = pcr_y2,
             fill = COL_PCR, color = NA, alpha = 0.8) +
    annotate("text", x = 15.5, y = 31, label = tr$gse25_name,
             size = 3.3, fontface = "bold", color = "white") +
    annotate("text", x = 15.5, y = 28.5, label = tr$gse25_plat,
             size = 2.5, color = "#C5E1A5") +
    annotate("text", x = 15.5, y = 26, label = tr$gse25_n,
             size = 2.8, color = "white")

  # GSE20194
  p <- p +
    annotate("rect", xmin = 28, xmax = 50, ymin = pcr_y1, ymax = pcr_y2,
             fill = COL_PCR, color = NA, alpha = 0.8) +
    annotate("text", x = 39, y = 31, label = tr$gse20194_name,
             size = 3.3, fontface = "bold", color = "white") +
    annotate("text", x = 39, y = 28.5, label = tr$gse20194_plat,
             size = 2.5, color = "#C5E1A5") +
    annotate("text", x = 39, y = 26, label = tr$gse20194_n,
             size = 2.8, color = "white")

  # GSE32646
  p <- p +
    annotate("rect", xmin = 51, xmax = 73, ymin = pcr_y1, ymax = pcr_y2,
             fill = COL_PCR, color = NA, alpha = 0.8) +
    annotate("text", x = 62, y = 31, label = tr$gse32_name,
             size = 3.3, fontface = "bold", color = "white") +
    annotate("text", x = 62, y = 28.5, label = tr$gse32_plat,
             size = 2.5, color = "#C5E1A5") +
    annotate("text", x = 62, y = 26, label = tr$gse32_n,
             size = 2.8, color = "white")

  # I-SPY1
  p <- p +
    annotate("rect", xmin = 74, xmax = 97, ymin = pcr_y1, ymax = pcr_y2,
             fill = COL_PCR, color = NA, alpha = 0.8) +
    annotate("text", x = 85.5, y = 31, label = tr$ispy_name,
             size = 3.3, fontface = "bold", color = "white") +
    annotate("text", x = 85.5, y = 28.5, label = tr$ispy_plat,
             size = 2.5, color = "#C5E1A5") +
    annotate("text", x = 85.5, y = 26, label = tr$ispy_n,
             size = 2.8, color = "white")

  # ============================================================
  # ROW 5: FROZEN PARAMETERS (y 3-15)
  # ============================================================
  p <- p +
    annotate("rect", xmin = 3, xmax = 97, ymin = 3, ymax = 15,
             fill = COL_FREEZE, color = "#E0A050", linewidth = 0.8) +
    annotate("text", x = 50, y = 13, label = tr$freeze_title,
             size = 3.8, fontface = "bold", color = "#4E342E") +
    annotate("text", x = 50, y = 9.5, label = tr$freeze_l1,
             size = 3, color = "#5D4037") +
    annotate("text", x = 50, y = 6, label = tr$freeze_l2,
             size = 3, color = "#5D4037")

  p
}

# --------------------------------------------------------------------------
# Save EN and PT (full-page: 200 × 280mm)
# --------------------------------------------------------------------------
FIG_W <- 200   # mm — A4 text width
FIG_H <- 280   # mm — A4 text height (nearly full page)

dir.create("figures/main/en/pdf", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/main/en/png", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/main/pt/pdf", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/main/pt/png", showWarnings = FALSE, recursive = TRUE)

message(sprintf("[%s] Generating Fig1 EN...", SCRIPT_NAME))
p_en <- make_fig1("EN")
ggsave("figures/main/en/pdf/Fig1_StudyDesign_EN.pdf", p_en,
       width = FIG_W, height = FIG_H, units = "mm", device = cairo_pdf)
ggsave("figures/main/en/png/Fig1_StudyDesign_EN.png", p_en,
       width = FIG_W, height = FIG_H, units = "mm", dpi = 300, bg = "white")
message(sprintf("[%s] EN saved", SCRIPT_NAME))

message(sprintf("[%s] Generating Fig1 PT...", SCRIPT_NAME))
p_pt <- make_fig1("PT")
ggsave("figures/main/pt/pdf/Fig1_StudyDesign_PT.pdf", p_pt,
       width = FIG_W, height = FIG_H, units = "mm", device = cairo_pdf)
ggsave("figures/main/pt/png/Fig1_StudyDesign_PT.png", p_pt,
       width = FIG_W, height = FIG_H, units = "mm", dpi = 300, bg = "white")
message(sprintf("[%s] PT saved", SCRIPT_NAME))

# Registry
for (f in c("figures/main/en/pdf/Fig1_StudyDesign_EN.pdf",
            "figures/main/pt/pdf/Fig1_StudyDesign_PT.pdf")) {
  h  <- sha256_file(f)
  sz <- file.info(f)$size / 1e6
  registry_append("STUDY_DESIGN", "fig1_study_design", f, h, "ok", SCRIPT_NAME, sz, list())
}

message(sprintf("[%s] DONE — Fig1_StudyDesign EN+PT generated (full-page %dx%dmm)",
                SCRIPT_NAME, FIG_W, FIG_H))
