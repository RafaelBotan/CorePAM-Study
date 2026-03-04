# =============================================================================
# SCRIPT: 00_colors.R
# PURPOSE: Standardized color palette for all CorePAM project figures.
#          RULE: All figures must use exclusively these colors.
#          Palette: ColorBrewer RdBu (divergent) + Set1 qualitative (cohorts).
#          Colorblind-safe (deuteranopia/protanopia verified).
#          Works in grayscale (B&W printing).
#
# Usage:
#   source("scripts/00_colors.R")
#   geom_line(colour = COL$primary)
#
# PROJECT: Core-PAM (Memorial v6.1 / Freeze Core-PAM)
# =============================================================================

COL <- list(

  # ── Main analysis ─────────────────────────────────────────────────────────
  primary    = "#2166AC",   # Dark blue    — primary data curve/line
  selected   = "#D73027",   # Red          — selected panel / main result
  reference  = "#636363",   # Mid grey     — reference points (C_max, etc.)
  threshold  = "#D73027",   # Red          — threshold lines (same as selected)
  fill_warn  = "#FDAE61",   # Soft orange  — below-threshold zone (fill, reduced alpha)
  fill_ok    = "#ABD9E9",   # Soft blue    — above-threshold zone (fill, if needed)

  # ── Survival / KM ───────────────────────────────────────────────────────
  km_high    = "#D73027",   # Red   — high risk (worse prognosis)
  km_low     = "#1A9850",   # Green — low risk (better prognosis)
  km_q1      = "#1A9850",   # Q1 (best)
  km_q2      = "#91CF60",   # Q2
  km_q3      = "#FC8D59",   # Q3
  km_q4      = "#D73027",   # Q4 (worst)

  # ── OS cohorts ───────────────────────────────────────────────────────────
  scanb      = "#2166AC",   # Blue   — SCAN-B (training)
  tcga       = "#1A9850",   # Green  — TCGA-BRCA
  metabric   = "#D73027",   # Red    — METABRIC
  gse20685   = "#756BB1",   # Purple — GSE20685
  summary    = "#252525",   # Black  — meta-analysis diamond / summary

  # ── pCR cohorts (NACT block) ─────────────────────────────────────────────
  gse25066   = "#E6550D",   # Dark orange  — GSE25066 (Hatzis 2011, N=508)
  gse20194   = "#31A354",   # Mid green    — GSE20194 (Mina, N=278)
  gse32646   = "#3182BD",   # Mid blue     — GSE32646 (Tabchy, N=154)
  ispy1      = "#756BB1",   # Purple       — I-SPY1 / GSE22226 (N=149)
  ispy2      = "#BCBD22",   # Olive        — I-SPY2 (conditional)

  # ── PAM50 subtypes ──────────────────────────────────────────────────────
  lumA       = "#2166AC",   # Blue
  lumB       = "#74ADD1",   # Light blue
  her2       = "#D73027",   # Red
  basal      = "#FDAE61",   # Orange
  normal     = "#1A9850",   # Green

  # ── Utilities ────────────────────────────────────────────────────────────
  grey_light = "#F7F7F7",
  grey_mid   = "#BDBDBD",
  grey_dark  = "#636363",
  black      = "#252525",
  white      = "#FFFFFF"
)

# Named vectors for direct use in scale_colour_manual / scale_fill_manual
COL_COHORTS <- c(
  "SCANB"    = COL$scanb,
  "TCGA"     = COL$tcga,
  "METABRIC" = COL$metabric,
  "GSE20685" = COL$gse20685,
  "Summary"  = COL$summary
)

COL_RISK <- c("High" = COL$km_high, "Low" = COL$km_low)
