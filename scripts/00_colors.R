# =============================================================================
# 00_colors.R — Paleta de cores padronizada do projeto Core-PAM
#
# REGRA: Todas as figuras do projeto devem usar EXCLUSIVAMENTE estas cores.
# Paleta base: ColorBrewer RdBu (divergente) + Set1 qualitativo (coortes).
# Colorblind-safe (deuteranopia/protanopia verificado).
# Funciona em escala de cinza (impressão P&B).
#
# Como usar:
#   source("scripts/00_colors.R")
#   geom_line(colour = COL$primary)
# =============================================================================

COL <- list(

  # ── Análise principal ──────────────────────────────────────────────────────
  primary    = "#2166AC",   # Azul escuro  — curva/linha principal de dados
  selected   = "#D73027",   # Vermelho     — painel selecionado / resultado principal
  reference  = "#636363",   # Cinza médio  — pontos de referência (C_max, etc.)
  threshold  = "#D73027",   # Vermelho     — linhas de limiar (mesmo que selected)
  fill_warn  = "#FDAE61",   # Laranja suave — zona abaixo do limiar (fill, alpha reduzido)
  fill_ok    = "#ABD9E9",   # Azul suave   — zona acima do limiar (fill, se necessário)

  # ── Sobrevida / KM ────────────────────────────────────────────────────────
  km_high    = "#D73027",   # Vermelho — alto risco (pior prognóstico)
  km_low     = "#1A9850",   # Verde    — baixo risco (melhor prognóstico)
  km_q1      = "#1A9850",   # Q1 (melhor)
  km_q2      = "#91CF60",   # Q2
  km_q3      = "#FC8D59",   # Q3
  km_q4      = "#D73027",   # Q4 (pior)

  # ── Coortes ───────────────────────────────────────────────────────────────
  scanb      = "#2166AC",   # Azul     — SCAN-B (treino)
  tcga       = "#1A9850",   # Verde    — TCGA-BRCA
  metabric   = "#D73027",   # Vermelho — METABRIC
  gse20685   = "#756BB1",   # Roxo     — GSE20685
  summary    = "#252525",   # Preto    — diamante de meta-análise / resumo

  # ── Subtipos PAM50 ────────────────────────────────────────────────────────
  lumA       = "#2166AC",   # Azul
  lumB       = "#74ADD1",   # Azul claro
  her2       = "#D73027",   # Vermelho
  basal      = "#FDAE61",   # Laranja
  normal     = "#1A9850",   # Verde

  # ── Utilitários ───────────────────────────────────────────────────────────
  grey_light = "#F7F7F7",
  grey_mid   = "#BDBDBD",
  grey_dark  = "#636363",
  black      = "#252525",
  white      = "#FFFFFF"
)

# Vetores nomeados para uso direto em scale_colour_manual / scale_fill_manual
COL_COHORTS <- c(
  "SCANB"    = COL$scanb,
  "TCGA"     = COL$tcga,
  "METABRIC" = COL$metabric,
  "GSE20685" = COL$gse20685,
  "Summary"  = COL$summary
)

COL_RISK <- c("High" = COL$km_high, "Low" = COL$km_low)
