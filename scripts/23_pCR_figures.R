# =============================================================================
# SCRIPT: 23_pCR_figures.R
# PURPOSE: Generate all figures for the pCR (NACT response) analysis block:
#
#   Fig6:  Forest plot of CorePAM OR across pCR cohorts (main, bilingual)
#   FigS9: ROC curves per cohort (supp, bilingual)
#   FigS10: Score distribution by pCR status per cohort (supp, bilingual)
#
# INPUTS:
#   results/pcr/pCR_results_by_cohort.csv
#   results/pcr/meta_pCR_results.csv
#   results/pcr/meta_pCR_cohort_weights.csv
#   01_Base_Pura_CorePAM/PROCESSED/pCR/{COHORT}/analysis_ready.parquet
#
# PROJETO: Core-PAM (Memorial v6.1 / Freeze Core-PAM — pCR block)
# =============================================================================

source("scripts/00_setup.R")
source("scripts/00_colors.R")
SCRIPT_NAME <- "23_pCR_figures.R"

suppressPackageStartupMessages({
  library(ggplot2)
  library(pROC)
})

FORCE <- as.logical(Sys.getenv("FORCE_RERUN", "FALSE"))

fig_main <- PATHS$figures$main
fig_supp <- PATHS$figures$supp
dir.create(fig_main, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_supp, showWarnings = FALSE, recursive = TRUE)

# --------------------------------------------------------------------------
# Load results tables
# --------------------------------------------------------------------------
pcr_res_path    <- file.path(PATHS$results$pcr, "pCR_results_by_cohort.csv")
meta_res_path   <- file.path(PATHS$results$pcr, "meta_pCR_results.csv")
cohort_wt_path  <- file.path(PATHS$results$pcr, "meta_pCR_cohort_weights.csv")

for (p in c(pcr_res_path, meta_res_path, cohort_wt_path)) {
  if (!file.exists(p)) stop(sprintf("[%s] Required file missing: %s", SCRIPT_NAME, p))
}

pcr_df     <- strict_csv(pcr_res_path)
meta_df    <- strict_csv(meta_res_path)
cohort_wt  <- strict_csv(cohort_wt_path)

PCR_COHORTS <- pcr_df$cohort

message(sprintf("[%s] Loaded results for %d cohorts", SCRIPT_NAME, length(PCR_COHORTS)))

# RE pooled estimates
re_row    <- meta_df[meta_df$model == "RE", ]
re_logOR  <- re_row$pooled_logOR
re_se     <- re_row$pooled_se_logOR
re_OR     <- re_row$pooled_OR
re_lo     <- re_row$pooled_OR_lo95
re_hi     <- re_row$pooled_OR_hi95
re_p      <- re_row$p_pooled
re_I2     <- re_row$I2_pct
re_tau2   <- re_row$tau2
re_p_het  <- re_row$p_heterogeneity

# FE pooled estimates
fe_row  <- meta_df[meta_df$model == "FE", ]
fe_OR   <- fe_row$pooled_OR
fe_lo   <- fe_row$pooled_OR_lo95
fe_hi   <- fe_row$pooled_OR_hi95

# --------------------------------------------------------------------------
# pCR cohort colour map
# --------------------------------------------------------------------------
COL_PCR <- c(
  "GSE25066" = COL$gse25066,
  "GSE20194" = COL$gse20194,
  "GSE32646" = COL$gse32646,
  "ISPY1"    = COL$ispy1
)

# --------------------------------------------------------------------------
# FIGURE 6: Forest plot of ORs — bilingual
# --------------------------------------------------------------------------
make_forest_pcr <- function(lang = "EN") {
  xl    <- if (lang == "EN") "Odds Ratio (95% CI) per 1 SD CorePAM score_z" else
                              "Odds Ratio (IC 95%) por 1 DP do escore_z CorePAM"
  ttl   <- if (lang == "EN") "CorePAM score predicts pCR — meta-analysis" else
                              "Escore CorePAM prediz pCR — meta-análise"
  subt  <- if (lang == "EN")
    sprintf("Random-effects (DerSimonian-Laird) | I²=%.1f%% | τ²=%.4f | p_het=%.3g",
            re_I2, re_tau2, re_p_het)
  else
    sprintf("Efeitos aleatórios (DerSimonian-Laird) | I²=%.1f%% | τ²=%.4f | p_het=%.3g",
            re_I2, re_tau2, re_p_het)
  re_lbl <- if (lang == "EN") "Pooled RE" else "Combinado AE"
  fe_lbl <- if (lang == "EN") "Pooled FE" else "Combinado EF"

  # Individual cohort rows
  cohort_rows <- cohort_wt |>
    dplyr::mutate(
      label    = sprintf("%s\nN=%d, pCR=%d (%.0f%%)", cohort, n_total, n_pcr1,
                          100 * n_pcr1 / n_total),
      y_pos    = rev(seq_len(nrow(cohort_wt))),
      color_cat = cohort,
      is_pooled_re = FALSE,
      is_pooled_fe = FALSE,
      size_pt  = weight_re_pct / 10 + 1.5   # proportional to weight
    )

  # RE pooled row
  pooled_row_re <- tibble(
    cohort = "Pooled RE", label = re_lbl, y_pos = 0,
    or_uni = re_OR, or_lo95 = re_lo, or_hi95 = re_hi,
    n_total = NA, n_pcr1 = NA, p_uni = re_p, weight_re_pct = NA,
    color_cat = "Pooled RE", is_pooled_re = TRUE, is_pooled_fe = FALSE,
    size_pt = 5
  )

  # FE pooled row
  pooled_row_fe <- tibble(
    cohort = "Pooled FE", label = fe_lbl, y_pos = -1,
    or_uni = fe_OR, or_lo95 = fe_lo, or_hi95 = fe_hi,
    n_total = NA, n_pcr1 = NA, p_uni = fe_row$p_pooled, weight_re_pct = NA,
    color_cat = "Pooled FE", is_pooled_re = FALSE, is_pooled_fe = TRUE,
    size_pt = 5
  )

  forest_all <- dplyr::bind_rows(cohort_rows, pooled_row_re, pooled_row_fe)

  color_map <- c(COL_PCR,
                 "Pooled RE" = COL$black,
                 "Pooled FE" = COL$grey_dark)

  ggplot(forest_all, aes(x = or_uni, y = y_pos, colour = color_cat)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = COL$grey_mid, linewidth = 0.6) +
    geom_errorbarh(aes(xmin = or_lo95, xmax = or_hi95),
                   height = 0.25, linewidth = 0.7) +
    geom_point(aes(size = ifelse(is_pooled_re | is_pooled_fe, 5, size_pt),
                   shape = ifelse(is_pooled_re, "diamond", "circle")),
               stroke = 1.1) +
    scale_shape_identity() +
    scale_size_identity() +
    scale_colour_manual(values = color_map, guide = "none") +
    scale_x_log10(breaks = c(0.5, 0.75, 1, 1.25, 1.5, 2, 3),
                  labels = c("0.5", "0.75", "1", "1.25", "1.5", "2", "3")) +
    scale_y_continuous(
      breaks = forest_all$y_pos,
      labels = forest_all$label
    ) +
    geom_hline(yintercept = 0.5, linetype = "solid", colour = COL$grey_mid,
               linewidth = 0.4) +
    labs(title = ttl, subtitle = subt, x = xl, y = NULL,
         caption = sprintf(
           "OR per 1 SD score_z | Bootstrap B=1000 | %s",
           if (lang == "EN")
             sprintf("RE: OR=%.3f (%.3f–%.3f) p=%.3g | FE: OR=%.3f (%.3f–%.3f)",
                     re_OR, re_lo, re_hi, re_p, fe_OR, fe_lo, fe_hi)
           else
             sprintf("AE: OR=%.3f (%.3f–%.3f) p=%.3g | EF: OR=%.3f (%.3f–%.3f)",
                     re_OR, re_lo, re_hi, re_p, fe_OR, fe_lo, fe_hi)
         )) +
    theme_classic(base_size = 11) +
    theme(axis.text.y = element_text(size = 9),
          plot.caption = element_text(size = 7, colour = COL$grey_dark))
}

old_warn <- getOption("warn"); options(warn = 0)
for (lang in c("EN", "PT")) {
  p6 <- make_forest_pcr(lang)
  pdf_f <- file.path(fig_main, sprintf("Fig6_Forest_pCR_OR_%s.pdf", lang))
  png_f <- file.path(fig_main, sprintf("Fig6_Forest_pCR_OR_%s.png", lang))
  cairo_pdf(pdf_f, width = 9, height = 5); print(p6); dev.off()
  png(png_f, width = 9, height = 5, units = "in", res = 600); print(p6); dev.off()
  registry_append("META_PCR", sprintf("fig6_forest_pcr_%s", lang), pdf_f,
                  sha256_file(pdf_f), "ok", SCRIPT_NAME,
                  file.info(pdf_f)$size / 1e6)
  message(sprintf("[%s] [%s] Fig6 Forest saved: %s", SCRIPT_NAME, lang, pdf_f))
}
gc(); options(warn = old_warn)

# --------------------------------------------------------------------------
# FIGURE S9: ROC curves per cohort — bilingual
# --------------------------------------------------------------------------
# Load analysis_ready for each cohort to compute ROC
roc_list <- list()
for (cohort in PCR_COHORTS) {
  rp <- file.path(proc_pcr_cohort(cohort), "analysis_ready.parquet")
  if (!file.exists(rp)) next
  df <- strict_parquet(rp)
  df <- df[!is.na(df$pcr) & !is.na(df$score_z), ]
  old_warn2 <- getOption("warn"); options(warn = 0)
  roc_obj <- tryCatch(
    pROC::roc(df$pcr, df$score_z, direction = "<", quiet = TRUE),
    error = function(e) NULL
  )
  options(warn = old_warn2)
  if (!is.null(roc_obj)) {
    auc_val <- as.numeric(pROC::auc(roc_obj))
    roc_list[[cohort]] <- data.frame(
      cohort      = cohort,
      specificity = roc_obj$specificities,
      sensitivity = roc_obj$sensitivities,
      auc         = auc_val,
      stringsAsFactors = FALSE
    )
  }
}

if (length(roc_list) > 0) {
  roc_df <- bind_rows(roc_list)
  roc_df$auc_label <- sprintf("%s\nAUC=%.3f", roc_df$cohort, roc_df$auc)

  make_roc_fig <- function(lang = "EN") {
    xl <- if (lang == "EN") "1 – Specificity" else "1 – Especificidade"
    yl <- if (lang == "EN") "Sensitivity" else "Sensibilidade"
    tt <- if (lang == "EN") "CorePAM ROC — pCR prediction (NACT cohorts)" else
                             "CorePAM ROC — predição de pCR (coortes NACT)"
    ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity,
                        colour = cohort, group = cohort)) +
      geom_line(linewidth = 0.9) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  colour = COL$grey_mid, linewidth = 0.5) +
      scale_colour_manual(
        values  = COL_PCR,
        labels  = setNames(
          sprintf("%s (AUC=%.3f)", names(COL_PCR),
                  sapply(names(COL_PCR), function(c) {
                    v <- unique(roc_df$auc[roc_df$cohort == c]); if (length(v)) v[1] else NA
                  })),
          names(COL_PCR)
        ),
        name = if (lang == "EN") "Cohort" else "Coorte"
      ) +
      labs(title = tt, x = xl, y = yl) +
      theme_classic(base_size = 11) +
      theme(legend.text = element_text(size = 9))
  }

  old_warn <- getOption("warn"); options(warn = 0)
  for (lang in c("EN", "PT")) {
    ps9 <- make_roc_fig(lang)
    pdf_f <- file.path(fig_supp, sprintf("FigS9_ROC_pCR_%s.pdf", lang))
    png_f <- file.path(fig_supp, sprintf("FigS9_ROC_pCR_%s.png", lang))
    cairo_pdf(pdf_f, width = 7, height = 6); print(ps9); dev.off()
    png(png_f, width = 7, height = 6, units = "in", res = 600); print(ps9); dev.off()
    registry_append("META_PCR", sprintf("figS9_roc_%s", lang), pdf_f,
                    sha256_file(pdf_f), "ok", SCRIPT_NAME,
                    file.info(pdf_f)$size / 1e6)
    message(sprintf("[%s] [%s] FigS9 ROC saved: %s", SCRIPT_NAME, lang, pdf_f))
  }
  gc(); options(warn = old_warn)
}

# --------------------------------------------------------------------------
# FIGURE S10: Score distribution by pCR status per cohort — bilingual
# --------------------------------------------------------------------------
dist_list <- list()
for (cohort in PCR_COHORTS) {
  rp <- file.path(proc_pcr_cohort(cohort), "analysis_ready.parquet")
  if (!file.exists(rp)) next
  df <- strict_parquet(rp)
  df <- df[!is.na(df$pcr) & !is.na(df$score_z), ]
  df$cohort     <- cohort
  df$pcr_label  <- ifelse(df$pcr == 1, "pCR", "RD")
  dist_list[[cohort]] <- df[, c("cohort", "score_z", "pcr_label")]
}

if (length(dist_list) > 0) {
  dist_df <- bind_rows(dist_list)

  make_dist_fig <- function(lang = "EN") {
    xl <- if (lang == "EN") "CorePAM score_z (intra-cohort Z-score)" else
                             "Escore CorePAM score_z (Z-score intracoorte)"
    yl <- if (lang == "EN") "Density" else "Densidade"
    tt <- if (lang == "EN") "CorePAM score distribution by pCR status" else
                             "Distribuição do escore CorePAM por status de pCR"
    pcr_lbl <- if (lang == "EN") "pCR status" else "Status pCR"
    ggplot(dist_df, aes(x = score_z, fill = pcr_label, colour = pcr_label)) +
      geom_density(alpha = 0.35, linewidth = 0.8) +
      facet_wrap(~ cohort, scales = "free_y", ncol = 2) +
      scale_fill_manual(values = c("pCR" = COL$km_low, "RD" = COL$km_high),
                        name = pcr_lbl) +
      scale_colour_manual(values = c("pCR" = COL$km_low, "RD" = COL$km_high),
                          guide = "none") +
      labs(title = tt, x = xl, y = yl) +
      theme_classic(base_size = 11) +
      theme(strip.text = element_text(size = 10, face = "bold"))
  }

  old_warn <- getOption("warn"); options(warn = 0)
  for (lang in c("EN", "PT")) {
    ps10 <- make_dist_fig(lang)
    pdf_f <- file.path(fig_supp, sprintf("FigS10_ScoreDist_pCR_%s.pdf", lang))
    png_f <- file.path(fig_supp, sprintf("FigS10_ScoreDist_pCR_%s.png", lang))
    cairo_pdf(pdf_f, width = 8, height = 6); print(ps10); dev.off()
    png(png_f, width = 8, height = 6, units = "in", res = 600); print(ps10); dev.off()
    registry_append("META_PCR", sprintf("figS10_dist_%s", lang), pdf_f,
                    sha256_file(pdf_f), "ok", SCRIPT_NAME,
                    file.info(pdf_f)$size / 1e6)
    message(sprintf("[%s] [%s] FigS10 Dist saved: %s", SCRIPT_NAME, lang, pdf_f))
  }
  gc(); options(warn = old_warn)
}

message(sprintf("[%s] COMPLETED — pCR figures: Fig6, FigS9, FigS10", SCRIPT_NAME))
