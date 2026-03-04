# =============================================================================
# SCRIPT: 23_pCR_figures.R
# PURPOSE: Generate all figures for the pCR (NACT response) analysis block.
#          All figures saved to figures/pcr/ (separate from OS block).
#
# PRIMARY pCR COHORTS (frozen, Memorial §11):
#   GSE25066 (Hatzis 2011, N=508)
#   GSE20194 (Park 2006, N=278)
#   GSE32646 (Tabchy 2010, N=154)
#   ISPY1 / GSE22226 (Esserman 2012, N=149)
#
# FIGURES PRODUCED:
#   Fig_pCR1:  Forest plot of ORs — main figure, bilingual
#   Fig_pCR2:  ROC curves per cohort (small multiples) — supp, bilingual
#   Fig_pCR3:  pCR rate by score quartile per cohort — supp, bilingual
#   Fig_pCR4:  Score distribution by pCR status — supp, bilingual
#
# INPUTS:
#   results/pcr/pCR_results_by_cohort.csv
#   results/pcr/meta_pCR_results.csv
#   results/pcr/meta_pCR_cohort_weights.csv
#   01_Base_Pura_CorePAM/PROCESSED/pCR/{COHORT}/analysis_ready.parquet
#
# PROJECT: Core-PAM (Memorial v6.1 / Freeze Core-PAM — pCR block)
# =============================================================================

source("scripts/00_setup.R")
source("scripts/00_colors.R")
SCRIPT_NAME <- "23_pCR_figures.R"

suppressPackageStartupMessages({
  library(ggplot2)
  library(pROC)
  library(dplyr)
})

FORCE <- as.logical(Sys.getenv("FORCE_RERUN", "FALSE"))

# pcr_fig_path() — returns the correct path for a pCR figure given language and extension
# Replaces the flat fig_pcr variable used in earlier versions
pcr_fig_path <- function(name, lang, ext) {
  dir <- PATHS$figures[[paste0("pcr_", tolower(lang), "_", ext)]]
  file.path(dir, name)
}

# --------------------------------------------------------------------------
# Load results tables
# --------------------------------------------------------------------------
pcr_res_path   <- file.path(PATHS$results$pcr, "pCR_results_by_cohort.csv")
meta_res_path  <- file.path(PATHS$results$pcr, "meta_pCR_results.csv")
cohort_wt_path <- file.path(PATHS$results$pcr, "meta_pCR_cohort_weights.csv")

for (p in c(pcr_res_path, meta_res_path, cohort_wt_path)) {
  if (!file.exists(p)) stop(sprintf("[%s] Required file missing: %s", SCRIPT_NAME, p))
}

pcr_df    <- strict_csv(pcr_res_path)
meta_df   <- strict_csv(meta_res_path)
cohort_wt <- strict_csv(cohort_wt_path)

PCR_COHORTS <- pcr_df$cohort
message(sprintf("[%s] Loaded results for %d cohorts", SCRIPT_NAME, length(PCR_COHORTS)))

# RE pooled estimates
re_row   <- meta_df[meta_df$model == "RE", ]
re_logOR <- re_row$pooled_logOR
re_se    <- re_row$pooled_se_logOR
re_OR    <- re_row$pooled_OR
re_lo    <- re_row$pooled_OR_lo95
re_hi    <- re_row$pooled_OR_hi95
re_p     <- re_row$p_pooled
re_I2    <- re_row$I2_pct
re_tau2  <- re_row$tau2
re_p_het <- re_row$p_heterogeneity

# FE pooled estimates
fe_row <- meta_df[meta_df$model == "FE", ]
fe_OR  <- fe_row$pooled_OR
fe_lo  <- fe_row$pooled_OR_lo95
fe_hi  <- fe_row$pooled_OR_hi95

# --------------------------------------------------------------------------
# pCR cohort colour map (primary cohorts only)
# --------------------------------------------------------------------------
COL_PCR <- c(
  "GSE25066" = COL$gse25066,
  "GSE20194" = COL$gse20194,
  "GSE32646" = COL$gse32646,
  "ISPY1"    = COL$ispy1
)

# --------------------------------------------------------------------------
# FIG_pCR1: Forest plot of ORs — bilingual
# --------------------------------------------------------------------------
make_forest_pcr <- function(lang = "EN") {
  xl  <- if (lang == "EN") "Odds Ratio (95% CI) per 1 SD CorePAM score_z" else
                            "Odds Ratio (IC 95%) por 1 DP do escore_z CorePAM"
  ttl <- if (lang == "EN") "CorePAM score predicts pCR — meta-analysis" else
                            "Escore CorePAM prediz pCR — meta-análise"
  subt <- if (lang == "EN")
    sprintf("Random-effects (DerSimonian-Laird) | I\u00b2=%.1f%% | \u03c4\u00b2=%.4f | p_het=%.3g",
            re_I2, re_tau2, re_p_het)
  else
    sprintf("Efeitos aleat\u00f3rios (DerSimonian-Laird) | I\u00b2=%.1f%% | \u03c4\u00b2=%.4f | p_het=%.3g",
            re_I2, re_tau2, re_p_het)
  re_lbl <- if (lang == "EN") "Pooled RE" else "Combinado AE"
  fe_lbl <- if (lang == "EN") "Pooled FE" else "Combinado EF"

  cohort_rows <- cohort_wt |>
    dplyr::mutate(
      label        = sprintf("%s  N=%d, pCR=%d (%.0f%%)",
                             cohort, n_total, n_pcr1, 100 * n_pcr1 / n_total),
      or_label     = sprintf("%.2f (%.2f\u2013%.2f)",
                              or_uni, or_lo95, or_hi95),
      y_pos        = rev(seq_len(nrow(cohort_wt))),
      color_cat    = cohort,
      is_pooled_re = FALSE,
      is_pooled_fe = FALSE,
      size_pt      = weight_re_pct / 10 + 1.5
    )

  pooled_row_re <- tibble(
    cohort = "Pooled RE", label = re_lbl, y_pos = 0,
    or_uni = re_OR, or_lo95 = re_lo, or_hi95 = re_hi,
    or_label = sprintf("%.2f (%.2f\u2013%.2f)", re_OR, re_lo, re_hi),
    n_total = NA, n_pcr1 = NA, p_uni = re_p, weight_re_pct = NA,
    color_cat = "Pooled RE", is_pooled_re = TRUE, is_pooled_fe = FALSE,
    size_pt = 5
  )
  pooled_row_fe <- tibble(
    cohort = "Pooled FE", label = fe_lbl, y_pos = -1,
    or_uni = fe_OR, or_lo95 = fe_lo, or_hi95 = fe_hi,
    or_label = sprintf("%.2f (%.2f\u2013%.2f)", fe_OR, fe_lo, fe_hi),
    n_total = NA, n_pcr1 = NA, p_uni = fe_row$p_pooled, weight_re_pct = NA,
    color_cat = "Pooled FE", is_pooled_re = FALSE, is_pooled_fe = TRUE,
    size_pt = 5
  )

  forest_all <- dplyr::bind_rows(cohort_rows, pooled_row_re, pooled_row_fe)

  color_map <- c(COL_PCR,
                 "Pooled RE" = COL$black,
                 "Pooled FE" = COL$grey_dark)

  ggplot(forest_all, aes(x = or_uni, y = y_pos, colour = color_cat)) +
    geom_vline(xintercept = 1, linetype = "dashed",
               colour = COL$grey_mid, linewidth = 0.6) +
    geom_errorbarh(aes(xmin = or_lo95, xmax = or_hi95),
                   height = 0.25, linewidth = 0.7) +
    geom_point(aes(size = ifelse(is_pooled_re | is_pooled_fe, 5, size_pt),
                   shape = ifelse(is_pooled_re, "diamond", "circle")),
               stroke = 1.1) +
    geom_text(aes(x = max(forest_all$or_hi95, na.rm = TRUE) * 1.15,
                  label = or_label),
              hjust = 0, size = 2.8, color = "grey30", show.legend = FALSE) +
    scale_shape_identity() +
    scale_size_identity() +
    scale_colour_manual(values = color_map, guide = "none") +
    scale_x_log10(breaks = c(0.5, 0.75, 1, 1.25, 1.5, 2, 3),
                  labels = c("0.5", "0.75", "1", "1.25", "1.5", "2", "3")) +
    coord_cartesian(xlim = c(0.3, max(forest_all$or_hi95, na.rm = TRUE) * 2.5)) +
    scale_y_continuous(breaks = forest_all$y_pos, labels = forest_all$label) +
    geom_hline(yintercept = 0.5, linetype = "solid",
               colour = COL$grey_mid, linewidth = 0.4) +
    labs(title = ttl, subtitle = subt, x = xl, y = NULL,
         caption = sprintf(
           "OR per 1 SD score_z | Bootstrap B=1000 | %s",
           if (lang == "EN")
             sprintf("RE: OR=%.3f (%.3f\u2013%.3f) p=%.3g | FE: OR=%.3f (%.3f\u2013%.3f)",
                     re_OR, re_lo, re_hi, re_p, fe_OR, fe_lo, fe_hi)
           else
             sprintf("AE: OR=%.3f (%.3f\u2013%.3f) p=%.3g | EF: OR=%.3f (%.3f\u2013%.3f)",
                     re_OR, re_lo, re_hi, re_p, fe_OR, fe_lo, fe_hi)
         )) +
    theme_classic(base_size = 11) +
    theme(axis.text.y  = element_text(size = 9),
          plot.caption = element_text(size = 7, colour = COL$grey_dark))
}

old_warn <- getOption("warn"); options(warn = 0)
for (lang in c("EN", "PT")) {
  p1 <- make_forest_pcr(lang)
  pdf_f <- pcr_fig_path(sprintf("Fig_pCR1_Forest_OR_%s.pdf", lang), lang, "pdf")
  png_f <- pcr_fig_path(sprintf("Fig_pCR1_Forest_OR_%s.png", lang), lang, "png")
  cairo_pdf(pdf_f, width = 9, height = 5); print(p1); dev.off()
  png(png_f, width = 9, height = 5, units = "in", res = 600); print(p1); dev.off()
  registry_append("META_PCR", sprintf("fig_pcr1_forest_%s", lang), pdf_f,
                  sha256_file(pdf_f), "ok", SCRIPT_NAME,
                  file.info(pdf_f)$size / 1e6)
  message(sprintf("[%s] [%s] Fig_pCR1 Forest saved: %s", SCRIPT_NAME, lang, pdf_f))
}
gc(); options(warn = old_warn)

# --------------------------------------------------------------------------
# Load per-sample data for figures 2-4
# --------------------------------------------------------------------------
sample_list <- list()
for (cohort in PCR_COHORTS) {
  rp <- file.path(proc_pcr_cohort(cohort), "analysis_ready.parquet")
  if (!file.exists(rp)) next
  df <- strict_parquet(rp)
  df <- df[!is.na(df$pcr) & !is.na(df$score_z), ]
  df$cohort <- cohort
  sample_list[[cohort]] <- df
}

# --------------------------------------------------------------------------
# FIG_pCR2: ROC curves per cohort — faceted small multiples, bilingual
# --------------------------------------------------------------------------
roc_list <- list()
for (cohort in names(sample_list)) {
  df <- sample_list[[cohort]]
  old_warn2 <- getOption("warn"); options(warn = 0)
  roc_obj <- tryCatch(
    pROC::roc(df$pcr, df$score_z, direction = "<", quiet = TRUE),
    error = function(e) NULL
  )
  options(warn = old_warn2)
  if (!is.null(roc_obj)) {
    auc_val <- as.numeric(pROC::auc(roc_obj))
    auc_ci  <- as.numeric(pROC::ci.auc(roc_obj, conf.level = 0.95))
    roc_list[[cohort]] <- data.frame(
      cohort      = sprintf("%s\nAUC=%.3f (%.3f\u2013%.3f)", cohort,
                            auc_val, auc_ci[1], auc_ci[3]),
      cohort_raw  = cohort,
      specificity = roc_obj$specificities,
      sensitivity = roc_obj$sensitivities,
      auc         = auc_val,
      stringsAsFactors = FALSE
    )
  }
}

if (length(roc_list) > 0) {
  roc_df <- bind_rows(roc_list)

  make_roc_fig <- function(lang = "EN") {
    xl  <- if (lang == "EN") "1 \u2013 Specificity" else "1 \u2013 Especificidade"
    yl  <- if (lang == "EN") "Sensitivity"          else "Sensibilidade"
    tt  <- if (lang == "EN") "CorePAM ROC \u2014 pCR prediction (NACT cohorts)" else
                              "CorePAM ROC \u2014 predi\u00e7\u00e3o de pCR (coortes NACT)"
    ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity,
                       colour = cohort_raw, group = cohort)) +
      geom_line(linewidth = 0.9) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  colour = COL$grey_mid, linewidth = 0.5) +
      facet_wrap(~ cohort, ncol = 2) +
      scale_colour_manual(values = COL_PCR, guide = "none") +
      labs(title = tt, x = xl, y = yl) +
      theme_classic(base_size = 11) +
      theme(strip.text = element_text(size = 9, face = "bold"),
            strip.background = element_blank())
  }

  old_warn <- getOption("warn"); options(warn = 0)
  for (lang in c("EN", "PT")) {
    p2 <- make_roc_fig(lang)
    pdf_f <- pcr_fig_path(sprintf("Fig_pCR2_ROC_%s.pdf", lang), lang, "pdf")
    png_f <- pcr_fig_path(sprintf("Fig_pCR2_ROC_%s.png", lang), lang, "png")
    cairo_pdf(pdf_f, width = 8, height = 7); print(p2); dev.off()
    png(png_f, width = 8, height = 7, units = "in", res = 600); print(p2); dev.off()
    registry_append("META_PCR", sprintf("fig_pcr2_roc_%s", lang), pdf_f,
                    sha256_file(pdf_f), "ok", SCRIPT_NAME,
                    file.info(pdf_f)$size / 1e6)
    message(sprintf("[%s] [%s] Fig_pCR2 ROC saved: %s", SCRIPT_NAME, lang, pdf_f))
  }
  gc(); options(warn = old_warn)
}

# --------------------------------------------------------------------------
# FIG_pCR3: pCR rate by score quartile per cohort — bilingual
# --------------------------------------------------------------------------
quartile_list <- list()
for (cohort in names(sample_list)) {
  df <- sample_list[[cohort]]
  df$quartile <- cut(df$score_z, breaks = quantile(df$score_z, probs = 0:4/4,
                                                    na.rm = TRUE),
                     labels = c("Q1", "Q2", "Q3", "Q4"), include.lowest = TRUE)
  df_q <- df[!is.na(df$quartile), ]
  qdf <- do.call(rbind, lapply(levels(df_q$quartile), function(q) {
    sub_df <- df_q[df_q$quartile == q, ]
    {
      n_q   <- nrow(sub_df)
      k_q   <- sum(sub_df$pcr == 1L, na.rm = TRUE)
      # Wilson score interval (asymmetric, better for proportions)
      z95   <- 1.96
      denom <- n_q + z95^2
      p_hat <- k_q / n_q
      center <- (p_hat + z95^2 / (2 * n_q)) / (1 + z95^2 / n_q)
      half   <- z95 * sqrt(p_hat * (1 - p_hat) / n_q + z95^2 / (4 * n_q^2)) / (1 + z95^2 / n_q)
      data.frame(
        cohort   = cohort,
        quartile = q,
        pcr_rate = p_hat,
        n        = n_q,
        pcr1     = k_q,
        ci_lo    = pmax(0, center - half),
        ci_hi    = pmin(1, center + half),
        stringsAsFactors = FALSE
      )
    }
  }))
  quartile_list[[cohort]] <- qdf
}

if (length(quartile_list) > 0) {
  quart_df <- bind_rows(quartile_list)

  make_quartile_fig <- function(lang = "EN") {
    xl  <- if (lang == "EN") "CorePAM score quartile" else
                              "Quartil do escore CorePAM"
    yl  <- if (lang == "EN") "pCR rate" else "Taxa de pCR"
    tt  <- if (lang == "EN") "pCR rate by CorePAM score quartile" else
                              "Taxa de pCR por quartil do escore CorePAM"
    ggplot(quart_df, aes(x = quartile, y = pcr_rate, fill = cohort,
                         group = cohort)) +
      geom_col(position = "dodge", colour = COL$black, linewidth = 0.3, alpha = 0.85) +
      geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                    position = position_dodge(width = 0.9),
                    width = 0.25, linewidth = 0.6) +
      geom_text(aes(label = sprintf("%d/%d", pcr1, n)),
                position = position_dodge(width = 0.9),
                vjust = -0.5, size = 2.5) +
      scale_fill_manual(values = COL_PCR,
                        name = if (lang == "EN") "Cohort" else "Coorte") +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                         limits = c(0, 1)) +
      labs(title = tt, x = xl, y = yl) +
      theme_classic(base_size = 11) +
      theme(legend.position = "top")
  }

  old_warn <- getOption("warn"); options(warn = 0)
  for (lang in c("EN", "PT")) {
    p3 <- make_quartile_fig(lang)
    pdf_f <- pcr_fig_path(sprintf("Fig_pCR3_QuartileRate_%s.pdf", lang), lang, "pdf")
    png_f <- pcr_fig_path(sprintf("Fig_pCR3_QuartileRate_%s.png", lang), lang, "png")
    cairo_pdf(pdf_f, width = 8, height = 5); print(p3); dev.off()
    png(png_f, width = 8, height = 5, units = "in", res = 600); print(p3); dev.off()
    registry_append("META_PCR", sprintf("fig_pcr3_quartile_%s", lang), pdf_f,
                    sha256_file(pdf_f), "ok", SCRIPT_NAME,
                    file.info(pdf_f)$size / 1e6)
    message(sprintf("[%s] [%s] Fig_pCR3 QuartileRate saved: %s", SCRIPT_NAME, lang, pdf_f))
  }
  gc(); options(warn = old_warn)
}

# --------------------------------------------------------------------------
# FIG_pCR4: Score distribution by pCR status per cohort — bilingual
# --------------------------------------------------------------------------
if (length(sample_list) > 0) {
  dist_df <- bind_rows(lapply(names(sample_list), function(coh) {
    df <- sample_list[[coh]]
    df$pcr_label <- ifelse(df$pcr == 1, "pCR", "RD")
    df[, c("cohort", "score_z", "pcr_label")]
  }))

  make_dist_fig <- function(lang = "EN") {
    xl      <- if (lang == "EN") "CorePAM score_z (intra-cohort Z-score)" else
                                  "Escore CorePAM score_z (Z-score intracoorte)"
    yl      <- if (lang == "EN") "Density" else "Densidade"
    tt      <- if (lang == "EN") "CorePAM score distribution by pCR status" else
                                  "Distribui\u00e7\u00e3o do escore CorePAM por status de pCR"
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
      theme(strip.text = element_text(size = 10, face = "bold"),
            strip.background = element_blank())
  }

  old_warn <- getOption("warn"); options(warn = 0)
  for (lang in c("EN", "PT")) {
    p4 <- make_dist_fig(lang)
    pdf_f <- pcr_fig_path(sprintf("Fig_pCR4_ScoreDist_%s.pdf", lang), lang, "pdf")
    png_f <- pcr_fig_path(sprintf("Fig_pCR4_ScoreDist_%s.png", lang), lang, "png")
    cairo_pdf(pdf_f, width = 8, height = 6); print(p4); dev.off()
    png(png_f, width = 8, height = 6, units = "in", res = 600); print(p4); dev.off()
    registry_append("META_PCR", sprintf("fig_pcr4_dist_%s", lang), pdf_f,
                    sha256_file(pdf_f), "ok", SCRIPT_NAME,
                    file.info(pdf_f)$size / 1e6)
    message(sprintf("[%s] [%s] Fig_pCR4 ScoreDist saved: %s", SCRIPT_NAME, lang, pdf_f))
  }
  gc(); options(warn = old_warn)
}

# --------------------------------------------------------------------------
# FIG_pCR_S1: PR curves (precision-recall) — supplementary, bilingual
# --------------------------------------------------------------------------
pr_list <- list()
for (cohort in PCR_COHORTS) {
  rp <- file.path(proc_pcr_cohort(cohort), "analysis_ready.parquet")
  if (!file.exists(rp)) next
  df_pr <- strict_parquet(rp)
  df_pr <- df_pr[!is.na(df_pr$pcr) & !is.na(df_pr$score_z), ]
  old_warn2 <- getOption("warn"); options(warn = 0)
  fit_pr <- tryCatch(
    glm(pcr ~ score_z, data = df_pr, family = binomial),
    error = function(e) NULL,
    warning = function(w) glm(pcr ~ score_z, data = df_pr, family = binomial)
  )
  options(warn = old_warn2)
  if (is.null(fit_pr)) next
  probs_pr <- predict(fit_pr, type = "response")
  n_pos    <- sum(df_pr$pcr == 1L)
  base_rate <- n_pos / nrow(df_pr)
  # Precision-Recall curve
  thresh <- sort(unique(probs_pr), decreasing = TRUE)
  pr_pts <- lapply(thresh, function(t) {
    tp <- sum(probs_pr >= t & df_pr$pcr == 1L)
    pp <- sum(probs_pr >= t)
    c(prec = if (pp > 0) tp / pp else 1, rec = tp / n_pos)
  })
  pr_mat <- do.call(rbind, pr_pts)
  # Add boundary points
  pr_mat <- rbind(c(1, 0), pr_mat, c(base_rate, 1))
  # Trapezoid PR-AUC
  pr_auc <- sum(diff(pr_mat[, "rec"]) *
                (pr_mat[-nrow(pr_mat), "prec"] + pr_mat[-1, "prec"]) / 2)
  pr_list[[cohort]] <- data.frame(
    cohort    = sprintf("%s\nPR-AUC=%.3f", cohort, pr_auc),
    cohort_raw = cohort,
    precision = pr_mat[, "prec"],
    recall    = pr_mat[, "rec"],
    pr_auc    = pr_auc,
    base_rate = base_rate,
    stringsAsFactors = FALSE
  )
}

if (length(pr_list) > 0) {
  pr_df_all <- bind_rows(pr_list)

  make_pr_fig <- function(lang = "EN") {
    xl  <- if (lang == "EN") "Recall (Sensitivity)" else "Recall (Sensibilidade)"
    yl  <- if (lang == "EN") "Precision (PPV)"      else "Precisão (VPP)"
    tt  <- if (lang == "EN")
              "CorePAM precision-recall curve (pCR prediction)" else
              "CorePAM curva precisão-recall (predição de pCR)"
    ggplot(pr_df_all, aes(x = recall, y = precision,
                          colour = cohort_raw, group = cohort)) +
      geom_line(linewidth = 0.9) +
      geom_hline(aes(yintercept = base_rate, colour = cohort_raw),
                 linetype = "dashed", linewidth = 0.4, alpha = 0.6) +
      facet_wrap(~ cohort, ncol = 2) +
      scale_colour_manual(values = COL_PCR, guide = "none") +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      labs(title = tt, x = xl, y = yl,
           caption = if (lang == "EN")
             "Dashed line: prevalence (no-skill baseline)"
           else
             "Linha tracejada: prevalência (baseline sem discriminação)") +
      theme_classic(base_size = 11) +
      theme(strip.text = element_text(size = 9, face = "bold"),
            strip.background = element_blank())
  }

  old_warn <- getOption("warn"); options(warn = 0)
  for (lang in c("EN", "PT")) {
    ps1 <- make_pr_fig(lang)
    pdf_f <- pcr_fig_path(sprintf("FigS_pCR_PR_%s.pdf", lang), lang, "pdf")
    png_f <- pcr_fig_path(sprintf("FigS_pCR_PR_%s.png", lang), lang, "png")
    cairo_pdf(pdf_f, width = 8, height = 7); print(ps1); dev.off()
    png(png_f, width = 8, height = 7, units = "in", res = 600); print(ps1); dev.off()
    registry_append("META_PCR", sprintf("figS_pcr_pr_%s", lang), pdf_f,
                    sha256_file(pdf_f), "ok", SCRIPT_NAME,
                    file.info(pdf_f)$size / 1e6)
    message(sprintf("[%s] [%s] FigS_pCR PR-curve saved: %s", SCRIPT_NAME, lang, pdf_f))
  }
  gc(); options(warn = old_warn)
}

# --------------------------------------------------------------------------
# Artifact hash manifest
# --------------------------------------------------------------------------
fig_pcr_dir <- PATHS$figures$pcr_en_pdf
fig_files <- list.files(fig_pcr_dir, pattern = "Fig_pCR.*\\.(pdf|png)$",
                        full.names = TRUE)
if (length(fig_files) > 0) {
  manifest_df <- tibble(
    file     = basename(fig_files),
    path     = fig_files,
    sha256   = sapply(fig_files, sha256_file),
    size_mb  = round(file.info(fig_files)$size / 1e6, 3),
    script   = SCRIPT_NAME,
    created  = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  manifest_path <- file.path(PATHS$results$pcr, "artifact_hash_manifest.csv")
  readr::write_csv(manifest_df, manifest_path)
  message(sprintf("[%s] artifact_hash_manifest.csv saved (%d files)",
                  SCRIPT_NAME, nrow(manifest_df)))
}

message(sprintf("[%s] COMPLETED — pCR figures: Fig_pCR1-4 saved to %s",
                SCRIPT_NAME, fig_pcr_dir))
