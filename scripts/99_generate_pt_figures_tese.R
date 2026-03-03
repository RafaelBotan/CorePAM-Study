# =============================================================================
# SCRIPT: 99_generate_pt_figures_tese.R
# PURPOSE: Generate Portuguese (PT) versions of all figures that exist only in
#          English. For the doctoral thesis (tese).
#
# FIGURES GENERATED (15 total):
#   MAIN (figures/main/pt/png/):
#     1. Fig3_KM_SCANB_OS_CorePAM_PT.png       — KM SCAN-B training
#     2. Fig5_Calibration_60m_Panels_CorePAM_PT.png — Calibration plots
#     3. Fig5_DeltaCindex_COREA_vs_COREAplus_CorePAM_PT.png — Delta C-index
#
#   SUPP (figures/supp/pt/png/):
#     4.  FigS3_Correlation_OffDiagonal_PT.png
#     5.  FigS4_METABRIC_PCA_forensics_PT.png
#     6.  FigS_Bootstrap_GeneFreq_PT.png
#     7.  FigS_COREA_Sensitivity_PT.png
#     8.  FigS_Cindex_ByCohort_PT.png
#     9.  FigS_CorePAM_vs_PAM50full_HR_PT.png
#     10. FigS_DropGene_Tornado_PT.png
#     11. FigS_Forest_HR_ValidationCohorts_PT.png
#     12. FigS_Heatmap_GSE20685_PT.png
#     13. FigS_Heatmap_METABRIC_PT.png
#     14. FigS_Heatmap_TCGA_BRCA_PT.png
#     15. FigS_Weights_CorePAM_Lollipop_PT.png
#
# PROJETO: Core-PAM (Memorial v6.1 / Freeze Core-PAM)
# =============================================================================

source("scripts/00_setup.R")
options(warn = 0)
SCRIPT_NAME <- "99_generate_pt_figures_tese.R"

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(arrow)
  library(patchwork)
  library(ggrepel)
  library(glmnet)
})
source("scripts/00_colors.R")

set.seed(42)

message(sprintf("[%s] Starting PT figure generation for thesis", SCRIPT_NAME))

# --------------------------------------------------------------------------
# Helper: save PT figure (PDF + PNG)
# --------------------------------------------------------------------------
save_pt <- function(p, name, w = 8, h = 6, dpi = 300,
                    section = "supp", device_type = "ggplot") {
  pdf_dir <- PATHS$figures[[paste0(section, "_pt_pdf")]]
  png_dir <- PATHS$figures[[paste0(section, "_pt_png")]]
  dir.create(pdf_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(png_dir, showWarnings = FALSE, recursive = TRUE)
  pdf_path <- file.path(pdf_dir, paste0(name, ".pdf"))
  png_path <- file.path(png_dir, paste0(name, ".png"))
  tryCatch({
    if (device_type == "ggplot") {
      cairo_pdf(pdf_path, width = w, height = h); print(p); dev.off()
      png(png_path, width = round(w * dpi), height = round(h * dpi), res = dpi)
      print(p); dev.off()
    } else if (device_type == "base") {
      # For base R / grid-based plots (survminer, pheatmap, gridExtra)
      pdf(pdf_path, width = w, height = h); print(p); dev.off()
      png(png_path, width = round(w * dpi), height = round(h * dpi), res = dpi)
      print(p); dev.off()
    }
    message(sprintf("[%s] Saved PT: %s", SCRIPT_NAME, png_path))
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    message(sprintf("[%s] ERROR saving %s: %s", SCRIPT_NAME, name, e$message))
  })
  invisible(list(pdf = pdf_path, png = png_path))
}

# --------------------------------------------------------------------------
# Helper: bootstrap C-index (reused across figures)
# --------------------------------------------------------------------------
bootstrap_cindex <- function(time, event, score_z, n_boot = 1000, seed = 42) {
  set.seed(seed)
  n <- length(time)
  cvals <- numeric(n_boot)
  for (i in seq_len(n_boot)) {
    idx    <- sample(n, n, replace = TRUE)
    cx_raw <- tryCatch(
      concordance(Surv(time[idx], event[idx]) ~ score_z[idx])$concordance,
      error = function(e) NA_real_
    )
    cvals[i] <- max(cx_raw, 1 - cx_raw, na.rm = TRUE)
  }
  c_raw <- concordance(Surv(time, event) ~ score_z)$concordance
  list(
    c_index = max(c_raw, 1 - c_raw),
    ci_low  = quantile(cvals, 0.025, na.rm = TRUE),
    ci_high = quantile(cvals, 0.975, na.rm = TRUE)
  )
}

# ==========================================================================
# FIGURE 1: Fig3_KM_SCANB_OS_CorePAM_PT
# ==========================================================================
message(sprintf("[%s] === 1/15: KM SCAN-B (PT) ===", SCRIPT_NAME))

tryCatch({
  df_scanb <- strict_parquet(file.path(proc_cohort("SCANB"), "analysis_ready.parquet"))
  df_scanb <- df_scanb[!is.na(df_scanb$os_time_months) & df_scanb$os_time_months > 0 &
                          !is.na(df_scanb$os_event) & !is.na(df_scanb$score_z), ]

  # Cox for HR annotation
  cox_uni <- coxph(Surv(os_time_months, os_event) ~ score_z, data = df_scanb)
  sm_uni  <- summary(cox_uni)
  hr_uni  <- sm_uni$conf.int[1, "exp(coef)"]
  lo_uni  <- sm_uni$conf.int[1, "lower .95"]
  hi_uni  <- sm_uni$conf.int[1, "upper .95"]
  p_uni   <- sm_uni$coefficients[1, "Pr(>|z|)"]

  # KM median split
  km_cutpoint <- median(df_scanb$score_z, na.rm = TRUE)
  df_scanb$risk_group_median <- ifelse(df_scanb$score_z >= km_cutpoint, "Alto", "Baixo")

  km_fit <- survfit(Surv(os_time_months, os_event) ~ risk_group_median, data = df_scanb)

  km_plot <- ggsurvplot(
    km_fit,
    data          = df_scanb,
    risk.table    = TRUE,
    pval          = TRUE,
    conf.int      = TRUE,
    palette       = c("#E74C3C", "#2980B9"),
    title         = "KM CorePAM \u2014 SCANB | OS | Mediana intra-coorte (ponto de corte pre-especificado)",
    xlab          = "Tempo (meses)",
    ylab          = "Sobrevida global",
    legend.labs   = c("Alto", "Baixo"),
    risk.table.title = "N\u00famero em risco",
    ggtheme       = theme_classic()
  )

  km_hr_lbl <- sprintf("HR = %.2f (%.2f\u2013%.2f), p = %s",
                        hr_uni, lo_uni, hi_uni,
                        formatC(p_uni, format = "e", digits = 1))
  km_plot$plot <- km_plot$plot +
    ggplot2::annotate("text",
                      x = 5, y = 0.10,
                      label = km_hr_lbl,
                      hjust = 0, size = 3.0, colour = "grey20")

  pdf_path <- file.path(PATHS$figures$main_pt_pdf, "Fig3_KM_SCANB_OS_CorePAM_PT.pdf")
  png_path <- file.path(PATHS$figures$main_pt_png, "Fig3_KM_SCANB_OS_CorePAM_PT.png")

  pdf(pdf_path, width = 8, height = 6); print(km_plot); dev.off()
  png(png_path, width = 2400, height = 1800, res = 300); print(km_plot); dev.off()
  message(sprintf("[%s] Saved: %s", SCRIPT_NAME, png_path))
}, error = function(e) {
  message(sprintf("[%s] ERROR Fig3 KM SCANB PT: %s", SCRIPT_NAME, e$message))
})

# ==========================================================================
# FIGURE 2: Fig5_Calibration_60m_Panels_CorePAM_PT
# ==========================================================================
message(sprintf("[%s] === 2/15: Calibration PT ===", SCRIPT_NAME))

tryCatch({
  cal_files <- list.files(PATHS$results$supp, pattern = "^calibration_\\d+m_.+\\.csv$",
                           full.names = TRUE)
  cal_all <- bind_rows(lapply(cal_files, function(f) {
    tryCatch(read_csv(f, show_col_types = FALSE), error = function(e) NULL)
  }))

  if (!is.null(cal_all) && nrow(cal_all) > 0) {
    if (!"horizon" %in% names(cal_all)) cal_all$horizon <- 60L
    cal_all$facet_label <- sprintf("%s (%dm)", cal_all$cohort, cal_all$horizon)

    p_cal <- ggplot(cal_all, aes(x = pred_mean, y = obs_surv)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
      geom_point(aes(color = cohort), size = 2) +
      geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 0.7,
                  formula = y ~ x) +
      facet_wrap(~facet_label) +
      labs(
        title    = "Calibra\u00e7\u00e3o CorePAM \u2014 Sobrevida Predita vs Observada",
        subtitle = "Horizonte por coorte: TCGA-BRCA=24m; SCANB/METABRIC/GSE20685=60m",
        x        = "Probabilidade de sobrevida predita no horizonte",
        y        = "Sobrevida observada no horizonte (Kaplan-Meier)"
      ) +
      theme_classic(base_size = 11) +
      theme(
        legend.position = "none",
        plot.subtitle   = element_text(size = 9)
      )
    save_pt(p_cal, "Fig5_Calibration_60m_Panels_CorePAM_PT", w = 10, h = 6, section = "main")
  }
}, error = function(e) {
  message(sprintf("[%s] ERROR Calibration PT: %s", SCRIPT_NAME, e$message))
})

# ==========================================================================
# FIGURE 3: Fig5_DeltaCindex_COREA_vs_COREAplus_CorePAM_PT
# ==========================================================================
message(sprintf("[%s] === 3/15: Delta C-index PT ===", SCRIPT_NAME))

tryCatch({
  incr_df <- read_csv(file.path(PATHS$results$main, "incremental_value_by_cohort.csv"),
                       show_col_types = FALSE)
  df_plot <- incr_df[!is.na(incr_df$delta_cindex), ]

  if (nrow(df_plot) > 0) {
    df_plot$corea_label <- dplyr::case_when(
      df_plot$corea_vars_used == "age+er_status" ~ "CORE-A: idade+RE",
      df_plot$corea_vars_used == "age"           ~ "CORE-A: somente idade",
      TRUE                                       ~ paste0("CORE-A: ", df_plot$corea_vars_used)
    )
    df_plot$y_label <- sprintf("%s\n(%s)", df_plot$cohort, df_plot$corea_label)

    p_delta <- ggplot(df_plot, aes(x = delta_cindex, y = reorder(y_label, delta_cindex))) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      geom_errorbar(aes(xmin = delta_ci_lo95, xmax = delta_ci_hi95),
                    orientation = "y", width = 0.3, linewidth = 0.8, color = "#2980B9") +
      geom_point(size = 4, color = "#2980B9") +
      labs(
        title    = "Valor Incremental do CorePAM: Varia\u00e7\u00e3o no C-index",
        subtitle = "CORE-A + CorePAM vs CORE-A isolado | Bootstrap 1.000 itera\u00e7\u00f5es | SCANB=idade+RE; demais=somente idade",
        x        = "Varia\u00e7\u00e3o no C-index (IC 95%)",
        y        = NULL
      ) +
      theme_classic(base_size = 12) +
      theme(
        plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9)
      )
    save_pt(p_delta, "Fig5_DeltaCindex_COREA_vs_COREAplus_CorePAM_PT", w = 8, h = 5, section = "main")
  }
}, error = function(e) {
  message(sprintf("[%s] ERROR Delta Cindex PT: %s", SCRIPT_NAME, e$message))
})

# ==========================================================================
# FIGURE 4: FigS3_Correlation_OffDiagonal_PT
# ==========================================================================
message(sprintf("[%s] === 4/15: Correlation off-diagonal PT ===", SCRIPT_NAME))

tryCatch({
  cor_df <- read_csv(file.path(PATHS$results$supp, "qc_score_correlations_offdiag.csv"),
                      show_col_types = FALSE)
  cohort_names <- setdiff(names(cor_df), "cohort_row")
  cor_mat <- as.matrix(cor_df[, cohort_names])
  rownames(cor_mat) <- cor_df$cohort_row

  cor_long <- do.call(rbind, lapply(cohort_names, function(i) {
    do.call(rbind, lapply(cohort_names, function(j) {
      if (i != j) {
        data.frame(Coorte_A = i, Coorte_B = j,
                   Spearman_rho = cor_mat[i, j],
                   stringsAsFactors = FALSE)
      }
    }))
  }))

  p_cor <- ggplot(cor_long, aes(x = Coorte_A, y = Coorte_B, fill = Spearman_rho)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", Spearman_rho)),
              color = "black", size = 4, na.rm = TRUE) +
    scale_fill_gradient2(
      low  = "#2166AC", mid = "white", high = "#D6604D",
      midpoint = 0, limits = c(-1, 1), na.value = "gray90",
      name = "Spearman rho\n(percentis)"
    ) +
    labs(
      title    = "Correla\u00e7\u00f5es do escore CorePAM (fora da diagonal, via percentis)",
      subtitle = "Diagonal exclu\u00edda (sem autocorrela\u00e7\u00e3o pr\u00f3pria = 1)",
      x        = NULL,
      y        = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1),
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(color = "gray40")
    )

  save_pt(p_cor, "FigS3_Correlation_OffDiagonal_PT", w = 7, h = 6)
}, error = function(e) {
  message(sprintf("[%s] ERROR Correlation PT: %s", SCRIPT_NAME, e$message))
})

# ==========================================================================
# FIGURE 5: FigS4_METABRIC_PCA_forensics_PT
# ==========================================================================
message(sprintf("[%s] === 5/15: METABRIC PCA forensics PT ===", SCRIPT_NAME))

tryCatch({
  COHORT_PCA <- "METABRIC"
  expr_path <- file.path(proc_cohort(COHORT_PCA), "expression_genelevel_preZ.parquet")
  expr_df   <- strict_parquet(expr_path)
  genes     <- expr_df[[1]]
  expr_mat  <- t(as.matrix(expr_df[, -1]))
  colnames(expr_mat) <- genes

  gene_vars <- apply(expr_mat, 2, var, na.rm = TRUE)
  top_genes <- order(gene_vars, decreasing = TRUE)[seq_len(min(500L, ncol(expr_mat)))]
  expr_top  <- expr_mat[, top_genes]
  col_means <- colMeans(expr_top, na.rm = TRUE)
  for (j in seq_len(ncol(expr_top))) {
    na_idx <- is.na(expr_top[, j])
    if (any(na_idx)) expr_top[na_idx, j] <- col_means[j]
  }

  pca_fit  <- prcomp(expr_top, center = TRUE, scale. = TRUE)
  pct_var  <- 100 * pca_fit$sdev^2 / sum(pca_fit$sdev^2)
  pc_scores <- as.data.frame(pca_fit$x[, 1:min(10L, ncol(pca_fit$x))])
  pc_scores$sample_id <- rownames(pca_fit$x)

  # Load clinical
  clin_df   <- strict_parquet(file.path(proc_cohort(COHORT_PCA), "clinical_FINAL.parquet"))
  pc_scores$sample_norm <- normalize_id(pc_scores$sample_id)
  annotated <- merge(pc_scores, clin_df, by.x = "sample_norm", by.y = "sample_id", all.x = TRUE)

  # Outliers
  pc_matrix <- as.matrix(annotated[, c("PC1", "PC2", "PC3")])
  cov_mat   <- cov(pc_matrix, use = "complete.obs")
  center    <- colMeans(pc_matrix, na.rm = TRUE)
  maha_dist <- tryCatch(
    mahalanobis(pc_matrix, center = center, cov = cov_mat),
    error = function(e) rep(NA_real_, nrow(pc_matrix))
  )
  maha_threshold <- qchisq(0.999, df = 3)
  annotated$mahal_dist <- maha_dist
  annotated$is_outlier <- !is.na(maha_dist) & maha_dist > maha_threshold
  n_outliers <- sum(annotated$is_outlier, na.rm = TRUE)

  # ER from raw cBioPortal
  raw_clin_path <- "01_Base_Pura_CorePAM/RAW/METABRIC/brca_metabric/data_clinical_patient.txt"
  er_col <- NULL
  if (file.exists(raw_clin_path)) {
    raw_clin <- tryCatch(
      readr::read_tsv(raw_clin_path, skip = 4, show_col_types = FALSE),
      error = function(e) NULL
    )
    if (!is.null(raw_clin) && "ER_IHC" %in% names(raw_clin)) {
      er_map <- raw_clin[, c("PATIENT_ID", "ER_IHC")]
      names(er_map) <- c("sample_norm", "er_ihc_raw")
      er_map$sample_norm <- toupper(trimws(er_map$sample_norm))
      er_map$er_label_raw <- dplyr::case_when(
        er_map$er_ihc_raw == "Positve"  ~ "RE+",
        er_map$er_ihc_raw == "Negative" ~ "RE\u2212",
        TRUE                            ~ NA_character_
      )
      annotated <- dplyr::left_join(annotated, er_map, by = "sample_norm")
      er_col <- "er_label_raw"
    }
  }

  p1 <- ggplot(annotated, aes(x = PC1, y = PC2)) +
    {if (!is.null(er_col)) aes(color = as.factor(get(er_col))) else NULL} +
    geom_point(alpha = 0.5, size = 1.5) +
    {if (n_outliers > 0 && n_outliers <= 20) {
      geom_label_repel(
        data  = annotated[annotated$is_outlier, ],
        aes(label = sample_id),
        size  = 2.5, color = "red", max.overlaps = 20
      )
    }} +
    labs(
      title   = "PCA Forense METABRIC \u2014 PC1 vs PC2 (500 genes mais vari\u00e1veis)",
      x       = sprintf("PC1 (%.1f%% vari\u00e2ncia)", pct_var[1]),
      y       = sprintf("PC2 (%.1f%% vari\u00e2ncia)", pct_var[2]),
      color   = "Status RE",
      caption = sprintf("Outliers (Mahalanobis p<0,001): n=%d", n_outliers)
    ) +
    theme_classic(base_size = 11)

  p2 <- ggplot(annotated, aes(x = PC1, y = PC3)) +
    {if (!is.null(er_col)) aes(color = as.factor(get(er_col))) else NULL} +
    geom_point(alpha = 0.5, size = 1.5) +
    labs(
      title = "PCA Forense METABRIC \u2014 PC1 vs PC3",
      x     = sprintf("PC1 (%.1f%%)", pct_var[1]),
      y     = sprintf("PC3 (%.1f%%)", pct_var[3]),
      color = "Status RE"
    ) +
    theme_classic(base_size = 11)

  p3 <- ggplot(
    data.frame(pc = seq_len(10), pct = pct_var[seq_len(10)]),
    aes(x = pc, y = pct)
  ) +
    geom_col(fill = "#2980B9") +
    geom_line(color = "red", linewidth = 0.8) +
    geom_point(color = "red", size = 2) +
    labs(
      title = "Scree plot \u2014 Vari\u00e2ncia explicada pela PCA",
      x     = "Componente Principal",
      y     = "Vari\u00e2ncia explicada (%)"
    ) +
    theme_classic(base_size = 11)

  # Composite via gridExtra
  library(gridExtra)

  pdf_path <- file.path(PATHS$figures$supp_pt_pdf, "FigS4_METABRIC_PCA_forensics_PT.pdf")
  png_path <- file.path(PATHS$figures$supp_pt_png, "FigS4_METABRIC_PCA_forensics_PT.png")

  pdf(pdf_path, width = 14, height = 5)
  gridExtra::grid.arrange(p1, p2, p3, ncol = 3)
  dev.off()

  png(png_path, width = 4200, height = 1500, res = 300)
  gridExtra::grid.arrange(p1, p2, p3, ncol = 3)
  dev.off()
  message(sprintf("[%s] Saved: %s", SCRIPT_NAME, png_path))
}, error = function(e) {
  message(sprintf("[%s] ERROR PCA forensics PT: %s", SCRIPT_NAME, e$message))
})

# ==========================================================================
# FIGURE 6: FigS_Bootstrap_GeneFreq_PT
# Regenerate bootstrap gene selection frequency from the model
# ==========================================================================
message(sprintf("[%s] === 6/15: Bootstrap Gene Frequency PT ===", SCRIPT_NAME))

tryCatch({
  # Load training data
  path_expr <- file.path(proc_cohort("SCANB"), "expression_genelevel_preZ.parquet")
  path_clin <- file.path(proc_cohort("SCANB"), "clinical_FINAL.parquet")

  expr_full <- arrow::read_parquet(path_expr)
  PAM50_GENES <- c(
    "ACTR3B","ANLN","BAG1","BCL2","BIRC5","BLVRA","CCNB1","CCNE1","CDC20","CDC6",
    "CDH3","CENPF","CEP55","CXXC5","EGFR","ERBB2","ESR1","EXO1","FGFR4","FOXA1",
    "FOXC1","GPR160","GRB7","KIF2C","KRT14","KRT17","KRT5","MAPT","MDM2",
    "MELK","MIA","MKI67","MLPH","MMP11","MYBL2","MYC","NAT1","NDC80","NUF2",
    "ORC6","PGR","PHGDH","PTTG1","RRM2","SFRP1","SLC39A6","TMEM45B","TYMS","UBE2C","UBE2T"
  )

  pam50_df <- expr_full[expr_full$gene %in% PAM50_GENES, ]
  rm(expr_full)
  gene_names_found <- pam50_df$gene
  X_genes_x_samples <- as.matrix(pam50_df[, !names(pam50_df) %in% "gene"])
  rownames(X_genes_x_samples) <- gene_names_found
  X_raw <- t(X_genes_x_samples)

  clin <- strict_parquet(path_clin)
  common_ids <- intersect(clin$patient_id, rownames(X_raw))
  clin_sub <- clin[match(common_ids, clin$patient_id), ]
  X_sub <- X_raw[common_ids, , drop = FALSE]
  X_scaled <- scale(X_sub)
  y_surv <- Surv(clin_sub$os_time_months, as.integer(clin_sub$os_event))

  # CorePAM genes
  core_genes <- read_csv(file.path(PATHS$results$corepam, "CorePAM_weights.csv"),
                          show_col_types = FALSE)$gene

  ALPHA_EN <- FREEZE$alpha

  # Bootstrap: B=200 refits
  B <- 200
  gene_count <- setNames(integer(ncol(X_scaled)), colnames(X_scaled))

  for (b in seq_len(B)) {
    set.seed(b)
    idx <- sample(nrow(X_scaled), nrow(X_scaled), replace = TRUE)
    fit_b <- tryCatch(
      glmnet::glmnet(
        x = X_scaled[idx, , drop = FALSE],
        y = y_surv[idx, ],
        family = "cox",
        alpha = ALPHA_EN,
        standardize = FALSE
      ),
      error = function(e) NULL
    )
    if (!is.null(fit_b)) {
      # Use lambda with max df
      max_df_idx <- which.max(fit_b$df)
      b_vec <- as.matrix(coef(fit_b, s = fit_b$lambda[max_df_idx]))
      selected <- rownames(b_vec)[b_vec[, 1] != 0]
      gene_count[selected] <- gene_count[selected] + 1L
    }
    if (b %% 50 == 0) message(sprintf("[%s] Bootstrap %d/%d done", SCRIPT_NAME, b, B))
  }

  freq_df <- tibble(
    gene = names(gene_count),
    freq = gene_count / B,
    panel = ifelse(names(gene_count) %in% core_genes, "CorePAM (24)", "Outro PAM50")
  ) |>
    arrange(desc(freq)) |>
    mutate(gene = factor(gene, levels = rev(gene)))

  p_boot <- ggplot(freq_df, aes(x = freq, y = gene, fill = panel)) +
    geom_col(width = 0.75) +
    geom_vline(xintercept = 0.70, linetype = "dashed", color = "red", linewidth = 0.6) +
    scale_fill_manual(
      values = c("CorePAM (24)" = COL$primary, "Outro PAM50" = COL$grey_mid),
      name = "Painel"
    ) +
    labs(
      title    = sprintf("Frequ\u00eancia de Sele\u00e7\u00e3o Bootstrap por Gene (B=%d)", B),
      subtitle = "VALIDA\u00c7\u00c3O T\u00c9CNICA",
      x        = "Frequ\u00eancia de Sele\u00e7\u00e3o",
      y        = "Gene"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold"),
      legend.position = "bottom"
    )

  save_pt(p_boot, "FigS_Bootstrap_GeneFreq_PT", w = 8, h = 10)
}, error = function(e) {
  message(sprintf("[%s] ERROR Bootstrap GeneFreq PT: %s", SCRIPT_NAME, e$message))
})

# ==========================================================================
# FIGURE 7: FigS_COREA_Sensitivity_PT
# Re-generate CORE-A sensitivity analysis (age-only vs age+ER)
# ==========================================================================
message(sprintf("[%s] === 7/15: CORE-A Sensitivity PT ===", SCRIPT_NAME))

tryCatch({
  COHORTS_SENS <- c("SCANB", "TCGA_BRCA", "METABRIC", "GSE20685")
  ENDPOINT_MAP <- list(
    SCANB     = list(time = "os_time_months",  event = "os_event"),
    TCGA_BRCA = list(time = "os_time_months",  event = "os_event"),
    METABRIC  = list(time = "dss_time_months", event = "dss_event"),
    GSE20685  = list(time = "os_time_months",  event = "os_event")
  )

  sens_rows <- list()
  for (coh in COHORTS_SENS) {
    df_coh <- strict_parquet(file.path(proc_cohort(coh), "analysis_ready.parquet"))
    ep     <- ENDPOINT_MAP[[coh]]
    tc     <- ep$time; ec <- ep$event

    if (!all(c(tc, ec, "score_z") %in% names(df_coh))) next
    df_coh <- df_coh[!is.na(df_coh[[tc]]) & df_coh[[tc]] > 0 &
                       !is.na(df_coh[[ec]]) & !is.na(df_coh$score_z), ]

    # Age-only model
    if ("age" %in% names(df_coh) && mean(!is.na(df_coh$age)) >= 0.8) {
      df_age <- df_coh[!is.na(df_coh$age), ]
      fml_age <- as.formula(paste0("Surv(", tc, ",", ec, ") ~ score_z + age"))
      cox_age <- tryCatch(coxph(fml_age, data = df_age), error = function(e) NULL)
      if (!is.null(cox_age)) {
        sm <- summary(cox_age)
        sens_rows[[paste0(coh, "_age")]] <- tibble(
          cohort = coh,
          baseline = "FIXO_somente_idade",
          hr = sm$conf.int["score_z", "exp(coef)"],
          lo = sm$conf.int["score_z", "lower .95"],
          hi = sm$conf.int["score_z", "upper .95"]
        )
      }
    }

    # Age+ER model
    er_avail <- "er_status" %in% names(df_coh) && mean(!is.na(df_coh$er_status)) >= 0.8
    if (er_avail) {
      df_er <- df_coh[!is.na(df_coh$age) & !is.na(df_coh$er_status), ]
      fml_er <- as.formula(paste0("Surv(", tc, ",", ec, ") ~ score_z + age + er_status"))
      cox_er <- tryCatch(coxph(fml_er, data = df_er), error = function(e) NULL)
      if (!is.null(cox_er)) {
        sm <- summary(cox_er)
        sens_rows[[paste0(coh, "_er")]] <- tibble(
          cohort = coh,
          baseline = "ATUAL_idade+RE",
          hr = sm$conf.int["score_z", "exp(coef)"],
          lo = sm$conf.int["score_z", "lower .95"],
          hi = sm$conf.int["score_z", "upper .95"]
        )
      }
    } else {
      sens_rows[[paste0(coh, "_er")]] <- tibble(
        cohort = coh,
        baseline = "ATUAL_somente_idade_(RE_indisponivel)",
        hr = sens_rows[[paste0(coh, "_age")]]$hr,
        lo = sens_rows[[paste0(coh, "_age")]]$lo,
        hi = sens_rows[[paste0(coh, "_age")]]$hi
      )
    }
  }

  sens_df <- bind_rows(sens_rows)
  sens_df$y_label <- sprintf("%s\n%s", sens_df$cohort, sens_df$baseline)

  # Color by baseline type
  sens_df$color_cat <- dplyr::case_when(
    grepl("FIXO", sens_df$baseline)                           ~ "Fixo (somente idade)",
    grepl("ATUAL.*idade\\+RE", sens_df$baseline)              ~ "Atual (idade+RE)",
    grepl("RE_indisponivel|RE_indispon", sens_df$baseline)    ~ "Atual (RE indispon\u00edvel)",
    TRUE                                                      ~ "Outro"
  )

  p_sens <- ggplot(sens_df, aes(x = hr, y = reorder(y_label, hr),
                                 xmin = lo, xmax = hi, color = color_cat)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_pointrange(size = 0.6, linewidth = 0.7) +
    scale_color_manual(values = c(
      "Fixo (somente idade)" = COL$primary,
      "Atual (idade+RE)"     = COL$selected,
      "Atual (RE indispon\u00edvel)" = COL$grey_mid
    ), name = "Linha de base") +
    labs(
      title    = "Sensibilidade da Linha de Base CORE-A: Somente Idade vs Idade+RE",
      subtitle = "VALIDA\u00c7\u00c3O T\u00c9CNICA",
      x        = "Raz\u00e3o de risco (escore CorePAM, por 1 DP)",
      y        = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold"),
      legend.position = "bottom"
    )

  save_pt(p_sens, "FigS_COREA_Sensitivity_PT", w = 10, h = 6)
}, error = function(e) {
  message(sprintf("[%s] ERROR COREA Sensitivity PT: %s", SCRIPT_NAME, e$message))
})

# ==========================================================================
# FIGURE 8: FigS_Cindex_ByCohort_PT
# ==========================================================================
message(sprintf("[%s] === 8/15: C-index by cohort PT ===", SCRIPT_NAME))

tryCatch({
  surv_files <- c(
    file.path(PATHS$results$supp, "survival_results_TCGA_BRCA.csv"),
    file.path(PATHS$results$supp, "survival_results_METABRIC.csv"),
    file.path(PATHS$results$supp, "survival_results_GSE20685.csv")
  )
  tab_surv <- bind_rows(lapply(surv_files, read_csv, show_col_types = FALSE))

  cindex_df <- tab_surv |>
    mutate(
      cohort_label = case_when(
        cohort == "TCGA_BRCA" ~ "TCGA-BRCA\n(OS, 32m SG)",
        cohort == "METABRIC"  ~ "METABRIC\n(SLE, 159m SG)",
        cohort == "GSE20685"  ~ "GSE20685\n(OS, 113m SG)"
      ),
      cohort_label = factor(cohort_label,
                            levels = c("TCGA-BRCA\n(OS, 32m SG)",
                                       "METABRIC\n(SLE, 159m SG)",
                                       "GSE20685\n(OS, 113m SG)"))
    )

  p_cindex <- ggplot(cindex_df,
                     aes(x = cohort_label, y = c_index,
                         ymin = c_index_lo95, ymax = c_index_hi95,
                         fill = cohort_label)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    geom_col(alpha = 0.8, width = 0.55) +
    geom_errorbar(width = 0.18, linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.3f\n(%.3f\u2013%.3f)", c_index, c_index_lo95, c_index_hi95),
                  y = c_index_hi95 + 0.01),
              size = 3.2, fontface = "bold", vjust = 0) +
    scale_fill_manual(values = c(
      "TCGA-BRCA\n(OS, 32m SG)"   = "#2980B9",
      "METABRIC\n(SLE, 159m SG)"  = "#E67E22",
      "GSE20685\n(OS, 113m SG)"   = "#27AE60"
    ), guide = "none") +
    scale_y_continuous(breaks = seq(0.45, 0.75, 0.05)) +
    coord_cartesian(ylim = c(0.45, 0.75)) +
    labs(
      x        = NULL,
      y        = "C-index de Harrell (ajustado, IC 95% bootstrap)",
      title    = "Desempenho Discriminativo do CorePAM por Coorte",
      subtitle = "C_adj = max(C_bruto, 1\u2013C_bruto) | 1.000 reamostras bootstrap | Linha tracejada = 0,5 (aleat\u00f3rio)"
    ) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  save_pt(p_cindex, "FigS_Cindex_ByCohort_PT", w = 9, h = 5)
}, error = function(e) {
  message(sprintf("[%s] ERROR Cindex PT: %s", SCRIPT_NAME, e$message))
})

# ==========================================================================
# FIGURE 9: FigS_CorePAM_vs_PAM50full_HR_PT
# ==========================================================================
message(sprintf("[%s] === 9/15: CorePAM vs PAM50full HR forest PT ===", SCRIPT_NAME))

tryCatch({
  comp_df <- read_csv(file.path(PATHS$results$supp, "pam50full_comparison.csv"),
                       show_col_types = FALSE)

  corepam_df <- read_csv(file.path(PATHS$results$corepam, "CorePAM_weights.csv"),
                          show_col_types = FALSE)
  pam50full_df <- read_csv(file.path(PATHS$results$corepam, "PAM50full_weights.csv"),
                            show_col_types = FALSE)

  # Build pooled FE estimate
  pool_row <- function(loghr_col, se_col) {
    rows <- comp_df |> select(cohort, all_of(c(loghr_col, se_col))) |>
      filter(!is.na(.data[[se_col]]), .data[[se_col]] > 0)
    w  <- 1 / rows[[se_col]]^2
    lw <- sum(w * rows[[loghr_col]]) / sum(w)
    se <- 1 / sqrt(sum(w))
    tibble(cohort = "Combinado (EF)", loghr = lw, se_loghr = se,
           hr = exp(lw), hr_lo = exp(lw - 1.96 * se), hr_hi = exp(lw + 1.96 * se))
  }

  forest_df_cp <- comp_df |>
    transmute(cohort, model = "CorePAM",
              hr = hr_corepam, hr_lo = hr_corepam_lo95, hr_hi = hr_corepam_hi95,
              n = n_samples, endpoint)
  forest_df_pf <- comp_df |>
    transmute(cohort, model = "PAM50-completo",
              hr = hr_pam50full, hr_lo = hr_pam50full_lo95, hr_hi = hr_pam50full_hi95,
              n = n_samples, endpoint)

  pool_cp <- pool_row("loghr_corepam", "se_loghr_corepam") |>
    mutate(model = "CorePAM", n = sum(comp_df$n_samples), endpoint = "combinado")
  pool_pf <- pool_row("loghr_pam50full", "se_loghr_pam50full") |>
    mutate(model = "PAM50-completo", n = sum(comp_df$n_samples), endpoint = "combinado")

  forest_df <- bind_rows(forest_df_cp, forest_df_pf, pool_cp, pool_pf)

  cohort_order <- c("SCAN-B", "TCGA-BRCA", "METABRIC", "GSE20685", "Combinado (EF)")
  forest_df$cohort <- factor(forest_df$cohort, levels = rev(cohort_order))
  forest_df$model  <- factor(forest_df$model, levels = c("CorePAM", "PAM50-completo"))

  forest_df <- forest_df |>
    mutate(
      cohort_label = case_when(
        cohort == "Combinado (EF)" ~ "Combinado (EF)",
        TRUE ~ sprintf("%s\n(n=%d, %s)", cohort, n, endpoint)
      )
    )
  forest_df$cohort_label <- factor(forest_df$cohort_label,
                                    levels = rev(c(
                                      sprintf("SCAN-B\n(n=%d, OS)",     comp_df$n_samples[comp_df$cohort == "SCAN-B"]),
                                      sprintf("TCGA-BRCA\n(n=%d, OS)",  comp_df$n_samples[comp_df$cohort == "TCGA-BRCA"]),
                                      sprintf("METABRIC\n(n=%d, DSS)",  comp_df$n_samples[comp_df$cohort == "METABRIC"]),
                                      sprintf("GSE20685\n(n=%d, OS)",   comp_df$n_samples[comp_df$cohort == "GSE20685"]),
                                      "Combinado (EF)"
                                    )))

  p_forest_hr <- ggplot(
    forest_df,
    aes(x = hr, y = cohort_label, color = model, shape = model)
  ) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.7) +
    geom_pointrange(
      aes(xmin = hr_lo, xmax = hr_hi),
      position = position_dodge(width = 0.5),
      linewidth = 0.7, size = 0.5
    ) +
    scale_color_manual(
      name   = "Escore",
      values = c("CorePAM" = "#2166AC", "PAM50-completo" = "#888888")
    ) +
    scale_shape_manual(
      name   = "Escore",
      values = c("CorePAM" = 16, "PAM50-completo" = 17)
    ) +
    scale_x_log10(
      breaks = c(0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0),
      labels = c("0,5","0,8","1,0","1,2","1,5","2,0","2,5","3,0")
    ) +
    labs(
      title    = "HR por 1 DP: CorePAM vs PAM50-completo",
      subtitle = sprintf("CorePAM: %d genes | PAM50-completo: %d genes | Cox univariado",
                         nrow(corepam_df), nrow(pam50full_df)),
      x        = "Raz\u00e3o de risco por 1 DP (IC 95%, escala log)",
      y        = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position    = "bottom",
      plot.title         = element_text(face = "bold"),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(linetype = "dotted", color = "grey80"),
      axis.text.y        = element_text(size = 9)
    )

  save_pt(p_forest_hr, "FigS_CorePAM_vs_PAM50full_HR_PT", w = 9, h = 6)
}, error = function(e) {
  message(sprintf("[%s] ERROR CorePAM vs PAM50full HR PT: %s", SCRIPT_NAME, e$message))
})

# ==========================================================================
# FIGURE 10: FigS_DropGene_Tornado_PT
# Leave-one-gene-out C-index impact
# ==========================================================================
message(sprintf("[%s] === 10/15: Drop Gene Tornado PT ===", SCRIPT_NAME))

tryCatch({
  core_weights <- read_csv(file.path(PATHS$results$corepam, "CorePAM_weights.csv"),
                            show_col_types = FALSE)
  w_named <- setNames(core_weights$weight, core_weights$gene)

  COHORTS_DROP <- c("SCANB", "TCGA_BRCA", "METABRIC", "GSE20685")
  ENDPOINT_MAP_DROP <- list(
    SCANB     = list(time = "os_time_months",  event = "os_event"),
    TCGA_BRCA = list(time = "os_time_months",  event = "os_event"),
    METABRIC  = list(time = "dss_time_months", event = "dss_event"),
    GSE20685  = list(time = "os_time_months",  event = "os_event")
  )

  compute_score_from_ready <- function(df, weights, genes_to_use) {
    w_sub <- weights[genes_to_use]
    # score_z is pre-computed; we need to reconstruct from expression
    # Instead, compute delta-C by comparing full vs leave-one-out scores
    # We'll use the score formula: sum(w_i * z_i) / sum(|w_i|)
    # Since we only have score_z, approximate using the analysis_ready's score_z
    # Actually, we need to use the raw approach
    NULL
  }

  drop_rows <- list()
  for (coh in COHORTS_DROP) {
    ready <- strict_parquet(file.path(proc_cohort(coh), "analysis_ready.parquet"))
    ep    <- ENDPOINT_MAP_DROP[[coh]]
    tc    <- ep$time; ec <- ep$event

    # Load expression for score recomputation
    expr_path <- file.path(proc_cohort(coh), "expression_genelevel_preZ.parquet")
    expr_df   <- arrow::read_parquet(expr_path)
    gene_names_e <- expr_df$gene
    sample_ids   <- setdiff(names(expr_df), "gene")
    expr_vals    <- t(as.matrix(expr_df[, sample_ids]))
    colnames(expr_vals) <- gene_names_e

    # Z-score intra-cohort
    gene_means <- colMeans(expr_vals, na.rm = TRUE)
    gene_sds   <- apply(expr_vals, 2, sd, na.rm = TRUE)
    zero_sd <- gene_sds == 0 | is.na(gene_sds)
    gene_sds[zero_sd] <- 1
    z_mat <- sweep(sweep(expr_vals, 2, gene_means, "-"), 2, gene_sds, "/")
    z_mat[, zero_sd] <- 0

    # Match to clinical
    id_col <- "patient_id"
    if (length(intersect(rownames(z_mat), ready[[id_col]])) == 0 &&
        "sample_id" %in% names(ready)) {
      id_col <- "sample_id"
    }
    common <- intersect(rownames(z_mat), ready[[id_col]])
    if (length(common) < 50) next

    z_sub  <- z_mat[common, , drop = FALSE]
    df_coh <- ready[match(common, ready[[id_col]]), ]
    df_coh <- df_coh[!is.na(df_coh[[tc]]) & df_coh[[tc]] > 0 &
                       !is.na(df_coh[[ec]]), ]
    common2 <- intersect(common, df_coh[[id_col]])
    z_sub  <- z_sub[common2, , drop = FALSE]
    df_coh <- df_coh[match(common2, df_coh[[id_col]]), ]

    time_v  <- df_coh[[tc]]
    event_v <- df_coh[[ec]]

    # Determine score direction from ready data
    score_dir <- if ("score_direction" %in% names(df_coh)) {
      unique(df_coh$score_direction)[1]
    } else { "original" }
    dir_sign <- if (score_dir == "inverted") -1 else 1

    genes_present <- intersect(names(w_named), colnames(z_sub))

    # Full score
    w_full <- w_named[genes_present]
    score_full <- as.vector(z_sub[, genes_present, drop = FALSE] %*% w_full) / sum(abs(w_full))
    score_full_z <- as.vector(scale(dir_sign * score_full))
    c_full <- tryCatch({
      c_raw <- concordance(Surv(time_v, event_v) ~ score_full_z)$concordance
      max(c_raw, 1 - c_raw)
    }, error = function(e) NA_real_)

    # Leave-one-gene-out
    for (g in genes_present) {
      genes_loo <- setdiff(genes_present, g)
      w_loo     <- w_named[genes_loo]
      score_loo <- as.vector(z_sub[, genes_loo, drop = FALSE] %*% w_loo) / sum(abs(w_loo))
      score_loo_z <- as.vector(scale(dir_sign * score_loo))
      c_loo <- tryCatch({
        c_raw <- concordance(Surv(time_v, event_v) ~ score_loo_z)$concordance
        max(c_raw, 1 - c_raw)
      }, error = function(e) NA_real_)

      drop_rows[[paste0(coh, "_", g)]] <- tibble(
        cohort    = coh,
        gene      = g,
        c_full    = c_full,
        c_loo     = c_loo,
        delta_c   = c_loo - c_full
      )
    }
    message(sprintf("[%s] Drop-gene %s done (%d genes)", SCRIPT_NAME, coh, length(genes_present)))
  }

  drop_df <- bind_rows(drop_rows)

  # Order genes by median delta across cohorts
  gene_order <- drop_df |>
    group_by(gene) |>
    summarise(med_delta = median(delta_c, na.rm = TRUE), .groups = "drop") |>
    arrange(med_delta) |>
    pull(gene)

  drop_df$gene <- factor(drop_df$gene, levels = gene_order)
  drop_df$direction <- ifelse(drop_df$delta_c >= 0, "Melhora (gene dispens\u00e1vel)",
                               "Piora (gene importante)")

  p_tornado <- ggplot(drop_df, aes(x = delta_c, y = gene, fill = direction)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
    facet_wrap(~ cohort, scales = "free_x") +
    scale_fill_manual(values = c(
      "Melhora (gene dispens\u00e1vel)" = COL$primary,
      "Piora (gene importante)" = COL$selected
    ), name = "Dire\u00e7\u00e3o") +
    labs(
      title    = "Remo\u00e7\u00e3o de Um Gene por Vez: Impacto no C-index",
      subtitle = "VALIDA\u00c7\u00c3O T\u00c9CNICA",
      x        = "\u0394C-index",
      y        = "Gene Removido"
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title    = element_text(face = "bold"),
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey90"),
      strip.text       = element_text(face = "bold")
    )

  save_pt(p_tornado, "FigS_DropGene_Tornado_PT", w = 10, h = 10)
}, error = function(e) {
  message(sprintf("[%s] ERROR DropGene Tornado PT: %s", SCRIPT_NAME, e$message))
})

# ==========================================================================
# FIGURE 11: FigS_Forest_HR_ValidationCohorts_PT
# ==========================================================================
message(sprintf("[%s] === 11/15: Forest HR Validation Cohorts PT ===", SCRIPT_NAME))

tryCatch({
  surv_files <- c(
    file.path(PATHS$results$supp, "survival_results_TCGA_BRCA.csv"),
    file.path(PATHS$results$supp, "survival_results_METABRIC.csv"),
    file.path(PATHS$results$supp, "survival_results_GSE20685.csv")
  )
  tab_surv <- bind_rows(lapply(surv_files, read_csv, show_col_types = FALSE)) |>
    select(cohort, endpoint, n_samples, n_events, hr_uni, hr_uni_lo95, hr_uni_hi95,
           p_uni, loghr_uni, se_loghr_uni, c_index, c_index_lo95, c_index_hi95)

  # Random-effects pooling (DerSimonian-Laird)
  k    <- nrow(tab_surv)
  yi   <- tab_surv$loghr_uni
  vi   <- tab_surv$se_loghr_uni^2
  wi_f <- 1 / vi
  q_stat   <- sum(wi_f * (yi - sum(wi_f * yi) / sum(wi_f))^2)
  df_q     <- k - 1
  c_const  <- sum(wi_f) - sum(wi_f^2) / sum(wi_f)
  tau2_dl  <- max(0, (q_stat - df_q) / c_const)
  wi_re    <- 1 / (vi + tau2_dl)
  pooled_loghr <- sum(wi_re * yi) / sum(wi_re)
  pooled_se    <- sqrt(1 / sum(wi_re))
  pooled_hr    <- exp(pooled_loghr)
  pooled_lo    <- exp(pooled_loghr - 1.96 * pooled_se)
  pooled_hi    <- exp(pooled_loghr + 1.96 * pooled_se)
  pooled_p     <- 2 * pnorm(-abs(pooled_loghr / pooled_se))
  i2_pct   <- max(0, 100 * (q_stat - df_q) / q_stat)
  p_het    <- pchisq(q_stat, df = df_q, lower.tail = FALSE)

  # Fixed-effects
  pooled_loghr_fe <- sum(wi_f * yi) / sum(wi_f)
  pooled_se_fe    <- sqrt(1 / sum(wi_f))
  pooled_hr_fe    <- exp(pooled_loghr_fe)
  pooled_lo_fe    <- exp(pooled_loghr_fe - 1.96 * pooled_se_fe)
  pooled_hi_fe    <- exp(pooled_loghr_fe + 1.96 * pooled_se_fe)
  pooled_p_fe     <- 2 * pnorm(-abs(pooled_loghr_fe / pooled_se_fe))

  forest_df <- tab_surv |>
    mutate(
      label = sprintf("%s (%s)\nN=%d, Eventos=%d", cohort, endpoint, n_samples, n_events),
      p_label = ifelse(p_uni < 0.001,
                       formatC(p_uni, format = "e", digits = 1),
                       sprintf("%.4f", p_uni)),
      c_label = sprintf("%.3f (%.3f\u2013%.3f)", c_index, c_index_lo95, c_index_hi95),
      right_label = sprintf("HR=%.2f (%.2f\u2013%.2f) | p=%s | C=%.3f",
                            hr_uni, hr_uni_lo95, hr_uni_hi95, p_label, c_index),
      y_pos = rev(seq_len(n()))
    )

  pooled_row_re <- tibble(
    cohort = "Combinado (EA/DL)", endpoint = "",
    n_samples = sum(tab_surv$n_samples), n_events = sum(tab_surv$n_events),
    hr_uni = pooled_hr, hr_uni_lo95 = pooled_lo, hr_uni_hi95 = pooled_hi,
    p_uni = pooled_p, loghr_uni = pooled_loghr, se_loghr_uni = pooled_se,
    c_index = NA_real_, c_index_lo95 = NA_real_, c_index_hi95 = NA_real_,
    label = sprintf("Combinado EA (DerSimonian-Laird)\nI\u00b2=%.0f%%, \u03c4\u00b2=%.4f, p_het=%.3f",
                    i2_pct, tau2_dl, p_het),
    p_label = formatC(pooled_p, format = "e", digits = 1),
    c_label = NA_character_,
    right_label = sprintf("HR=%.2f (%.2f\u2013%.2f) | p=%s | EA",
                          pooled_hr, pooled_lo, pooled_hi,
                          formatC(pooled_p, format = "e", digits = 1)),
    y_pos = 0
  )

  pooled_row_fe <- tibble(
    cohort = "Combinado (EF/IVP)", endpoint = "",
    n_samples = sum(tab_surv$n_samples), n_events = sum(tab_surv$n_events),
    hr_uni = pooled_hr_fe, hr_uni_lo95 = pooled_lo_fe, hr_uni_hi95 = pooled_hi_fe,
    p_uni = pooled_p_fe, loghr_uni = pooled_loghr_fe, se_loghr_uni = pooled_se_fe,
    c_index = NA_real_, c_index_lo95 = NA_real_, c_index_hi95 = NA_real_,
    label = "Combinado EF (vari\u00e2ncia inversa)\n(assume homogeneidade)",
    p_label = formatC(pooled_p_fe, format = "e", digits = 1),
    c_label = NA_character_,
    right_label = sprintf("HR=%.2f (%.2f\u2013%.2f) | p=%s | EF",
                          pooled_hr_fe, pooled_lo_fe, pooled_hi_fe,
                          formatC(pooled_p_fe, format = "e", digits = 1)),
    y_pos = -1
  )

  forest_all <- bind_rows(forest_df, pooled_row_re, pooled_row_fe) |>
    mutate(
      is_pooled_re = (cohort == "Combinado (EA/DL)"),
      is_pooled_fe = (cohort == "Combinado (EF/IVP)"),
      pt_shape = dplyr::case_when(is_pooled_re ~ 18, is_pooled_fe ~ 23, TRUE ~ 15),
      pt_size  = dplyr::case_when(is_pooled_re ~ 5, is_pooled_fe ~ 4, TRUE ~ 3.5),
      color_cat = dplyr::case_when(
        cohort == "TCGA_BRCA"          ~ "TCGA-BRCA",
        cohort == "METABRIC"           ~ "METABRIC",
        cohort == "GSE20685"           ~ "GSE20685",
        cohort == "Combinado (EA/DL)"  ~ "Combinado EA",
        cohort == "Combinado (EF/IVP)" ~ "Combinado EF"
      )
    )

  x_right <- max(forest_all$hr_uni_hi95, na.rm = TRUE) * 1.05
  x_xlim  <- max(forest_all$hr_uni_hi95, na.rm = TRUE) * 2.5

  p_forest <- ggplot(forest_all,
                     aes(x = hr_uni, xmin = hr_uni_lo95, xmax = hr_uni_hi95,
                         y = y_pos)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    geom_errorbar(aes(color = color_cat),
                  width = 0.25, linewidth = 0.7, na.rm = TRUE, orientation = "y") +
    geom_point(aes(shape = pt_shape, size = pt_size, color = color_cat), na.rm = TRUE) +
    geom_text(aes(label = right_label),
              x = x_right, hjust = 0, size = 2.8, family = "sans") +
    scale_shape_identity() +
    scale_size_identity() +
    scale_color_manual(values = c(
      "TCGA-BRCA"   = "#2980B9",
      "METABRIC"    = "#E67E22",
      "GSE20685"    = "#27AE60",
      "Combinado EA" = "#2C3E50",
      "Combinado EF" = "#7F8C8D"
    )) +
    scale_x_continuous(
      name   = "Raz\u00e3o de Risco (por 1 DP do escore CorePAM) com IC 95%",
      limits = c(0.85, x_xlim),
      breaks = c(0.9, 1.0, 1.2, 1.4, 1.6, 1.8)
    ) +
    scale_y_continuous(
      breaks = forest_all$y_pos,
      labels = forest_all$label
    ) +
    annotate("segment", x = -Inf, xend = -Inf,
             y = 0.4, yend = nrow(forest_df) + 0.4, color = "grey40") +
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = -1.5, ymax = 0.5, fill = "#ECF0F1", alpha = 0.4) +
    labs(
      title    = "Escore CorePAM \u2014 Raz\u00f5es de Risco nas Coortes de Valida\u00e7\u00e3o",
      subtitle = "HR por 1 DP de score_z (Z-score intra-coorte) | TCGA-BRCA/GSE20685=OS, METABRIC=SLE | Cox univariado",
      caption  = sprintf("Heterogeneidade: I\u00b2=%.0f%%, \u03c4\u00b2=%.4f, Q=%.2f (df=%d, p_het=%.3f) | EA=efeitos aleat\u00f3rios (DL); EF=efeitos fixos (IVP; assume homogeneidade)",
                         i2_pct, tau2_dl, q_stat, df_q, p_het),
      color    = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.title.y    = element_blank(),
      axis.line.y     = element_blank(),
      axis.ticks.y    = element_blank(),
      legend.position = "none",
      plot.title      = element_text(face = "bold", size = 12),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3)
    )

  save_pt(p_forest, "FigS_Forest_HR_ValidationCohorts_PT", w = 12, h = 6)
}, error = function(e) {
  message(sprintf("[%s] ERROR Forest HR PT: %s", SCRIPT_NAME, e$message))
})

# ==========================================================================
# FIGURES 12-14: Heatmaps (TCGA_BRCA, METABRIC, GSE20685) PT
# ==========================================================================
message(sprintf("[%s] === 12-14/15: Heatmaps PT ===", SCRIPT_NAME))

core_genes <- c("ACTR3B","BCL2","BLVRA","CENPF","CXXC5","ERBB2","ESR1","EXO1",
                "FGFR4","FOXC1","GPR160","GRB7","KRT17","KRT5","MDM2","MIA",
                "MLPH","MYBL2","MYC","NAT1","PGR","PHGDH","PTTG1","SFRP1")
weights_df <- read_csv(file.path(PATHS$results$corepam, "CorePAM_weights.csv"),
                        show_col_types = FALSE)

heatmap_configs <- list(
  list(name = "TCGA_BRCA",
       title = "Mapa de Calor de Express\u00e3o G\u00eanica CorePAM \u2014 TCGA-BRCA (OS, n=200 aleat\u00f3rios, semente=42)",
       fig_name = "FigS_Heatmap_TCGA_BRCA_PT"),
  list(name = "METABRIC",
       title = "Mapa de Calor de Express\u00e3o G\u00eanica CorePAM \u2014 METABRIC (SLE, n=200 aleat\u00f3rios, semente=42)",
       fig_name = "FigS_Heatmap_METABRIC_PT"),
  list(name = "GSE20685",
       title = "Mapa de Calor de Express\u00e3o G\u00eanica CorePAM \u2014 GSE20685 (OS, todas as amostras)",
       fig_name = "FigS_Heatmap_GSE20685_PT")
)

for (cfg in heatmap_configs) {
  tryCatch({
    message(sprintf("[%s] Heatmap %s PT", SCRIPT_NAME, cfg$name))

    expr_p <- file.path(PATHS$processed, cfg$name, "expression_genelevel_preZ.parquet")
    rdy_p  <- file.path(PATHS$processed, cfg$name, "analysis_ready.parquet")

    expr_mat <- arrow::read_parquet(expr_p)
    rdy_df   <- arrow::read_parquet(rdy_p)

    genes_in <- intersect(core_genes, expr_mat$gene)
    expr_sub <- expr_mat[expr_mat$gene %in% genes_in, ]
    gene_names   <- expr_sub$gene
    sample_ids   <- setdiff(names(expr_sub), "gene")
    expr_vals    <- t(as.matrix(expr_sub[, sample_ids]))
    colnames(expr_vals) <- gene_names

    z_mat <- scale(expr_vals)

    id_col <- "patient_id"
    if (length(intersect(rownames(z_mat), rdy_df[[id_col]])) == 0 &&
        "sample_id" %in% names(rdy_df)) {
      id_col <- "sample_id"
    }
    common_samp <- intersect(rownames(z_mat), rdy_df[[id_col]])
    z_sub   <- z_mat[common_samp, , drop = FALSE]
    rdy_sub <- rdy_df[match(common_samp, rdy_df[[id_col]]), ]

    set.seed(42)
    n_show   <- min(200, nrow(z_sub))
    idx      <- sample(nrow(z_sub), n_show)
    z_show   <- z_sub[idx, ]
    score_z_show <- rdy_sub$score_z[idx]
    ann_df <- data.frame(
      grupo_risco = ifelse(score_z_show >= 0, "Alto", "Baixo"),
      score_z     = score_z_show,
      row.names   = rownames(z_show)
    )

    # Sort columns by weight
    g_order <- weights_df$gene[weights_df$gene %in% colnames(z_show)]
    g_order <- g_order[order(weights_df$weight[match(g_order, weights_df$gene)])]
    z_show  <- z_show[, g_order, drop = FALSE]

    # Sort rows by group then score
    row_ord <- order(ann_df$grupo_risco, ann_df$score_z)
    z_show  <- z_show[row_ord, ]
    ann_plot <- ann_df[row_ord, "grupo_risco", drop = FALSE]

    if (requireNamespace("pheatmap", quietly = TRUE)) {
      ann_colors <- list(grupo_risco = c(Alto = COL$km_high, Baixo = COL$km_low))
      heat_args  <- list(
        t(z_show),
        annotation_col    = ann_plot,
        annotation_colors = ann_colors,
        cluster_rows      = FALSE,
        cluster_cols      = FALSE,
        show_colnames     = FALSE,
        color             = colorRampPalette(c("#2C3E50","white","#C0392B"))(100),
        breaks            = seq(-3, 3, length.out = 101),
        main              = cfg$title,
        fontsize_row      = 8,
        border_color      = NA
      )

      png_path <- file.path(PATHS$figures$supp_pt_png, paste0(cfg$fig_name, ".png"))
      tryCatch({
        do.call(pheatmap::pheatmap,
                c(heat_args, list(filename = png_path, width = 10, height = 6)))
        message(sprintf("[%s] Heatmap PNG saved: %s", SCRIPT_NAME, png_path))
      }, error = function(e) {
        message(sprintf("[%s] Heatmap PNG error (%s): %s", SCRIPT_NAME, cfg$fig_name, e$message))
      })

      pdf_path <- file.path(PATHS$figures$supp_pt_pdf, paste0(cfg$fig_name, ".pdf"))
      tryCatch({
        pdf(pdf_path, width = 10, height = 6)
        do.call(pheatmap::pheatmap, heat_args)
        dev.off()
        message(sprintf("[%s] Heatmap PDF saved: %s", SCRIPT_NAME, pdf_path))
      }, error = function(e) {
        try(dev.off(), silent = TRUE)
        message(sprintf("[%s] Heatmap PDF error (%s): %s", SCRIPT_NAME, cfg$fig_name, e$message))
      })
    } else {
      message(sprintf("[%s] pheatmap not available; skipping heatmap %s", SCRIPT_NAME, cfg$name))
    }
  }, error = function(e) {
    message(sprintf("[%s] ERROR Heatmap %s PT: %s", SCRIPT_NAME, cfg$name, e$message))
  })
}

# ==========================================================================
# FIGURE 15: FigS_Weights_CorePAM_Lollipop_PT
# ==========================================================================
message(sprintf("[%s] === 15/15: Lollipop weights PT ===", SCRIPT_NAME))

tryCatch({
  wt_df <- weights_df |>
    mutate(
      direction = ifelse(weight > 0, "Risco maior (+)", "Risco menor (\u2212)"),
      gene      = reorder(gene, weight)
    )

  pal_dir <- c("Risco maior (+)" = "#C0392B", "Risco menor (\u2212)" = "#2980B9")

  p_lollipop <- ggplot(wt_df, aes(x = gene, y = weight, color = direction)) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.5) +
    geom_segment(aes(x = gene, xend = gene, y = 0, yend = weight),
                 linewidth = 0.8) +
    geom_point(aes(size = abs(weight)), shape = 16) +
    ggrepel::geom_text_repel(
      aes(label = gene),
      size = 2.8, color = "black",
      direction = "y", nudge_x = 0.3,
      segment.size = 0.2, segment.color = "grey60",
      box.padding = 0.2, point.padding = 0.3,
      min.segment.length = 0.1,
      max.overlaps = 30
    ) +
    scale_color_manual(values = pal_dir, name = "Dire\u00e7\u00e3o do efeito") +
    scale_size_continuous(range = c(2, 6), guide = "none") +
    coord_flip() +
    labs(
      x        = NULL,
      y        = "Coeficiente elastic-net (peso)",
      title    = "Painel CorePAM \u2014 Pesos dos Genes",
      subtitle = sprintf("%d genes | N\u00e3o-inferioridade OOF (\u0394C=0,010) | Treinamento SCAN-B | Pesos do reajuste no conjunto completo de treinamento",
                         nrow(wt_df))
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.y     = element_blank(),
      axis.ticks.y    = element_blank(),
      legend.position = "bottom",
      plot.title      = element_text(face = "bold")
    )

  save_pt(p_lollipop, "FigS_Weights_CorePAM_Lollipop_PT", w = 8, h = 6)
}, error = function(e) {
  message(sprintf("[%s] ERROR Lollipop PT: %s", SCRIPT_NAME, e$message))
})

# ==========================================================================
# DONE
# ==========================================================================
message(sprintf("[%s] === ALL 15 PT FIGURES COMPLETED ===", SCRIPT_NAME))
message(sprintf("[%s] Main PT figures: %s", SCRIPT_NAME, PATHS$figures$main_pt_png))
message(sprintf("[%s] Supp PT figures: %s", SCRIPT_NAME, PATHS$figures$supp_pt_png))
