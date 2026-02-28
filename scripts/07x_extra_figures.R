# =============================================================================
# SCRIPT: 07x_extra_figures.R
# PURPOSE: Supplementary and enrichment figures beyond KM curves:
#   FigS_Forest       — Forest plot HR across validation cohorts
#   FigS_ScoreDist    — Score_z distribution per cohort (density overlay)
#   FigS_Heatmap      — Heatmap of Core-PAM gene expression z-scores (validation)
#   FigS_Weights      — Lollipop chart of Core-PAM gene weights
#   FigS_ScoreBySubtype — Score by PAM50/ER subtype (if available)
# PROJETO: Core-PAM (Memorial v6.1)
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "07x_extra_figures.R"

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(arrow)
  library(patchwork)
})
source("Y:/Phd-Genomic-claude/scripts/00_colors.R")

fig_main <- PATHS$figures$main
fig_supp <- PATHS$figures$supp
dir.create(fig_main, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_supp, showWarnings = FALSE, recursive = TRUE)

# --------------------------------------------------------------------------
# Helper: save PDF + PNG
# --------------------------------------------------------------------------
save_fig <- function(p, name, w = 8, h = 6, dir = fig_supp, dpi = 300) {
  pdf_path <- file.path(dir, paste0(name, ".pdf"))
  png_path <- file.path(dir, paste0(name, ".png"))
  old_warn <- getOption("warn"); options(warn = 0)
  tryCatch({
    cairo_pdf(pdf_path, width = w, height = h); print(p); dev.off()
    png(png_path, width = w, height = h, units = "in", res = dpi); print(p); dev.off()
    registry_append("ALL", name, pdf_path, sha256_file(pdf_path), "ok",
                    SCRIPT_NAME, file.info(pdf_path)$size / 1e6)
    message(sprintf("[%s] Saved: %s | %s", SCRIPT_NAME, pdf_path, png_path))
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    message(sprintf("[%s] ERROR saving %s: %s", SCRIPT_NAME, name, e$message))
  })
  options(warn = old_warn)
  invisible(list(pdf = pdf_path, png = png_path))
}

# --------------------------------------------------------------------------
# 1) FOREST PLOT — HR across validation cohorts
# --------------------------------------------------------------------------
message(sprintf("[%s] === 1) Forest plot ===", SCRIPT_NAME))

surv_files <- c(
  file.path(PATHS$results$supp, "survival_results_TCGA_BRCA.csv"),
  file.path(PATHS$results$supp, "survival_results_METABRIC.csv"),
  file.path(PATHS$results$supp, "survival_results_GSE20685.csv")
)

tab_surv <- bind_rows(lapply(surv_files, read_csv, show_col_types = FALSE)) |>
  select(cohort, endpoint, n_samples, n_events, hr_uni, hr_uni_lo95, hr_uni_hi95,
         p_uni, loghr_uni, se_loghr_uni, c_index, c_index_lo95, c_index_hi95)

# Pooled (inverse-variance random-effects, simplified: weighted mean)
# Use logHR and SE for proper pooling
wi      <- 1 / tab_surv$se_loghr_uni^2
pooled_loghr <- sum(wi * tab_surv$loghr_uni) / sum(wi)
pooled_se    <- sqrt(1 / sum(wi))
pooled_hr    <- exp(pooled_loghr)
pooled_lo    <- exp(pooled_loghr - 1.96 * pooled_se)
pooled_hi    <- exp(pooled_loghr + 1.96 * pooled_se)
pooled_p     <- 2 * pnorm(-abs(pooled_loghr / pooled_se))

# Build data for forest plot
forest_df <- tab_surv |>
  mutate(
    label       = sprintf("%s (%s)\nN=%d, Events=%d", cohort, endpoint, n_samples, n_events),
    p_label     = ifelse(p_uni < 0.001,
                         formatC(p_uni, format = "e", digits = 1),
                         sprintf("%.4f", p_uni)),
    c_label     = sprintf("%.3f (%.3f–%.3f)", c_index, c_index_lo95, c_index_hi95),
    right_label = sprintf("HR=%.2f (%.2f–%.2f) | p=%s | C=%.3f",
                          hr_uni, hr_uni_lo95, hr_uni_hi95, p_label, c_index),
    y_pos       = rev(seq_len(n()))
  )

# Add pooled row
pooled_row <- tibble(
  cohort    = "Pooled (IVW)",
  endpoint  = "",
  n_samples = sum(tab_surv$n_samples),
  n_events  = sum(tab_surv$n_events),
  hr_uni    = pooled_hr,
  hr_uni_lo95 = pooled_lo,
  hr_uni_hi95 = pooled_hi,
  p_uni     = pooled_p,
  loghr_uni = pooled_loghr,
  se_loghr_uni = pooled_se,
  c_index   = NA_real_, c_index_lo95 = NA_real_, c_index_hi95 = NA_real_,
  label     = sprintf("Pooled IVW (fixed)\nN=%d, Events=%d",
                      sum(tab_surv$n_samples), sum(tab_surv$n_events)),
  p_label   = formatC(pooled_p, format = "e", digits = 1),
  c_label   = NA_character_,
  right_label = sprintf("HR=%.2f (%.2f–%.2f) | p=%s",
                        pooled_hr, pooled_lo, pooled_hi,
                        formatC(pooled_p, format = "e", digits = 1)),
  y_pos     = 0
)

forest_all <- bind_rows(forest_df, pooled_row) |>
  mutate(
    is_pooled = (cohort == "Pooled (IVW)"),
    pt_shape  = ifelse(is_pooled, 18, 15),  # diamond for pooled
    pt_size   = ifelse(is_pooled, 5, 3.5),
    color_cat = case_when(
      cohort == "TCGA_BRCA"       ~ "TCGA-BRCA",
      cohort == "METABRIC"        ~ "METABRIC",
      cohort == "GSE20685"        ~ "GSE20685",
      cohort == "Pooled (IVW)"   ~ "Pooled"
    )
  )

# Y-axis labels
y_labels <- setNames(forest_all$label, forest_all$y_pos)
x_right  <- max(forest_all$hr_uni_hi95, na.rm = TRUE) * 1.05
x_xlim   <- max(forest_all$hr_uni_hi95, na.rm = TRUE) * 2.5

p_forest <- ggplot(forest_all,
                   aes(x = hr_uni, xmin = hr_uni_lo95, xmax = hr_uni_hi95,
                       y = y_pos)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.6) +
  geom_errorbar(aes(color = color_cat),
                width = 0.25, linewidth = 0.7, na.rm = TRUE,
                orientation = "y") +
  geom_point(aes(shape = pt_shape, size = pt_size, color = color_cat), na.rm = TRUE) +
  geom_text(aes(label = right_label),
            x = x_right, hjust = 0, size = 2.8, family = "sans") +
  scale_shape_identity() +
  scale_size_identity() +
  scale_color_manual(values = c(
    "TCGA-BRCA" = "#2980B9",
    "METABRIC"  = "#E67E22",
    "GSE20685"  = "#27AE60",
    "Pooled"    = "#2C3E50"
  )) +
  scale_x_continuous(
    name   = "Hazard Ratio (per 1 SD Core-PAM score) with 95% CI",
    limits = c(0.85, x_xlim),
    breaks = c(0.9, 1.0, 1.2, 1.4, 1.6, 1.8)
  ) +
  scale_y_continuous(
    breaks = forest_all$y_pos,
    labels = forest_all$label
  ) +
  annotate("segment", x = -Inf, xend = -Inf,
           y = 0.4, yend = nrow(forest_df) + 0.4,
           color = "grey40") +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = -0.5, ymax = 0.5,
           fill = "#ECF0F1", alpha = 0.5) +
  labs(
    title    = "Core-PAM Score — Hazard Ratios Across Validation Cohorts",
    subtitle = "Primary endpoint per cohort | Cox univariate | per 1 SD score",
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

save_fig(p_forest, "FigS_Forest_HR_ValidationCohorts", w = 12, h = 5)

# --------------------------------------------------------------------------
# 2) SCORE DISTRIBUTION — density overlay per cohort
# --------------------------------------------------------------------------
message(sprintf("[%s] === 2) Score distribution ===", SCRIPT_NAME))

load_score <- function(cohort_name) {
  p <- file.path(PATHS$processed, cohort_name, "analysis_ready.parquet")
  df <- arrow::read_parquet(p, col_select = c("score_z"))
  df$cohort <- cohort_name
  df
}

score_dist <- bind_rows(lapply(c("TCGA_BRCA", "METABRIC", "GSE20685"), load_score)) |>
  mutate(
    cohort_label = case_when(
      cohort == "TCGA_BRCA" ~ "TCGA-BRCA (OS, n=1072)",
      cohort == "METABRIC"  ~ "METABRIC (DSS, n=1978)",
      cohort == "GSE20685"  ~ "GSE20685 (OS, n=327)"
    )
  )

pal_dist <- c(
  "TCGA-BRCA (OS, n=1072)" = "#2980B9",
  "METABRIC (DSS, n=1978)" = "#E67E22",
  "GSE20685 (OS, n=327)"  = "#27AE60"
)

p_dist <- ggplot(score_dist, aes(x = score_z, color = cohort_label, fill = cohort_label)) +
  geom_density(alpha = 0.15, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  scale_color_manual(values = pal_dist) +
  scale_fill_manual(values  = pal_dist) +
  labs(
    x        = "Core-PAM Score (z-standardized)",
    y        = "Density",
    title    = "Core-PAM Score Distribution by Cohort",
    subtitle = "Intra-cohort z-standardized; median split line shown",
    color    = NULL, fill = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title      = element_text(face = "bold")
  )

save_fig(p_dist, "FigS_ScoreDistribution_ByCohort", w = 7, h = 5)

# --------------------------------------------------------------------------
# 3) HEATMAP — Core-PAM gene z-scores across validation cohorts
# --------------------------------------------------------------------------
message(sprintf("[%s] === 3) Heatmap ===", SCRIPT_NAME))

core_genes <- c("ACTR3B","BCL2","BLVRA","CENPF","CXXC5","ERBB2","ESR1","EXO1",
                "FGFR4","FOXC1","GPR160","GRB7","KRT17","KRT5","MDM2","MIA",
                "MLPH","MYBL2","MYC","NAT1","PGR","PHGDH","PTTG1","SFRP1")
weights_df <- read_csv(file.path(PATHS$results$corepam, "CorePAM_weights.csv"),
                       show_col_types = FALSE)

load_heatmap_data <- function(cohort_name, endpoint_col, event_col) {
  expr_p <- file.path(PATHS$processed, cohort_name, "expression_genelevel_preZ.parquet")
  rdy_p  <- file.path(PATHS$processed, cohort_name, "analysis_ready.parquet")

  expr_mat <- arrow::read_parquet(expr_p)
  rdy_df   <- arrow::read_parquet(rdy_p)

  genes_in <- intersect(core_genes, expr_mat$gene)
  expr_sub <- expr_mat[expr_mat$gene %in% genes_in, ]

  gene_names   <- expr_sub$gene
  sample_ids   <- setdiff(names(expr_sub), "gene")
  expr_vals    <- t(as.matrix(expr_sub[, sample_ids]))
  colnames(expr_vals) <- gene_names

  # z-score intra-cohort
  old_warn <- getOption("warn"); options(warn = 0)
  z_mat <- scale(expr_vals)
  options(warn = old_warn)

  # Determine which column in rdy_df matches expression sample IDs
  # Try patient_id first, then sample_id (GSE20685 uses sample_id = GSM...)
  id_col <- "patient_id"
  if (length(intersect(rownames(z_mat), rdy_df[[id_col]])) == 0 &&
      "sample_id" %in% names(rdy_df) &&
      length(intersect(rownames(z_mat), rdy_df$sample_id)) > 0) {
    id_col <- "sample_id"
  }

  common_samp <- intersect(rownames(z_mat), rdy_df[[id_col]])
  if (length(common_samp) == 0) {
    stop(sprintf("No matching samples between expression and clinical data (tried %s)", id_col))
  }
  z_sub   <- z_mat[common_samp, , drop = FALSE]
  rdy_sub <- rdy_df[match(common_samp, rdy_df[[id_col]]), ]

  # take a sample of max 200 patients for heatmap visibility
  set.seed(42)
  n_show <- min(200, nrow(z_sub))
  idx    <- sample(nrow(z_sub), n_show)
  z_show <- z_sub[idx, ]
  ann_df <- data.frame(
    score_group = ifelse(rdy_sub$score_z[idx] >= 0, "High", "Low"),
    row.names   = rownames(z_show)
  )
  list(z = z_show, ann = ann_df, cohort = cohort_name)
}

# Attempt heatmap with pheatmap (or ggplot tile if unavailable)
make_heatmap <- function(hdata, title_str, fig_name) {
  z_m <- hdata$z
  ann <- hdata$ann

  # Sort columns by weight
  g_order <- weights_df$gene[weights_df$gene %in% colnames(z_m)]
  g_order <- g_order[order(weights_df$weight[match(g_order, weights_df$gene)])]
  z_m <- z_m[, g_order, drop = FALSE]

  # Sort rows by score group then score
  row_ord <- order(ann$score_group, rownames(z_m))
  z_m <- z_m[row_ord, ]
  ann <- ann[row_ord, , drop = FALSE]

  if (requireNamespace("pheatmap", quietly = TRUE)) {
    ann_colors    <- list(score_group = c(High = COL$km_high, Low = COL$km_low))
    heat_args     <- list(
      t(z_m),
      annotation_col    = ann,
      annotation_colors = ann_colors,
      cluster_rows      = FALSE,
      cluster_cols      = FALSE,
      show_colnames     = FALSE,
      color             = colorRampPalette(c("#2C3E50","white","#C0392B"))(100),
      breaks            = seq(-3, 3, length.out = 101),
      main              = title_str,
      fontsize_row      = 8,
      border_color      = NA
    )
    # PNG via pheatmap's filename arg (avoids S7 dev issues)
    old_warn <- getOption("warn"); options(warn = 0)
    tryCatch({
      do.call(pheatmap::pheatmap,
              c(heat_args, list(filename = file.path(fig_supp, paste0(fig_name, ".png")),
                                width = 10, height = 6)))
    }, error = function(e) {
      message(sprintf("[%s] Heatmap PNG error (%s): %s", SCRIPT_NAME, fig_name, e$message))
    })
    # PDF via grDevices (no ggplot2 S7 involved)
    tryCatch({
      pdf(file.path(fig_supp, paste0(fig_name, ".pdf")), width = 10, height = 6)
      do.call(pheatmap::pheatmap, heat_args)
      dev.off()
    }, error = function(e) {
      try(dev.off(), silent = TRUE)
      message(sprintf("[%s] Heatmap PDF error (%s): %s", SCRIPT_NAME, fig_name, e$message))
    })
    options(warn = old_warn)
    registry_append("ALL", fig_name,
                    file.path(fig_supp, paste0(fig_name, ".png")),
                    sha256_file(file.path(fig_supp, paste0(fig_name, ".png"))),
                    "ok", SCRIPT_NAME,
                    file.info(file.path(fig_supp, paste0(fig_name, ".png")))$size / 1e6)
    message(sprintf("[%s] Heatmap saved: %s", SCRIPT_NAME, fig_name))
  } else {
    # Fallback: ggplot tile heatmap
    z_long <- as.data.frame(z_m) |>
      mutate(sample_id = rownames(z_m),
             score_group = ann$score_group) |>
      pivot_longer(-c(sample_id, score_group), names_to = "gene", values_to = "z")

    z_long$gene <- factor(z_long$gene, levels = g_order)
    z_long$sample_id <- factor(z_long$sample_id, levels = rownames(z_m))

    p_heat <- ggplot(z_long, aes(x = sample_id, y = gene, fill = pmin(pmax(z, -3), 3))) +
      geom_tile() +
      scale_fill_gradient2(low = "#2C3E50", mid = "white", high = "#C0392B",
                           midpoint = 0, limits = c(-3, 3),
                           name = "z-score") +
      facet_grid(. ~ score_group, scales = "free_x", space = "free_x") +
      labs(title = title_str, x = "Samples", y = NULL) +
      theme_classic(base_size = 9) +
      theme(axis.text.x  = element_blank(),
            axis.ticks.x = element_blank(),
            strip.background = element_rect(fill = "#2C3E50"),
            strip.text       = element_text(color = "white", face = "bold"))
    save_fig(p_heat, fig_name, w = 10, h = 6)
  }
}

cohort_configs <- list(
  list(name = "TCGA_BRCA", endpoint = "os_time_months", event = "os_event",
       title = "Core-PAM Gene Expression Heatmap — TCGA-BRCA (OS, n=200 random)"),
  list(name = "METABRIC",  endpoint = "dss_time_months", event = "dss_event",
       title = "Core-PAM Gene Expression Heatmap — METABRIC (DSS, n=200 random)"),
  list(name = "GSE20685",  endpoint = "os_time_months",  event = "os_event",
       title = "Core-PAM Gene Expression Heatmap — GSE20685 (OS, all samples)")
)

for (cfg in cohort_configs) {
  tryCatch({
    hdata <- load_heatmap_data(cfg$name, cfg$endpoint, cfg$event)
    make_heatmap(
      hdata,
      cfg$title,
      sprintf("FigS_Heatmap_%s", cfg$name)
    )
  }, error = function(e) {
    message(sprintf("[%s] Heatmap %s error: %s", SCRIPT_NAME, cfg$name, e$message))
  })
}

# --------------------------------------------------------------------------
# 4) LOLLIPOP CHART — Core-PAM gene weights
# --------------------------------------------------------------------------
message(sprintf("[%s] === 4) Weight lollipop ===", SCRIPT_NAME))

# Assign biological annotation
gene_bio <- tibble(
  gene = c("EXO1","PTTG1","FOXC1","PHGDH","KRT5","MYC","CENPF","ERBB2","GRB7",
           "KRT17","FGFR4","CXXC5",
           "NAT1","BLVRA","ACTR3B","MIA","MYBL2","MDM2","SFRP1","GPR160",
           "ESR1","PGR","BCL2","MLPH"),
  group = c(rep("Proliferation / basal", 12), rep("Luminal / stromal / mixed", 4),
            rep("Cell cycle / apoptosis", 2), "Wnt antagonist", "GPCR",
            rep("ER/PGR signaling", 3), "Luminal marker")
)

wt_df <- weights_df |>
  left_join(gene_bio, by = "gene") |>
  mutate(
    direction = ifelse(weight > 0, "Higher risk (+)", "Lower risk (−)"),
    gene      = reorder(gene, weight)
  )

pal_dir <- c("Higher risk (+)" = "#C0392B", "Lower risk (−)" = "#2980B9")

p_lollipop <- ggplot(wt_df, aes(x = gene, y = weight, color = direction)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.5) +
  geom_segment(aes(x = gene, xend = gene, y = 0, yend = weight),
               linewidth = 0.8) +
  geom_point(aes(size = abs(weight)), shape = 16) +
  geom_text(aes(label = gene, hjust = ifelse(weight > 0, -0.15, 1.15),
                y = weight),
            size = 3, color = "black") +
  scale_color_manual(values = pal_dir, name = "Effect direction") +
  scale_size_continuous(range = c(2, 6), guide = "none") +
  coord_flip() +
  labs(
    x        = NULL,
    y        = "Elastic-net coefficient (weight)",
    title    = "Core-PAM Panel — Gene Weights",
    subtitle = sprintf("%d genes selected by OOF non-inferiority (ΔC=0.010) | SCAN-B training",
                       nrow(wt_df))
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    legend.position = "bottom",
    plot.title      = element_text(face = "bold")
  )

save_fig(p_lollipop, "FigS_Weights_CorePAM_Lollipop", w = 8, h = 6, dir = fig_supp)

# --------------------------------------------------------------------------
# 5) SCORE BY ER STATUS — box/violin per cohort
# --------------------------------------------------------------------------
message(sprintf("[%s] === 5) Score by ER status ===", SCRIPT_NAME))

load_er_score <- function(cohort_name) {
  p <- file.path(PATHS$processed, cohort_name, "analysis_ready.parquet")
  df <- arrow::read_parquet(p)
  cols_keep <- intersect(c("patient_id","score_z","er_status","pam50_subtype"), names(df))
  df <- df[, cols_keep]
  df$cohort <- cohort_name
  df
}

er_data <- bind_rows(lapply(c("TCGA_BRCA","METABRIC","GSE20685"), load_er_score))

if ("er_status" %in% names(er_data)) {
  er_plot <- er_data |>
    filter(!is.na(er_status)) |>
    mutate(
      er_label   = ifelse(as.character(er_status) %in% c("1","Positive","positive","ER+"),
                          "ER+", "ER-"),
      cohort_label = case_when(
        cohort == "TCGA_BRCA" ~ "TCGA-BRCA",
        cohort == "METABRIC"  ~ "METABRIC",
        cohort == "GSE20685"  ~ "GSE20685"
      )
    )

  p_er <- ggplot(er_plot, aes(x = er_label, y = score_z, fill = er_label)) +
    geom_violin(alpha = 0.5, trim = TRUE, linewidth = 0.4) +
    geom_boxplot(width = 0.12, outlier.size = 0.4, linewidth = 0.4, fill = "white") +
    scale_fill_manual(values = c("ER+" = "#2980B9", "ER-" = "#C0392B"), guide = "none") +
    facet_wrap(~ cohort_label, scales = "free_y") +
    labs(
      x        = "ER status",
      y        = "Core-PAM Score (z)",
      title    = "Core-PAM Score by ER Status",
      subtitle = "Validation cohorts | Higher score = higher predicted risk"
    ) +
    theme_classic(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "#2C3E50"),
      strip.text       = element_text(color = "white", face = "bold"),
      plot.title       = element_text(face = "bold")
    )

  save_fig(p_er, "FigS_ScoreByER_ByCohort", w = 9, h = 5, dir = fig_supp)
} else {
  message(sprintf("[%s] er_status column not found — skipping ER violin", SCRIPT_NAME))
}

# --------------------------------------------------------------------------
# 6) C-INDEX COMPARISON — bar chart across cohorts
# --------------------------------------------------------------------------
message(sprintf("[%s] === 6) C-index comparison bar ===", SCRIPT_NAME))

cindex_df <- tab_surv |>
  mutate(
    cohort_label = case_when(
      cohort == "TCGA_BRCA" ~ "TCGA-BRCA\n(OS, 32mo FU)",
      cohort == "METABRIC"  ~ "METABRIC\n(DSS, 159mo FU)",
      cohort == "GSE20685"  ~ "GSE20685\n(OS, 113mo FU)"
    ),
    cohort_label = factor(cohort_label,
                          levels = c("TCGA-BRCA\n(OS, 32mo FU)",
                                     "METABRIC\n(DSS, 159mo FU)",
                                     "GSE20685\n(OS, 113mo FU)"))
  )

p_cindex <- ggplot(cindex_df,
                   aes(x = cohort_label, y = c_index,
                       ymin = c_index_lo95, ymax = c_index_hi95,
                       fill = cohort_label)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50", linewidth = 0.6) +
  geom_col(alpha = 0.8, width = 0.55) +
  geom_errorbar(width = 0.18, linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f\n(%.3f–%.3f)", c_index, c_index_lo95, c_index_hi95),
                y = c_index_hi95 + 0.01),
            size = 3.2, fontface = "bold", vjust = 0) +
  scale_fill_manual(values = c(
    "TCGA-BRCA\n(OS, 32mo FU)"   = "#2980B9",
    "METABRIC\n(DSS, 159mo FU)"  = "#E67E22",
    "GSE20685\n(OS, 113mo FU)"   = "#27AE60"
  ), guide = "none") +
  scale_y_continuous(breaks = seq(0.45, 0.75, 0.05)) +
  coord_cartesian(ylim = c(0.45, 0.75)) +
  labs(
    x        = NULL,
    y        = "Harrell C-index (adjusted, bootstrap 95% CI)",
    title    = "Core-PAM Discriminative Performance by Cohort",
    subtitle = "C_adj = max(C_raw, 1–C_raw) | 1,000 bootstrap resamples | Dashed line = 0.5 (random)"
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

save_fig(p_cindex, "FigS_Cindex_ByCohort", w = 7, h = 5, dir = fig_supp)

# --------------------------------------------------------------------------
# 7) PARETO CURVE — df vs C-index OOF (derivation figure)
# --------------------------------------------------------------------------
message(sprintf("[%s] === 7) Pareto curve ===", SCRIPT_NAME))

pareto_path <- file.path(PATHS$results$corepam, "pareto_df_cindex_oof.csv")
if (file.exists(pareto_path)) {
  pareto <- read_csv(pareto_path, show_col_types = FALSE)

  c_max   <- max(pareto$c_adj, na.rm = TRUE)
  delta_c <- FREEZE$delta_c
  thresh  <- c_max - delta_c

  sel_df  <- fromJSON(file.path(PATHS$results$corepam, "selected_CorePAM_summary.json"))
  n_sel   <- sel_df$selected$df

  for (lang in c("EN", "PT")) {
    x_lbl <- if (lang == "EN") "Number of genes (df)" else "Número de genes (df)"
    y_lbl <- if (lang == "EN") "OOF C-index (adjusted)" else "C-index OOF (ajustado)"
    cap   <- if (lang == "EN") {
      sprintf("Red dashed line: C_max − ΔC = %.3f (non-inferiority threshold)\nSelected: %d genes (▲)",
              thresh, n_sel)
    } else {
      sprintf("Linha vermelha: C_máx − ΔC = %.3f (limiar de não-inferioridade)\nSelecionado: %d genes (▲)",
              thresh, n_sel)
    }

    selected_row <- pareto[pareto$df == n_sel, ]

    p_pareto <- ggplot(pareto, aes(x = df, y = c_adj)) +
      geom_line(color = "steelblue", linewidth = 0.8) +
      geom_point(size = 1.5, color = "steelblue") +
      geom_hline(yintercept = thresh, linetype = "dashed",
                 color = "#C0392B", linewidth = 0.7) +
      geom_hline(yintercept = c_max, linetype = "dotted",
                 color = "grey40", linewidth = 0.5) +
      geom_point(data = selected_row,
                 aes(x = df, y = c_adj),
                 shape = 17, size = 5, color = "#C0392B") +
      annotate("text",
               x = n_sel + 1, y = selected_row$c_adj + 0.003,
               label = sprintf("Selected\ndf=%d", n_sel),
               hjust = 0, size = 3, color = "#C0392B") +
      scale_x_continuous(breaks = seq(0, 50, 5)) +
      labs(x = x_lbl, y = y_lbl,
           title    = if (lang == "EN") "Core-PAM Derivation: Pareto Curve (df vs C-index OOF)"
                      else "Derivação Core-PAM: Curva Pareto (df × C-index OOF)",
           subtitle = if (lang == "EN") "SCAN-B training cohort (n=3069) | Elastic-net Cox (α=0.5, K=10)"
                      else "Coorte de treinamento SCAN-B (n=3.069) | Cox elastic-net (α=0,5, K=10)",
           caption  = cap) +
      theme_classic(base_size = 11) +
      theme(plot.title   = element_text(face = "bold"),
            plot.caption = element_text(size = 8))

    fig_name <- sprintf("FigS_Pareto_CorePAM_%s", lang)
    save_fig(p_pareto, fig_name, w = 8, h = 5, dir = fig_supp)
  }
} else {
  message(sprintf("[%s] Pareto CSV not found — skip", SCRIPT_NAME))
}

message(sprintf("[%s] COMPLETED", SCRIPT_NAME))
