# =============================================================================
# SCRIPT: 28_dropgene_sensitivity.R
# PURPOSE: Drop-one-gene sensitivity analysis — for each of the 24 CorePAM
#          genes, removes it from the panel, recomputes the score using the
#          remaining 23 genes (original weights, renormalized), and measures
#          delta-C-index in each validation cohort.
#          Generates FigS_DropGene_Tornado (lollipop/tornado plot).
#
# INPUTS:
#   results/corepam/CorePAM_weights.csv
#   01_Base_Pura_CorePAM/PROCESSED/{COHORT}/analysis_ready.parquet  (for baseline C)
#   01_Base_Pura_CorePAM/PROCESSED/{COHORT}/expression_genelevel_preZ.parquet
#   01_Base_Pura_CorePAM/PROCESSED/{COHORT}/clinical_FINAL.parquet
#
# OUTPUTS:
#   results/supp/dropgene_sensitivity.csv
#   figures/supp/en/png/FigS_DropGene_Tornado_EN.png
#   figures/supp/en/pdf/FigS_DropGene_Tornado_EN.pdf
#
# PROJECT: Core-PAM (Memorial v6.1 / Freeze Core-PAM)
# =============================================================================

source("scripts/00_setup.R")

SCRIPT_NAME <- "28_dropgene_sensitivity.R"

# Skip/force pattern
FORCE <- as.logical(Sys.getenv("FORCE_RERUN", "FALSE"))
out_csv <- file.path(PATHS$results$supp, "dropgene_sensitivity.csv")
if (!FORCE && file.exists(out_csv)) {
  message(sprintf("[%s] Output exists — skipping. Set FORCE_RERUN=TRUE.", SCRIPT_NAME))
  quit(save = "no", status = 0)
}

old_warn_pkg <- getOption("warn"); options(warn = 0)
suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
})
options(warn = old_warn_pkg)

# =============================================================================
# 1) LOAD CorePAM WEIGHTS
# =============================================================================
weights_df <- strict_csv(file.path(PATHS$results$corepam, "CorePAM_weights.csv"))
panel_genes   <- weights_df$gene
panel_weights <- setNames(weights_df$weight, weights_df$gene)
n_panel <- length(panel_genes)
message(sprintf("[%s] CorePAM panel: %d genes", SCRIPT_NAME, n_panel))

# =============================================================================
# 2) DEFINE VALIDATION COHORTS
# =============================================================================
COHORTS <- c("TCGA_BRCA", "METABRIC", "GSE20685", "GSE1456")

# METABRIC uses DSS as primary endpoint, others use OS
get_endpoint <- function(cohort) {
  if (cohort == "METABRIC") list(time = "dss_time_months", event = "dss_event")
  else list(time = "os_time_months", event = "os_event")
}

# =============================================================================
# 3) HELPER: Compute CorePAM score for a given gene subset
# =============================================================================
compute_score <- function(z_mat, genes_use, weights_all) {
  # Score = sum(w_i * z_i) / sum(|w_i|) for genes_use subset
  w_use <- weights_all[genes_use]
  denom <- sum(abs(w_use))
  if (denom == 0) return(rep(NA_real_, nrow(z_mat)))

  z_sub <- z_mat[, genes_use, drop = FALSE]
  raw   <- as.numeric(z_sub %*% w_use) / denom
  as.numeric(scale(raw))  # standardize (score_z)
}

# =============================================================================
# 4) HELPER: Harrell C-index (direction-adjusted: always >= 0.5)
# =============================================================================
harrell_c <- function(time, status, lp) {
  old_w <- getOption("warn"); options(warn = 0)
  c_raw <- as.numeric(survival::concordance(survival::Surv(time, status) ~ lp)$concordance)
  options(warn = old_w)
  max(c_raw, 1 - c_raw)
}

# =============================================================================
# 5) MAIN LOOP: for each cohort, compute baseline C and drop-gene delta-C
# =============================================================================
results_list <- list()

for (cohort in COHORTS) {
  message(sprintf("[%s] Processing %s...", SCRIPT_NAME, cohort))

  ep <- get_endpoint(cohort)

  # Load expression (genes x samples format)
  expr_path <- file.path(proc_cohort(cohort), "expression_genelevel_preZ.parquet")
  old_w <- getOption("warn"); options(warn = 0)
  expr_raw <- arrow::read_parquet(expr_path)
  options(warn = old_w)

  gene_names <- expr_raw[["gene"]]
  sample_ids <- setdiff(names(expr_raw), "gene")

  # Transpose to samples x genes
  expr_vals <- t(as.matrix(expr_raw[, sample_ids]))
  colnames(expr_vals) <- gene_names

  # Intra-cohort Z-score
  old_w <- getOption("warn"); options(warn = 0)
  gene_means <- colMeans(expr_vals, na.rm = TRUE)
  gene_sds   <- apply(expr_vals, 2, sd, na.rm = TRUE)
  options(warn = old_w)
  zero_sd <- gene_sds == 0 | is.na(gene_sds)
  gene_sds[zero_sd] <- 1
  z_mat <- sweep(sweep(expr_vals, 2, gene_means, "-"), 2, gene_sds, "/")
  z_mat[, zero_sd] <- 0

  # Load clinical
  clin_path <- file.path(proc_cohort(cohort), "clinical_FINAL.parquet")
  clin <- strict_parquet(clin_path)

  common_ids <- intersect(clin$patient_id, rownames(z_mat))
  # Some cohorts may use sample_id as key
  if (length(common_ids) < 10 && "sample_id" %in% names(clin)) {
    common_ids <- intersect(clin$sample_id, rownames(z_mat))
    clin_sub <- clin[match(common_ids, clin$sample_id), ]
  } else {
    clin_sub <- clin[match(common_ids, clin$patient_id), ]
  }
  z_sub <- z_mat[common_ids, , drop = FALSE]

  time_col  <- clin_sub[[ep$time]]
  event_col <- as.integer(clin_sub[[ep$event]])

  # Filter valid observations
  valid <- !is.na(time_col) & !is.na(event_col) & time_col > 0
  time_col  <- time_col[valid]
  event_col <- event_col[valid]
  z_sub     <- z_sub[valid, , drop = FALSE]

  n_cohort <- length(time_col)
  n_events <- sum(event_col)
  message(sprintf("[%s]   %s: N=%d | events=%d", SCRIPT_NAME, cohort, n_cohort, n_events))

  # Genes available in this cohort
  genes_present_full <- intersect(panel_genes, colnames(z_sub))

  # Baseline C-index (all present genes) — c_adj (direction-adjusted)
  score_full <- compute_score(z_sub, genes_present_full, panel_weights)
  c_full     <- harrell_c(time_col, event_col, score_full)

  message(sprintf("[%s]   %s baseline C_adj = %.4f (%d genes present)",
                  SCRIPT_NAME, cohort, c_full, length(genes_present_full)))

  # Drop each gene
  for (gene_drop in panel_genes) {
    genes_reduced <- setdiff(genes_present_full, gene_drop)

    if (length(genes_reduced) == length(genes_present_full)) {
      # Gene was not present in this cohort anyway
      delta_c <- NA_real_
      c_drop  <- NA_real_
    } else {
      score_drop <- compute_score(z_sub, genes_reduced, panel_weights)
      c_drop  <- harrell_c(time_col, event_col, score_drop)
      delta_c <- c_drop - c_full
    }

    results_list[[length(results_list) + 1]] <- tibble::tibble(
      cohort     = cohort,
      gene       = gene_drop,
      c_full     = round(c_full, 4),
      c_drop     = round(c_drop, 4),
      delta_c    = round(delta_c, 4),
      n_samples  = n_cohort,
      n_events   = n_events,
      gene_present = gene_drop %in% colnames(z_sub)
    )
  }
}

# =============================================================================
# 6) SAVE RESULTS
# =============================================================================
result_df <- dplyr::bind_rows(results_list)
readr::write_csv(result_df, out_csv)
message(sprintf("[%s] Saved: %s (%d rows)", SCRIPT_NAME, basename(out_csv), nrow(result_df)))

# =============================================================================
# 7) FIGURE: Tornado/lollipop plot of delta-C per gene, faceted by cohort
# =============================================================================
plot_df <- result_df |>
  dplyr::filter(!is.na(delta_c)) |>
  dplyr::mutate(
    gene = factor(gene, levels = rev(sort(unique(gene)))),
    color = ifelse(delta_c < 0, "Decrease", "Increase"),
    # Clean cohort label
    cohort_label = dplyr::case_when(
      cohort == "TCGA_BRCA" ~ "TCGA-BRCA",
      cohort == "METABRIC"  ~ "METABRIC",
      cohort == "GSE20685"  ~ "GSE20685",
      cohort == "GSE1456"   ~ "GSE1456",
      TRUE ~ cohort
    )
  )

p <- ggplot(plot_df, aes(x = delta_c, y = gene, color = color)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey50", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = delta_c, y = gene, yend = gene), linewidth = 0.6) +
  geom_point(size = 2) +
  facet_wrap(~ cohort_label, ncol = 2, scales = "free_x") +
  scale_color_manual(
    values = c("Decrease" = "#D6604D", "Increase" = "#2166AC"),
    name   = NULL,
    guide  = "none"
  ) +
  labs(
    title    = "Drop-One-Gene Sensitivity Analysis",
    subtitle = "Change in C-index when removing each CorePAM gene",
    x        = expression(Delta * "C-index (drop-one vs full panel)"),
    y        = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    strip.text       = element_text(face = "bold", size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, color = "grey40"),
    axis.text.y      = element_text(size = 7)
  )

out_png <- file.path(PATHS$figures$supp_en_png, "FigS_DropGene_Tornado_EN.png")
out_pdf <- file.path(PATHS$figures$supp_en_pdf, "FigS_DropGene_Tornado_EN.pdf")

old_warn_gs <- getOption("warn"); options(warn = 0)
ggsave(out_png, p, width = 200, height = 200, units = "mm", dpi = 300, bg = "white")
ggsave(out_pdf, p, width = 200, height = 200, units = "mm", device = cairo_pdf)
options(warn = old_warn_gs)

message(sprintf("[%s] Figure saved: %s", SCRIPT_NAME, basename(out_png)))

# =============================================================================
# 8) REGISTRY
# =============================================================================
for (f in c(out_csv, out_png)) {
  h  <- sha256_file(f)
  sz <- file.info(f)$size / 1e6
  registry_append("MULTI", "dropgene_sensitivity", f, h, "ok", SCRIPT_NAME, sz, list())
}

# Summary
mean_delta <- mean(abs(result_df$delta_c), na.rm = TRUE)
max_delta  <- max(abs(result_df$delta_c), na.rm = TRUE)
message(sprintf("[%s] Mean |delta-C| = %.4f | Max |delta-C| = %.4f",
                SCRIPT_NAME, mean_delta, max_delta))
message(sprintf("[%s] DONE", SCRIPT_NAME))
