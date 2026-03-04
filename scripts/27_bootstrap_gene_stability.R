# =============================================================================
# SCRIPT: 27_bootstrap_gene_stability.R
# PURPOSE: Bootstrap gene stability analysis — refits Cox elastic-net B=200
#          times on SCAN-B resamples and reports how often each PAM50 gene
#          appears with non-zero coefficient at the frozen lambda.
#          Generates FigS_Bootstrap_GeneFreq (bar plot, 50 genes).
#
# INPUTS:
#   01_Base_Pura_CorePAM/PROCESSED/SCANB/expression_genelevel_preZ.parquet
#   01_Base_Pura_CorePAM/PROCESSED/SCANB/clinical_FINAL.parquet
#   results/corepam/CorePAM_weights.csv
#   results/corepam/selected_CorePAM_summary.json
#
# OUTPUTS:
#   results/supp/bootstrap_gene_stability.csv
#   figures/supp/en/png/FigS_Bootstrap_GeneFreq_EN.png
#   figures/supp/en/pdf/FigS_Bootstrap_GeneFreq_EN.pdf
#
# PROJECT: Core-PAM (Memorial v6.1 / Freeze Core-PAM)
# =============================================================================

source("scripts/00_setup.R")

SCRIPT_NAME <- "27_bootstrap_gene_stability.R"
COHORT      <- "SCANB"

# Skip/force pattern
FORCE <- as.logical(Sys.getenv("FORCE_RERUN", "FALSE"))
out_csv <- file.path(PATHS$results$supp, "bootstrap_gene_stability.csv")
if (!FORCE && file.exists(out_csv)) {
  message(sprintf("[%s] Output exists — skipping. Set FORCE_RERUN=TRUE.", SCRIPT_NAME))
  quit(save = "no", status = 0)
}

old_warn_pkg <- getOption("warn"); options(warn = 0)
suppressPackageStartupMessages({
  library(glmnet)
  library(survival)
  library(ggplot2)
})
options(warn = old_warn_pkg)

# =============================================================================
# 1) FROZEN PARAMETERS
# =============================================================================
ALPHA_EN <- FREEZE$alpha  # 0.5

# Load frozen lambda from training card
card <- jsonlite::fromJSON(
  file.path(PATHS$results$corepam, "selected_CorePAM_summary.json")
)
LAMBDA_FROZEN <- card$selected$lambda
message(sprintf("[%s] Frozen lambda = %.6g | alpha = %.1f", SCRIPT_NAME, LAMBDA_FROZEN, ALPHA_EN))

# =============================================================================
# 2) PAM50 GENE LIST (canonical 50 genes)
# =============================================================================
PAM50_GENES <- c(
  "ACTR3B","ANLN","BAG1","BCL2","BIRC5","BLVRA","CCNB1","CCNE1","CDC20","CDC6",
  "CDH3","CENPF","CEP55","CXXC5","EGFR","ERBB2","ESR1","EXO1","FGFR4","FOXA1",
  "FOXC1","GPR160","GRB7","KIF2C","KRT14","KRT17","KRT5","MAPT","MDM2",
  "MELK","MIA","MKI67","MLPH","MMP11","MYBL2","MYC","NAT1","NDC80","NUF2",
  "ORC6","PGR","PHGDH","PTTG1","RRM2","SFRP1","SLC39A6","TMEM45B","TYMS","UBE2C","UBE2T"
)

# =============================================================================
# 3) LOAD CorePAM WEIGHTS (the 24 selected genes)
# =============================================================================
weights_df <- strict_csv(file.path(PATHS$results$corepam, "CorePAM_weights.csv"))
corepam_genes <- weights_df$gene
message(sprintf("[%s] CorePAM panel: %d genes", SCRIPT_NAME, length(corepam_genes)))

# =============================================================================
# 4) LOAD EXPRESSION DATA (SCAN-B — PAM50 genes only)
# =============================================================================
message(sprintf("[%s] Loading expression data...", SCRIPT_NAME))
path_expr <- file.path(proc_cohort(COHORT), "expression_genelevel_preZ.parquet")
old_warn_arrow <- getOption("warn"); options(warn = 0)
expr_full <- arrow::read_parquet(path_expr)
pam50_df  <- expr_full[expr_full$gene %in% PAM50_GENES, ]
rm(expr_full)
options(warn = old_warn_arrow)

gene_names_found <- pam50_df$gene
sample_ids       <- setdiff(names(pam50_df), "gene")
X_raw <- t(as.matrix(pam50_df[, sample_ids]))
colnames(X_raw) <- gene_names_found
message(sprintf("[%s] Expression: %d samples x %d PAM50 genes", SCRIPT_NAME,
                nrow(X_raw), ncol(X_raw)))

# =============================================================================
# 5) LOAD CLINICAL + JOIN
# =============================================================================
path_clin  <- file.path(proc_cohort(COHORT), "clinical_FINAL.parquet")
clin       <- strict_parquet(path_clin)
common_ids <- intersect(clin$patient_id, rownames(X_raw))
clin_sub   <- clin[match(common_ids, clin$patient_id), ]
X_sub      <- X_raw[common_ids, , drop = FALSE]

os_time   <- clin_sub$os_time_months
os_event  <- as.integer(clin_sub$os_event)
n_samples <- length(common_ids)
n_events  <- sum(os_event)

message(sprintf("[%s] Training data: N=%d | events=%d (%.1f%%)",
                SCRIPT_NAME, n_samples, n_events, 100 * n_events / n_samples))

# Scale X (same as script 05)
old_warn_sc <- getOption("warn"); options(warn = 0)
X_scaled <- scale(X_sub)
options(warn = old_warn_sc)

y_surv <- survival::Surv(os_time, os_event)

# =============================================================================
# 6) BOOTSTRAP LOOP (B=200)
# =============================================================================
B <- 200
set.seed(42)

# Matrix to track which genes appear in each bootstrap
gene_presence <- matrix(FALSE, nrow = B, ncol = ncol(X_scaled),
                        dimnames = list(NULL, colnames(X_scaled)))

message(sprintf("[%s] Starting %d bootstrap refits...", SCRIPT_NAME, B))

for (b in seq_len(B)) {
  # Bootstrap resample (with replacement)
  idx <- sample(n_samples, n_samples, replace = TRUE)
  X_b <- X_scaled[idx, , drop = FALSE]
  y_b <- y_surv[idx]

  old_w <- getOption("warn"); options(warn = 0)
  fit_b <- tryCatch(
    glmnet::glmnet(
      x           = X_b,
      y           = y_b,
      family      = "cox",
      alpha       = ALPHA_EN,
      lambda      = LAMBDA_FROZEN,
      standardize = FALSE
    ),
    error = function(e) NULL
  )
  options(warn = old_w)

  if (!is.null(fit_b)) {
    coefs <- as.matrix(coef(fit_b, s = LAMBDA_FROZEN))
    nonzero <- rownames(coefs)[as.numeric(coefs) != 0]
    gene_presence[b, nonzero] <- TRUE
  }

  if (b %% 50 == 0 || b == B) {
    message(sprintf("[%s]   Bootstrap %d/%d done", SCRIPT_NAME, b, B))
  }
}

# =============================================================================
# 7) COMPUTE FREQUENCIES AND SAVE
# =============================================================================
gene_freq <- colSums(gene_presence)
freq_pct  <- 100 * gene_freq / B

# Build results for all 50 PAM50 genes (including those with freq=0)
all_genes_in_data <- colnames(X_scaled)
result_df <- tibble::tibble(
  gene          = all_genes_in_data,
  boot_count    = as.integer(gene_freq[all_genes_in_data]),
  boot_freq_pct = round(freq_pct[all_genes_in_data], 1),
  in_corepam    = all_genes_in_data %in% corepam_genes
) |> dplyr::arrange(dplyr::desc(boot_freq_pct))

readr::write_csv(result_df, out_csv)
message(sprintf("[%s] Saved: %s (%d genes)", SCRIPT_NAME, basename(out_csv), nrow(result_df)))

# =============================================================================
# 8) FIGURE: Bar plot of bootstrap frequency
# =============================================================================
# Order genes by frequency
result_df$gene <- factor(result_df$gene, levels = rev(result_df$gene))

p <- ggplot(result_df, aes(x = gene, y = boot_freq_pct,
                           fill = in_corepam)) +
  geom_col(width = 0.75) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  coord_flip() +
  scale_fill_manual(
    values = c("TRUE" = "#762A83", "FALSE" = "#B8B8B8"),
    labels = c("TRUE" = "CorePAM (24 genes)", "FALSE" = "Excluded PAM50"),
    name   = NULL
  ) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title    = "Bootstrap Gene Selection Stability (B=200)",
    subtitle = sprintf("Cox Elastic-Net (alpha=%.1f) refits on SCAN-B (N=%d)",
                        ALPHA_EN, n_samples),
    x        = NULL,
    y        = "Selection Frequency (%)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position  = "bottom",
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, color = "grey40"),
    axis.text.y      = element_text(size = 7)
  )

# Save
out_png <- file.path(PATHS$figures$supp_en_png, "FigS_Bootstrap_GeneFreq_EN.png")
out_pdf <- file.path(PATHS$figures$supp_en_pdf, "FigS_Bootstrap_GeneFreq_EN.pdf")

old_warn_gs <- getOption("warn"); options(warn = 0)
ggsave(out_png, p, width = 170, height = 200, units = "mm", dpi = 300, bg = "white")
ggsave(out_pdf, p, width = 170, height = 200, units = "mm", device = cairo_pdf)
options(warn = old_warn_gs)

message(sprintf("[%s] Figure saved: %s", SCRIPT_NAME, basename(out_png)))

# =============================================================================
# 9) REGISTRY
# =============================================================================
for (f in c(out_csv, out_png)) {
  h  <- sha256_file(f)
  sz <- file.info(f)$size / 1e6
  registry_append(COHORT, "bootstrap_stability", f, h, "ok", SCRIPT_NAME, sz, list())
}

# Summary
n_core_above_80 <- sum(result_df$boot_freq_pct[result_df$in_corepam] >= 80)
n_core_above_50 <- sum(result_df$boot_freq_pct[result_df$in_corepam] >= 50)
message(sprintf("[%s] CorePAM genes >= 80%% frequency: %d/%d",
                SCRIPT_NAME, n_core_above_80, sum(result_df$in_corepam)))
message(sprintf("[%s] CorePAM genes >= 50%% frequency: %d/%d",
                SCRIPT_NAME, n_core_above_50, sum(result_df$in_corepam)))
message(sprintf("[%s] DONE", SCRIPT_NAME))
