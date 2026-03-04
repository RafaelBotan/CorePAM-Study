# =============================================================================
# SCRIPT: 13_qc_correlations_offdiag.R
# PURPOSE: Correlations between CorePAM scores across cohorts EXCLUDING diagonal
#          (no self-autocorrelation self=1). Figure FigS3.
#          Follows Memorial v6.1 sec.1.2 (no pooling; score correlations only).
# PROJECT: Core-PAM (Memorial v6.1 / Freeze Core-PAM)
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "13_qc_correlations_offdiag.R"

suppressPackageStartupMessages({
  library(ggplot2)
})

message(sprintf("[%s] Starting QC of score correlations (off-diagonal)", SCRIPT_NAME))

COHORTS <- c("SCANB", "TCGA_BRCA", "METABRIC", "GSE20685", "GSE1456")

# --------------------------------------------------------------------------
# 1) Load score per cohort (score_z from analysis_ready.parquet)
# --------------------------------------------------------------------------
score_list <- list()

for (coh in COHORTS) {
  ready_path <- file.path(proc_cohort(coh), "analysis_ready.parquet")
  if (!file.exists(ready_path)) {
    message(sprintf("[%s] WARNING: %s not found. Skipping %s.", SCRIPT_NAME, ready_path, coh))
    next
  }
  df <- strict_parquet(ready_path)
  if (!"score_z" %in% names(df)) {
    message(sprintf("[%s] %s: score_z column missing. Skipping.", SCRIPT_NAME, coh))
    next
  }
  # Use sample_id as key
  id_col <- if ("sample_id" %in% names(df)) "sample_id" else if ("patient_id" %in% names(df)) "patient_id" else NULL
  if (is.null(id_col)) {
    message(sprintf("[%s] %s: no sample_id or patient_id. Skipping.", SCRIPT_NAME, coh))
    next
  }
  score_list[[coh]] <- df[, c(id_col, "score_z")]
  names(score_list[[coh]])[2] <- coh
  message(sprintf("[%s] %s: %d samples with score_z loaded", SCRIPT_NAME, coh, nrow(df)))
}

if (length(score_list) < 2) {
  stop(sprintf("[%s] Fewer than 2 cohorts with data. Cannot compute correlations.", SCRIPT_NAME))
}

# --------------------------------------------------------------------------
# 2) Compute pairwise correlations (off-diagonal; IDs do not overlap between cohorts)
#    Strategy: correlate score distributions within each cohort
#    (quantile-quantile transformation or Spearman between percentiles)
#    Note: samples do not overlap across cohorts (Zero-Pooling principle).
#    Correlate empirical quantiles (percentiles) of scores across cohorts.
# --------------------------------------------------------------------------
cohort_names <- names(score_list)
n_coh        <- length(cohort_names)

# Compute uniform percentiles for each cohort
n_pts   <- 100
pct_seq <- seq(0.01, 0.99, length.out = n_pts)

pct_mat <- matrix(NA_real_, nrow = n_pts, ncol = n_coh,
                  dimnames = list(NULL, cohort_names))

for (coh in cohort_names) {
  scores <- score_list[[coh]][[coh]]
  scores <- scores[!is.na(scores)]
  old_warn <- getOption("warn"); options(warn = 0)
  pct_mat[, coh] <- quantile(scores, probs = pct_seq, na.rm = TRUE)
  options(warn = old_warn)
}

# Spearman correlation matrix between cohort percentiles
cor_mat <- matrix(NA_real_, nrow = n_coh, ncol = n_coh,
                  dimnames = list(cohort_names, cohort_names))

for (i in seq_len(n_coh)) {
  for (j in seq_len(n_coh)) {
    if (i == j) {
      cor_mat[i, j] <- NA_real_  # Diagonal excluded (off-diagonal)
    } else {
      old_warn <- getOption("warn"); options(warn = 0)
      cor_val <- tryCatch(
        cor(pct_mat[, i], pct_mat[, j], method = "pearson", use = "complete.obs"),
        error = function(e) NA_real_
      )
      options(warn = old_warn)
      cor_mat[i, j] <- cor_val
    }
  }
}

message(sprintf("[%s] Off-diagonal correlation matrix computed:", SCRIPT_NAME))
print(round(cor_mat, 3))

# --------------------------------------------------------------------------
# 3) Save correlation matrix
# --------------------------------------------------------------------------
cor_df <- as.data.frame(cor_mat)
cor_df$cohort_row <- rownames(cor_mat)
cor_df <- cor_df[, c("cohort_row", cohort_names)]

supp_dir <- PATHS$results$supp
dir.create(supp_dir, showWarnings = FALSE, recursive = TRUE)

cor_path <- file.path(supp_dir, "qc_score_correlations_offdiag.csv")
readr::write_csv(cor_df, cor_path)
h_cor <- sha256_file(cor_path)
registry_append("ALL", "qc_correlations_offdiag", cor_path, h_cor, "ok", SCRIPT_NAME,
                file.info(cor_path)$size / 1e6)

# --------------------------------------------------------------------------
# 4) Figure FigS3 — off-diagonal correlation heatmap
# --------------------------------------------------------------------------
# Convert to long format for ggplot
cor_long <- reshape2::melt(cor_mat, varnames = c("Cohort_A", "Cohort_B"),
                            value.name = "Spearman_rho", na.rm = TRUE)

# Ensure reshape2 is available
if (!requireNamespace("reshape2", quietly = TRUE)) {
  cor_long <- do.call(rbind, lapply(cohort_names, function(i) {
    do.call(rbind, lapply(cohort_names, function(j) {
      if (i != j) {
        data.frame(Cohort_A = i, Cohort_B = j, Spearman_rho = cor_mat[i, j],
                   stringsAsFactors = FALSE)
      }
    }))
  }))
}

# figure dirs created by 00_setup.R

p_cor <- ggplot(cor_long, aes(x = Cohort_A, y = Cohort_B, fill = Spearman_rho)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", Spearman_rho)),
            color = "black", size = 3.5, na.rm = TRUE) +
  scale_fill_gradient2(
    low  = "#2166AC", mid = "white", high = "#D6604D",
    midpoint = 0, limits = c(-1, 1), na.value = "gray90",
    name = "Pearson r\n(quantile values)"
  ) +
  labs(
    title    = "CorePAM score correlations (off-diagonal, via quantile values)",
    subtitle = "Pearson on quantile values | Diagonal excluded",
    x        = NULL,
    y        = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    plot.title   = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray40")
  )

figs3_pdf <- file.path(PATHS$figures$supp_en_pdf, "FigS3_Correlation_OffDiagonal.pdf")
figs3_png <- file.path(PATHS$figures$supp_en_png, "FigS3_Correlation_OffDiagonal.png")

old_warn <- getOption("warn"); options(warn = 0)
pdf(figs3_pdf, width = 7, height = 6); print(p_cor); dev.off()
png(figs3_png, width = 700, height = 600, res = 100); print(p_cor); dev.off()
options(warn = old_warn)

h_s3_pdf <- sha256_file(figs3_pdf)
h_s3_png <- sha256_file(figs3_png)
registry_append("ALL", "figure_qc_corr_offdiag", figs3_pdf, h_s3_pdf, "ok", SCRIPT_NAME,
                file.info(figs3_pdf)$size / 1e6)
registry_append("ALL", "figure_qc_corr_offdiag_png", figs3_png, h_s3_png, "ok", SCRIPT_NAME,
                file.info(figs3_png)$size / 1e6)

message(sprintf("[%s] COMPLETED | FigS3 saved", SCRIPT_NAME))
