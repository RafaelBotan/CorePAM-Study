# =============================================================================
# SCRIPT: 13_qc_correlations_offdiag.R
# PURPOSE: Correlacoes entre scores CorePAM das coortes EXCLUINDO diagonal
#          (sem autocorrelacao self=1). Figura FigS3.
#          Segue Memorial v6.1 §1.2 (sem pooling; correlacao apenas de scores).
# PROJETO: Core-PAM (Memorial v6.1 / Freeze Core-PAM)
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "13_qc_correlations_offdiag.R"

suppressPackageStartupMessages({
  library(ggplot2)
})

message(sprintf("[%s] Iniciando QC de correlacoes entre scores (off-diagonal)", SCRIPT_NAME))

COHORTS <- c("SCANB", "GSE96058", "TCGA_BRCA", "METABRIC", "GSE20685")

# --------------------------------------------------------------------------
# 1) Carregar score por coorte (score_z de analysis_ready.parquet)
# --------------------------------------------------------------------------
score_list <- list()

for (coh in COHORTS) {
  ready_path <- file.path(proc_cohort(coh), "analysis_ready.parquet")
  if (!file.exists(ready_path)) {
    message(sprintf("[%s] AVISO: %s nao encontrado. Pulando %s.", SCRIPT_NAME, ready_path, coh))
    next
  }
  df <- strict_parquet(ready_path)
  if (!"score_z" %in% names(df)) {
    message(sprintf("[%s] %s: coluna score_z ausente. Pulando.", SCRIPT_NAME, coh))
    next
  }
  # Usar sample_id como chave
  if (!"sample_id" %in% names(df)) {
    message(sprintf("[%s] %s: sample_id ausente. Pulando.", SCRIPT_NAME, coh))
    next
  }
  score_list[[coh]] <- df[, c("sample_id", "score_z")]
  names(score_list[[coh]])[2] <- coh
  message(sprintf("[%s] %s: %d amostras com score_z carregadas", SCRIPT_NAME, coh, nrow(df)))
}

if (length(score_list) < 2) {
  stop(sprintf("[%s] Menos de 2 coortes com dados. Impossivel calcular correlacoes.", SCRIPT_NAME))
}

# --------------------------------------------------------------------------
# 2) Calcular correlacoes par-a-par (off-diagonal; os IDs nao se repetem entre coortes)
#    Estrategia: correlacao da distribuicao de scores dentro de cada coorte
#    (transformacao quantil-quantil ou Spearman entre percentis)
#    Nota: amostras nao se sobrepoem entre coortes (principio Zero-Pooling).
#    Correlamos os quantis empiricos (percentis) dos scores entre coortes.
# --------------------------------------------------------------------------
cohort_names <- names(score_list)
n_coh        <- length(cohort_names)

# Calcular percentis uniformes para cada coorte
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

# Matriz de correlacoes Spearman entre percentis das coortes
cor_mat <- matrix(NA_real_, nrow = n_coh, ncol = n_coh,
                  dimnames = list(cohort_names, cohort_names))

for (i in seq_len(n_coh)) {
  for (j in seq_len(n_coh)) {
    if (i == j) {
      cor_mat[i, j] <- NA_real_  # Diagonal excluida (off-diagonal)
    } else {
      old_warn <- getOption("warn"); options(warn = 0)
      cor_val <- tryCatch(
        cor(pct_mat[, i], pct_mat[, j], method = "spearman", use = "complete.obs"),
        error = function(e) NA_real_
      )
      options(warn = old_warn)
      cor_mat[i, j] <- cor_val
    }
  }
}

message(sprintf("[%s] Matriz de correlacao off-diagonal calculada:", SCRIPT_NAME))
print(round(cor_mat, 3))

# --------------------------------------------------------------------------
# 3) Salvar matriz de correlacoes
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
# 4) Figura FigS3 — heatmap de correlacoes off-diagonal
# --------------------------------------------------------------------------
# Converter para formato longo para ggplot
cor_long <- reshape2::melt(cor_mat, varnames = c("Coorte_A", "Coorte_B"),
                            value.name = "Spearman_rho", na.rm = TRUE)

# Garantir reshape2 disponivel
if (!requireNamespace("reshape2", quietly = TRUE)) {
  cor_long <- do.call(rbind, lapply(cohort_names, function(i) {
    do.call(rbind, lapply(cohort_names, function(j) {
      if (i != j) {
        data.frame(Coorte_A = i, Coorte_B = j, Spearman_rho = cor_mat[i, j],
                   stringsAsFactors = FALSE)
      }
    }))
  }))
}

supp_fig_dir <- PATHS$figures$supp
dir.create(supp_fig_dir, showWarnings = FALSE, recursive = TRUE)

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
    title    = "Correlacoes entre scores CorePAM (off-diagonal, via percentis)",
    subtitle = "Diagonal excluida (sem autocorrelacao self=1)",
    x        = NULL,
    y        = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    plot.title   = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray40")
  )

figs3_pdf <- file.path(supp_fig_dir, "FigS3_Correlation_OffDiagonal.pdf")
figs3_png <- file.path(supp_fig_dir, "FigS3_Correlation_OffDiagonal.png")

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

message(sprintf("[%s] CONCLUIDO | FigS3 salva", SCRIPT_NAME))
