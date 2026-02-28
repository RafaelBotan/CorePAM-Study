# =============================================================================
# SCRIPT: 06_zscore_and_score_GSE20685.R
# PURPOSE: Aplica Core-PAM score na coorte GSE20685 (Taiwan):
#          Z-score intra-coorte, score re-escalonado (Effective Gene Count),
#          direcao padronizada (HR>1), gera analysis_ready.parquet
# COORTE:  GSE20685 (validacao microarray; Taiwan)
# PROJETO: Core-PAM (Memorial v6.1 / Freeze Core-PAM)
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "06_zscore_and_score_GSE20685.R"
COHORT      <- "GSE20685"

suppressPackageStartupMessages({
  library(survival)
})

message(sprintf("[%s] Iniciando score CorePAM para %s", SCRIPT_NAME, COHORT))

# --------------------------------------------------------------------------
# 1) Carregar pesos congelados do Core-PAM
# --------------------------------------------------------------------------
weights_path <- file.path(PATHS$results$corepam, "CorePAM_weights.csv")
weights_df   <- strict_csv(weights_path)

stopifnot(all(c("gene", "weight") %in% names(weights_df)))
weights_df    <- weights_df[weights_df$weight != 0, ]
panel_genes   <- weights_df$gene
panel_weights <- setNames(weights_df$weight, weights_df$gene)
message(sprintf("[%s] Painel CorePAM: %d genes com peso != 0", SCRIPT_NAME, length(panel_genes)))

# --------------------------------------------------------------------------
# 2) Carregar expressao pre-Z
# --------------------------------------------------------------------------
expr_path <- file.path(proc_cohort(COHORT), "expression_genelevel_preZ.parquet")
expr_mat   <- strict_parquet(expr_path)

sample_col <- names(expr_mat)[1]
sample_ids  <- expr_mat[[sample_col]]
gene_cols   <- setdiff(names(expr_mat), sample_col)

message(sprintf("[%s] Expressao: %d amostras x %d genes", SCRIPT_NAME, length(sample_ids), length(gene_cols)))

# --------------------------------------------------------------------------
# 3) Z-score intra-coorte por gene
# --------------------------------------------------------------------------
expr_vals <- as.matrix(expr_mat[, gene_cols])
rownames(expr_vals) <- sample_ids

old_warn <- getOption("warn"); options(warn = 0)
gene_means <- colMeans(expr_vals, na.rm = TRUE)
gene_sds   <- apply(expr_vals, 2, sd, na.rm = TRUE)
options(warn = old_warn)

zero_sd <- gene_sds == 0 | is.na(gene_sds)
if (any(zero_sd)) {
  message(sprintf("[%s] AVISO: %d genes com SD=0 (z=0 para esses genes)", SCRIPT_NAME, sum(zero_sd)))
}
gene_sds[zero_sd] <- 1

z_mat <- sweep(sweep(expr_vals, 2, gene_means, "-"), 2, gene_sds, "/")
z_mat[, zero_sd] <- 0

# --------------------------------------------------------------------------
# 4) Calcular score = sum(w_i * z_i) / sum(|w_i|)  apenas genes presentes
# --------------------------------------------------------------------------
genes_present <- intersect(panel_genes, colnames(z_mat))
genes_missing <- setdiff(panel_genes, colnames(z_mat))

n_present    <- length(genes_present)
n_panel      <- length(panel_genes)
frac_present <- n_present / n_panel

message(sprintf("[%s] Genes presentes: %d / %d (%.1f%%)",
                SCRIPT_NAME, n_present, n_panel, frac_present * 100))
if (length(genes_missing) > 0) {
  message(sprintf("[%s] Genes ausentes: %s", SCRIPT_NAME, paste(genes_missing, collapse = ", ")))
}

if (frac_present < FREEZE$min_genes_fraction) {
  stop(sprintf(
    "[%s] COBERTURA INSUFICIENTE: %.1f%% < %.0f%% minimo exigido",
    SCRIPT_NAME, frac_present * 100, FREEZE$min_genes_fraction * 100
  ))
}

w_present      <- panel_weights[genes_present]
denom_sum_absw <- sum(abs(w_present))

z_sub     <- z_mat[, genes_present, drop = FALSE]
score_raw <- as.vector(z_sub %*% w_present) / denom_sum_absw

# --------------------------------------------------------------------------
# 5) score_z = scale(score) para HR por 1 SD
# --------------------------------------------------------------------------
old_warn <- getOption("warn"); options(warn = 0)
score_z_raw <- as.vector(scale(score_raw))
options(warn = old_warn)

# --------------------------------------------------------------------------
# 6) Carregar clinica para verificar direcao do score
# --------------------------------------------------------------------------
clin_path <- file.path(proc_cohort(COHORT), "clinical_harmonized.parquet")
clin_df   <- strict_parquet(clin_path)

score_df <- tibble(
  sample_id      = sample_ids,
  score          = score_raw,
  score_z        = score_z_raw,
  genes_present  = n_present,
  denom_sum_absw = denom_sum_absw
)

clin_score <- inner_join(clin_df, score_df, by = "sample_id")
n_join <- nrow(clin_score)
message(sprintf("[%s] Join clinica x score: %d amostras", SCRIPT_NAME, n_join))
if (n_join == 0) stop(sprintf("[%s] Join resultou em 0 linhas. Verificar chaves.", SCRIPT_NAME))

# --------------------------------------------------------------------------
# 7) Determinar direcao do score (endpoint primario GSE20685 = OS)
# --------------------------------------------------------------------------
endpoint_col <- "os_time"
event_col    <- "os_event"
stopifnot(all(c(endpoint_col, event_col) %in% names(clin_score)))

n_before_drop <- nrow(clin_score)
clin_score    <- clin_score[!is.na(clin_score[[endpoint_col]]) & clin_score[[endpoint_col]] > 0, ]
n_dropped     <- n_before_drop - nrow(clin_score)
if (n_dropped > 0) {
  message(sprintf("[%s] Amostras com time<=0 removidas: %d", SCRIPT_NAME, n_dropped))
}

old_warn <- getOption("warn"); options(warn = 0)
cox_dir <- tryCatch(
  coxph(Surv(clin_score[[endpoint_col]], clin_score[[event_col]]) ~ clin_score$score_z),
  error = function(e) NULL
)
options(warn = old_warn)

if (!is.null(cox_dir)) {
  hr_dir <- exp(coef(cox_dir)[1])
  message(sprintf("[%s] HR score_z (raw): %.4f", SCRIPT_NAME, hr_dir))
  if (hr_dir < 1) {
    clin_score$score   <- -clin_score$score
    clin_score$score_z <- -clin_score$score_z
    score_direction    <- "inverted"
    message(sprintf("[%s] Direcao invertida (HR < 1)", SCRIPT_NAME))
  } else {
    score_direction <- "original"
    message(sprintf("[%s] Direcao original mantida (HR >= 1)", SCRIPT_NAME))
  }
} else {
  score_direction <- "unknown"
  message(sprintf("[%s] AVISO: Cox de direcao falhou; direcao = unknown", SCRIPT_NAME))
}

clin_score$score_direction <- score_direction

# --------------------------------------------------------------------------
# 8) Salvar analysis_ready.parquet
# --------------------------------------------------------------------------
out_dir <- proc_cohort(COHORT)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_path <- file.path(out_dir, "analysis_ready.parquet")
arrow::write_parquet(clin_score, out_path)
h  <- sha256_file(out_path)
sz <- file.info(out_path)$size / 1e6
message(sprintf("[%s] Salvo: %s (%.2f MB | SHA256: %s)", SCRIPT_NAME, out_path, sz, h))

registry_append(
  cohort    = COHORT,
  file_type = "analysis_ready",
  file_path = out_path,
  sha256    = h,
  status    = "ok",
  script    = SCRIPT_NAME,
  size_mb   = sz,
  extra     = list(genes_present = n_present, score_direction = score_direction)
)

# --------------------------------------------------------------------------
# 9) Resumo suplementar
# --------------------------------------------------------------------------
supp_dir <- PATHS$results$supp
dir.create(supp_dir, showWarnings = FALSE, recursive = TRUE)

supp_df <- tibble(
  cohort              = COHORT,
  n_samples           = nrow(clin_score),
  n_panel_genes       = n_panel,
  genes_present       = n_present,
  genes_missing       = n_panel - n_present,
  frac_present        = round(frac_present, 4),
  denom_sum_absw      = round(denom_sum_absw, 6),
  score_direction     = score_direction,
  score_mean          = round(mean(clin_score$score, na.rm = TRUE), 6),
  score_sd            = round(sd(clin_score$score, na.rm = TRUE), 6),
  n_time_leq0_dropped = n_dropped
)

supp_path <- file.path(supp_dir, sprintf("risk_score_summary_%s.csv", COHORT))
readr::write_csv(supp_df, supp_path)
h2 <- sha256_file(supp_path)
registry_append(
  cohort    = COHORT,
  file_type = "risk_score_summary",
  file_path = supp_path,
  sha256    = h2,
  status    = "ok",
  script    = SCRIPT_NAME,
  size_mb   = file.info(supp_path)$size / 1e6
)

message(sprintf("[%s] CONCLUIDO para %s", SCRIPT_NAME, COHORT))
