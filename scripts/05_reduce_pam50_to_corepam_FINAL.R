# =============================================================================
# SCRIPT: 05_reduce_pam50_to_corepam_FINAL.R
# PURPOSE: Derivação do Core-PAM — menor subconjunto PAM50 não-inferior ao
#          modelo PAM50 completo por C-index OOF (Harrell) no SCAN-B.
# PROJETO: Core-PAM (Memorial v6.1 §5)
#
# ALGORITMO (congelado):
#   1. Preparar X (50 genes PAM50, Z-score intra-SCANB) + y (Surv OS)
#   2. Folds K=10 estratificados fixos (seed congelado)
#   3. cv.glmnet (alpha=0.5, keep=TRUE) → OOF linear predictors por lambda
#   4. Para cada ponto do path (lambda, df):
#        Cadj = max(Craw, 1-Craw)  — orientation-free
#   5. Cmax = max(Cadj); threshold = Cmax - delta_c
#   6. Selecionar menor df com Cadj >= threshold
#      Desempate: maior lambda (mais parcimonioso) dentro do mesmo df
#   7. Refit no full data: glmnet.fit no lambda escolhido
#   8. Freeze: CorePAM_weights.csv, CorePAM_model.rds,
#              CorePAM_training_card.json, pareto_df_cindex_oof.csv
#
# INPUTS:
#   01_Base_Pura_CorePAM/PROCESSED/SCANB/expression_genelevel_preZ.parquet
#   01_Base_Pura_CorePAM/PROCESSED/SCANB/clinical_FINAL.parquet
#
# OUTPUTS (todos em results/corepam/):
#   pareto_df_cindex_oof.csv
#   CorePAM_weights.csv
#   CorePAM_model.rds
#   CorePAM_training_card.json
#   artifact_hashes.csv
# =============================================================================

source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(glmnet)
  library(survival)
  library(jsonlite)
})

SCRIPT_NAME <- "05_reduce_pam50_to_corepam_FINAL.R"

PAM50_GENES <- c(
  "ACTR3B","ANLN","BAG1","BCL2","BIRC5","BLVRA","CCNB1","CCNE1","CDC20",
  "CDC6","CDH3","CENPF","CEP55","CXXC5","EGFR","ERBB2","ESR1","EXO1",
  "FGFR4","FOXA1","FOXC1","GPR160","GRB7","KIF2C","KRT14","KRT17","KRT5",
  "MAPT","MDM2","MELK","MIA","MKI67","MLPH","MMP11","MYBL2","MYC","NAT1",
  "NDC80","NUF2","ORC6","PGR","PHGDH","PTTG1","RRM2","SFRP1","SLC39A6",
  "TMEM45B","TYMS","UBE2C","UBE2T"
)

# =============================================================================
# 1) CARREGAR DADOS SCANB
# =============================================================================
message("[05] Carregando SCANB expressão e clínica...")

expr <- strict_parquet(file.path(proc_cohort("SCANB"),
                                 "expression_genelevel_preZ.parquet"))
clin <- strict_parquet(file.path(proc_cohort("SCANB"),
                                 "clinical_FINAL.parquet"))

message(sprintf("[05] Expressão: %d genes x %d amostras", nrow(expr), ncol(expr) - 1))
message(sprintf("[05] Clínica: %d amostras | OS eventos: %d (%.1f%%)",
                nrow(clin),
                sum(clin$os_event, na.rm = TRUE),
                100 * mean(clin$os_event, na.rm = TRUE)))

# =============================================================================
# 2) ALINHAR AMOSTRAS (expressão × clínica)
# =============================================================================
expr_mat  <- as.matrix(expr[, -1])
rownames(expr_mat) <- expr$gene
expr_mat  <- t(expr_mat)           # amostras × genes
rownames(expr_mat) <- normalize_id(rownames(expr_mat))

clin_ids  <- normalize_id(clin$sample_id)
common_ids <- intersect(rownames(expr_mat), clin_ids)

if (length(common_ids) == 0) {
  stop("[05] Nenhuma amostra em comum entre expressão e clínica SCANB.")
}
message(sprintf("[05] Amostras em comum: %d", length(common_ids)))

expr_mat <- expr_mat[common_ids, , drop = FALSE]
clin_sub <- clin[match(common_ids, clin_ids), ]

# =============================================================================
# 3) FILTRAR GENES PAM50 E APLICAR Z-SCORE INTRA-SCANB
#    Z-score aqui serve apenas para a derivação (coeficientes em unidade SD).
#    O Z-score de validação será intra-coorte em 06_zscore_and_score_<COHORT>.R
# =============================================================================
pam50_available <- intersect(PAM50_GENES, colnames(expr_mat))
pam50_missing   <- setdiff(PAM50_GENES, colnames(expr_mat))

if (length(pam50_missing) > 0) {
  message(sprintf("[05] AVISO: %d genes PAM50 ausentes no SCANB: %s",
                  length(pam50_missing), paste(pam50_missing, collapse = ",")))
}
if (length(pam50_available) < 40) {
  stop(sprintf("[05] Apenas %d genes PAM50 disponíveis — mínimo 40 exigido para derivação.",
               length(pam50_available)))
}
message(sprintf("[05] Genes PAM50 disponíveis para derivação: %d/50",
                length(pam50_available)))

X_raw <- expr_mat[, pam50_available, drop = FALSE]

# Z-score por gene intra-SCANB (para derivação)
X <- scale(X_raw)
# Remover genes sem variância (sd=0 → coluna NA após scale)
zero_var <- apply(X, 2, function(col) all(is.na(col)))
if (any(zero_var)) {
  message(sprintf("[05] Removendo %d genes sem variância.", sum(zero_var)))
  X <- X[, !zero_var, drop = FALSE]
}
message(sprintf("[05] Matriz X final: %d amostras x %d genes (Z-scored)",
                nrow(X), ncol(X)))

# Resposta de sobrevida
y <- survival::Surv(time  = clin_sub$os_time_months,
                    event = clin_sub$os_event)

# =============================================================================
# 4) FOLDS ESTRATIFICADOS FIXOS (K=10, seed congelado)
#    Estratificação por evento para balancear eventos/censuras nos folds.
# =============================================================================
set.seed(FREEZE$seed_folds)

# Estratificação manual por status de evento
event_idx   <- which(clin_sub$os_event == 1L)
censor_idx  <- which(clin_sub$os_event == 0L)

foldid <- integer(nrow(X))
foldid[event_idx]  <- sample(rep(1:FREEZE$k_folds,
                                  length.out = length(event_idx)))
foldid[censor_idx] <- sample(rep(1:FREEZE$k_folds,
                                  length.out = length(censor_idx)))

message(sprintf("[05] Folds criados: K=%d | seed=%d | eventos por fold: ~%d",
                FREEZE$k_folds, FREEZE$seed_folds,
                round(sum(clin_sub$os_event) / FREEZE$k_folds)))

# =============================================================================
# 5) CV.GLMNET COM keep=TRUE (armazena OOF linear predictors)
# =============================================================================
message("[05] Rodando cv.glmnet (alpha=", FREEZE$alpha,
        ", K=", FREEZE$k_folds, ")...")

old_warn <- getOption("warn"); options(warn = 0)
cv_fit <- glmnet::cv.glmnet(
  x       = X,
  y       = y,
  family  = "cox",
  alpha   = FREEZE$alpha,
  foldid  = foldid,
  keep    = TRUE,   # armazena fit.preval (OOF linear predictors)
  grouped = TRUE
)
options(warn = old_warn)

message(sprintf("[05] Lambda path: %d pontos | lambda.min=%.6f | lambda.1se=%.6f",
                length(cv_fit$lambda), cv_fit$lambda.min, cv_fit$lambda.1se))

# =============================================================================
# 6) CALCULAR C-INDEX OOF PARA CADA PONTO DO PATH
#    fit.preval: matriz n_amostras × n_lambdas com OOF linear predictors
# =============================================================================
message("[05] Calculando C-index OOF por lambda...")

oof_lp  <- cv_fit$fit.preval   # n × n_lambda
lambdas <- cv_fit$lambda
dfs     <- cv_fit$glmnet.fit$df  # número de coeficientes não-zero por lambda

pareto_list <- vector("list", length(lambdas))

for (j in seq_along(lambdas)) {
  lp <- oof_lp[, j]

  # Amostras com predição válida (algumas podem ser NA nos primeiros lambdas)
  valid_idx <- !is.na(lp)
  if (sum(valid_idx) < 20) {
    pareto_list[[j]] <- tibble(
      lambda = lambdas[j], df = dfs[j],
      craw = NA_real_, cadj = NA_real_, n_valid = sum(valid_idx)
    )
    next
  }

  old_warn <- getOption("warn"); options(warn = 0)
  craw <- tryCatch({
    conc <- survival::concordance(y[valid_idx] ~ lp[valid_idx])
    as.numeric(conc$concordance)
  }, error = function(e) NA_real_)
  options(warn = old_warn)

  cadj <- if (!is.na(craw)) max(craw, 1 - craw) else NA_real_

  pareto_list[[j]] <- tibble(
    lambda  = lambdas[j],
    df      = dfs[j],
    craw    = craw,
    cadj    = cadj,
    n_valid = sum(valid_idx)
  )
}

df_pareto <- bind_rows(pareto_list) |>
  filter(!is.na(cadj)) |>
  # Para cada df: manter apenas o maior lambda (mais parcimonioso)
  group_by(df) |>
  slice_max(lambda, n = 1, with_ties = FALSE) |>
  ungroup() |>
  arrange(df)

message(sprintf("[05] Pareto calculado: %d pontos únicos de df", nrow(df_pareto)))

# =============================================================================
# 7) SELEÇÃO POR NÃO-INFERIORIDADE
# =============================================================================
cmax      <- max(df_pareto$cadj, na.rm = TRUE)
threshold <- cmax - FREEZE$delta_c

message(sprintf("[05] Cmax=%.4f | threshold=%.4f (delta_c=%.3f)",
                cmax, threshold, FREEZE$delta_c))

# Menor df com Cadj >= threshold; desempate: maior lambda
selected <- df_pareto |>
  filter(cadj >= threshold) |>
  arrange(df, desc(lambda)) |>
  slice(1)

if (nrow(selected) == 0) {
  stop("[05] Nenhum modelo satisfaz o critério de não-inferioridade. ",
       "Verifique delta_c em analysis_freeze.csv.")
}

message(sprintf(
  "[05] Core-PAM SELECIONADO: df=%d genes | Cadj=%.4f | lambda=%.6f",
  selected$df, selected$cadj, selected$lambda
))

# =============================================================================
# 8) REFIT NO FULL DATA (glmnet.fit no lambda escolhido — Memorial §5.4)
#    Usar glmnet.fit do cv_fit (full-data path), não refit single-lambda
# =============================================================================
message("[05] Extraindo coeficientes do full-data path no lambda escolhido...")

old_warn <- getOption("warn"); options(warn = 0)
coefs_raw <- coef(cv_fit$glmnet.fit, s = selected$lambda)
options(warn = old_warn)

# Genes com coeficiente não-zero
nonzero_idx <- which(coefs_raw != 0)
corepam_genes  <- rownames(coefs_raw)[nonzero_idx]
corepam_weights <- as.numeric(coefs_raw[nonzero_idx])
n_genes <- length(corepam_genes)

if (n_genes != selected$df) {
  message(sprintf("[05] AVISO: df do path=%d mas coefs não-zero=%d (discrepância normal por arredondamento).",
                  selected$df, n_genes))
}

message(sprintf("[05] Core-PAM: %d genes | Pesos (primeiros 10):", n_genes))
weights_df <- tibble(gene = corepam_genes, weight = corepam_weights) |>
  arrange(desc(abs(weight)))
print(head(weights_df, 10))

# =============================================================================
# 9) DIREÇÃO DO SCORE (Memorial §5.5)
#    Se HR < 1 na coorte de treino → inverter sinal dos pesos
# =============================================================================
score_train <- as.numeric(X[, corepam_genes, drop = FALSE] %*% corepam_weights)
score_train_z <- scale(score_train)[, 1]

old_warn <- getOption("warn"); options(warn = 0)
cox_dir <- survival::coxph(y ~ score_train_z)
options(warn = old_warn)

hr_train <- exp(coef(cox_dir))
message(sprintf("[05] HR treino (por 1 SD): %.3f", hr_train))

if (hr_train < 1) {
  message("[05] HR < 1 → invertendo sinal dos pesos (score_direction = -1)")
  corepam_weights <- -corepam_weights
  score_direction <- -1L
} else {
  message("[05] HR >= 1 → direção original (score_direction = +1)")
  score_direction <- 1L
}

weights_final <- tibble(gene = corepam_genes, weight = corepam_weights) |>
  arrange(desc(abs(weight)))

# =============================================================================
# 10) SALVAR ARTEFATOS (results/corepam/)
# =============================================================================
dir.create(PATHS$results$corepam, showWarnings = FALSE, recursive = TRUE)

artifact_hashes <- list()

# --- 10a) pareto_df_cindex_oof.csv ---
pareto_path <- file.path(PATHS$results$corepam, "pareto_df_cindex_oof.csv")
write_csv(df_pareto, pareto_path)
h <- sha256_file(pareto_path)
registry_append("SCANB", "Pareto_OOF", pareto_path, h, "INTEGRO", SCRIPT_NAME,
                file.info(pareto_path)$size / 1024^2)
artifact_hashes$pareto_df_cindex_oof <- h
message("[05] Salvo: ", pareto_path)

# --- 10b) CorePAM_weights.csv ---
weights_path <- file.path(PATHS$results$corepam, "CorePAM_weights.csv")
write_csv(weights_final, weights_path)
h <- sha256_file(weights_path)
registry_append("SCANB", "CorePAM_Weights", weights_path, h, "INTEGRO", SCRIPT_NAME,
                file.info(weights_path)$size / 1024^2)
artifact_hashes$CorePAM_weights <- h
message("[05] Salvo: ", weights_path)

# --- 10c) CorePAM_model.rds (glmnet.fit completo) ---
model_path <- file.path(PATHS$results$corepam, "CorePAM_model.rds")
saveRDS(cv_fit$glmnet.fit, model_path)
h <- sha256_file(model_path)
registry_append("SCANB", "CorePAM_Model", model_path, h, "INTEGRO", SCRIPT_NAME,
                file.info(model_path)$size / 1024^2)
artifact_hashes$CorePAM_model <- h
message("[05] Salvo: ", model_path)

# --- 10d) CorePAM_training_card.json ---
training_card <- list(
  project          = "Core-PAM",
  memorial_version = "v6.1",
  script           = SCRIPT_NAME,
  timestamp        = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
  cohort_train     = "SCANB",
  n_samples_train  = nrow(X),
  n_events_train   = sum(clin_sub$os_event, na.rm = TRUE),
  freeze = list(
    alpha              = FREEZE$alpha,
    k_folds            = FREEZE$k_folds,
    seed_folds         = FREEZE$seed_folds,
    delta_c            = FREEZE$delta_c,
    min_genes_fraction = FREEZE$min_genes_fraction
  ),
  derivation = list(
    pam50_genes_available = length(pam50_available),
    pam50_genes_missing   = pam50_missing,
    lambda_chosen  = selected$lambda,
    df_chosen      = selected$df,
    cadj_chosen    = selected$cadj,
    cmax           = cmax,
    threshold      = threshold,
    score_direction = score_direction
  ),
  corepam = list(
    n_genes  = n_genes,
    genes    = corepam_genes,
    weights  = corepam_weights
  ),
  artifact_hashes = artifact_hashes
)

card_path <- file.path(PATHS$results$corepam, "CorePAM_training_card.json")
jsonlite::write_json(training_card, card_path, pretty = TRUE, auto_unbox = TRUE)
h <- sha256_file(card_path)
registry_append("SCANB", "CorePAM_Training_Card", card_path, h,
                "INTEGRO", SCRIPT_NAME,
                file.info(card_path)$size / 1024^2)
artifact_hashes$CorePAM_training_card <- h
message("[05] Salvo: ", card_path)

# --- 10e) artifact_hashes.csv ---
hashes_df   <- tibble(artifact = names(artifact_hashes),
                      sha256   = unlist(artifact_hashes))
hashes_path <- file.path(PATHS$results$corepam, "artifact_hashes.csv")
write_csv(hashes_df, hashes_path)
message("[05] Salvo: ", hashes_path)

# =============================================================================
# 11) RELATÓRIO FINAL (impresso no console)
# =============================================================================
message("\n", strrep("=", 60))
message("[05] CORE-PAM — RESULTADO FINAL (FREEZE)")
message(strrep("=", 60))
message(sprintf("  Genes no Core-PAM : %d", n_genes))
message(sprintf("  Cadj OOF          : %.4f", selected$cadj))
message(sprintf("  Cmax PAM50        : %.4f", cmax))
message(sprintf("  Margem usada      : %.3f (delta_c congelado)", FREEZE$delta_c))
message(sprintf("  Lambda escolhido  : %.6f", selected$lambda))
message(sprintf("  Score direction   : %+d", score_direction))
message(sprintf("  HR treino (1 SD)  : %.3f", exp(coef(cox_dir)) * score_direction^2))
message("\n  Genes e pesos:")
print(weights_final, n = Inf)
message(strrep("=", 60))
message("[05] Artefatos em: ", PATHS$results$corepam)
message("[05] Concluido: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message("Proximo passo: scripts/06_zscore_and_score_<COHORT>.R")
