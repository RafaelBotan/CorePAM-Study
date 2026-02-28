# =============================================================================
# SCRIPT: 00_setup.R
# PURPOSE: Ambiente reprodutível — paths, strict I/O, helpers, parâmetros
#          congelados. Deve ser sourced no início de todos os outros scripts.
# PROJETO: Core-PAM (Memorial v6.1 / Freeze Core-PAM)
# REGRA:   Nunca conter números de genes (ex: PAM29) em paths ou variáveis.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(digest)
  library(arrow)      # leitura/escrita parquet
  library(jsonlite)   # training card JSON
})

# --------------------------------------------------------------------------
# 1) PATHS — single source of truth (ajustar ROOT_REPO uma única vez)
# --------------------------------------------------------------------------
ROOT_REPO <- normalizePath("Y:/Phd-Genomic-claude", mustWork = TRUE)
DATA_LAKE <- file.path(ROOT_REPO, "01_Base_Pura_CorePAM")

PATHS <- list(
  raw            = file.path(DATA_LAKE, "RAW"),
  processed      = file.path(DATA_LAKE, "PROCESSED"),
  docs           = file.path(ROOT_REPO, "01_docs"),
  registry_docs  = file.path(ROOT_REPO, "01_docs", "registry"),
  scripts        = file.path(ROOT_REPO, "scripts"),
  results = list(
    corepam    = file.path(ROOT_REPO, "results", "corepam"),
    corepam_os = file.path(ROOT_REPO, "results", "corepam_os"),
    main       = file.path(ROOT_REPO, "results", "main"),
    supp       = file.path(ROOT_REPO, "results", "supp")
  ),
  figures = list(
    main   = file.path(ROOT_REPO, "figures", "main"),
    supp   = file.path(ROOT_REPO, "figures", "supp"),
    artigo = file.path(ROOT_REPO, "06_plots", "artigo"),
    tese   = file.path(ROOT_REPO, "06_plots", "tese")
  ),
  run_registry = file.path(ROOT_REPO, "registry", "study_registry.csv")
)

# Helpers de path por coorte
raw_cohort    <- function(cohort) file.path(PATHS$raw, cohort)
proc_cohort   <- function(cohort) file.path(PATHS$processed, cohort)

# --------------------------------------------------------------------------
# 2) STRICT I/O — warning = error (obrigatório; Memorial v6.1 §1.4 / §9.1)
# --------------------------------------------------------------------------
options(warn = 2)

# --------------------------------------------------------------------------
# 3) PARÂMETROS CONGELADOS (analysis_freeze.csv — não editar diretamente)
# --------------------------------------------------------------------------
freeze_path <- file.path(PATHS$registry_docs, "analysis_freeze.csv")

if (!file.exists(freeze_path)) {
  stop(
    "FREEZE NAO ENCONTRADO: ", freeze_path,
    "\nCrie 01_docs/registry/analysis_freeze.csv antes de qualquer analise."
  )
}

.freeze_raw <- read_csv(freeze_path, show_col_types = FALSE)
FREEZE <- setNames(
  lapply(.freeze_raw$value, function(x) {
    num <- suppressWarnings(as.numeric(x))
    if (!is.na(num)) num else x
  }),
  .freeze_raw$parameter
)
rm(.freeze_raw)

# Validação mínima
.required <- c("delta_c", "alpha", "k_folds", "seed_folds",
               "min_genes_fraction", "time_unit_divisor")
.missing  <- setdiff(.required, names(FREEZE))
if (length(.missing) > 0) {
  stop("Parametros ausentes em analysis_freeze.csv: ",
       paste(.missing, collapse = ", "))
}
rm(.required, .missing)

# --------------------------------------------------------------------------
# 4) HELPERS
# --------------------------------------------------------------------------

#' SHA-256 de um arquivo — verifica existência e tamanho > 0
sha256_file <- function(path) {
  if (!file.exists(path))       stop("Arquivo nao encontrado: ", path)
  if (file.info(path)$size == 0) stop("Arquivo vazio (0 bytes): ", path)
  digest::digest(path, algo = "sha256", file = TRUE)
}

#' Normalizar ID de paciente/amostra
normalize_id <- function(x) trimws(toupper(as.character(x)))

#' Leitura estrita de CSV
strict_csv <- function(path, ...) {
  if (!file.exists(path)) stop("CSV nao encontrado: ", path)
  readr::read_csv(path, show_col_types = FALSE, ...)
}

#' Leitura estrita de Parquet
strict_parquet <- function(path) {
  if (!file.exists(path)) stop("Parquet nao encontrado: ", path)
  if (file.info(path)$size == 0) stop("Parquet vazio (0 bytes): ", path)
  arrow::read_parquet(path)
}

#' Leitura estrita de RDS
strict_rds <- function(path) {
  if (!file.exists(path)) stop("RDS nao encontrado: ", path)
  if (file.info(path)$size == 0) stop("RDS vazio (0 bytes): ", path)
  readRDS(path)
}

#' Append ao registry (append-only; cria header se arquivo ainda não existe)
registry_append <- function(cohort, file_type, file_path, sha256,
                            status, script, size_mb = NA_real_,
                            extra = list()) {
  path <- PATHS$run_registry
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)

  row <- tibble(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    script    = script,
    cohort    = cohort,
    file_type = file_type,
    path      = file_path,
    sha256    = sha256,
    status    = status,
    size_mb   = round(size_mb, 3)
  )
  if (length(extra) > 0) row <- bind_cols(row, as_tibble(extra))

  write_csv(row, path,
            append   = file.exists(path),
            col_names = !file.exists(path))
  invisible(row)
}

# --------------------------------------------------------------------------
# 5) VALIDAÇÃO DE ESTRUTURA DE PASTAS (cria se ausente)
# --------------------------------------------------------------------------
.required_dirs <- c(
  PATHS$raw, PATHS$processed,
  PATHS$docs, PATHS$registry_docs,
  PATHS$scripts,
  PATHS$results$corepam, PATHS$results$corepam_os,
  PATHS$results$main,    PATHS$results$supp,
  PATHS$figures$main,    PATHS$figures$supp,
  PATHS$figures$artigo,  PATHS$figures$tese,
  dirname(PATHS$run_registry)
)
for (.d in .required_dirs) dir.create(.d, showWarnings = FALSE, recursive = TRUE)
rm(.required_dirs, .d)

# --------------------------------------------------------------------------
message(sprintf(
  "[00_setup] OK | ROOT: %s | delta_c=%.3f | alpha=%.1f | K=%d | seed=%d",
  ROOT_REPO, FREEZE$delta_c, FREEZE$alpha, FREEZE$k_folds, FREEZE$seed_folds
))
