# =============================================================================
# SCRIPT: 02_harmonize_clinical_SCANB.R
# PURPOSE: Harmonização clínica do SCAN-B (GSE96058) e split determinístico
#          treino (SCAN-B) vs validação (GSE96058).
# PROJETO: Core-PAM (Memorial v6.1)
#
# INPUT:
#   01_Base_Pura_CorePAM/RAW/GSE96058/GSE96058_clinical_raw.rds  (pData GEO)
#
# OUTPUTS:
#   01_Base_Pura_CorePAM/PROCESSED/SCANB/clinical_FINAL.parquet
#   01_Base_Pura_CorePAM/PROCESSED/GSE96058/clinical_FINAL.parquet
#   01_docs/endpoint_mapping_templates/endpoint_mapping_SCANB.csv
#   01_docs/endpoint_mapping_templates/endpoint_mapping_GSE96058.csv
#
# REGRAS (Memorial v6.1):
#   - Tempo em meses (dias / 30.4375).
#   - time <= 0: DROP (reportar contagem).
#   - event OS: 1 = óbito qualquer causa; 0 = censura.
#   - Sem imputação: reportar N efetivo por variável.
#   - Leakage-proof: IDs treino e validação são conjuntos disjuntos.
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "02_harmonize_clinical_SCANB.R"

# =============================================================================
# 1) LEITURA ESTRITA
# =============================================================================
raw_path <- file.path(raw_cohort("GSE96058"), "GSE96058_clinical_raw.rds")
pheno    <- strict_rds(raw_path)

message(sprintf("[02_SCANB] pData carregado: %d amostras x %d colunas",
                nrow(pheno), ncol(pheno)))

# Inspecionar colunas disponíveis (log para auditoria)
message("[02_SCANB] Colunas disponiveis:")
print(names(pheno))

# =============================================================================
# 2) IDENTIFICAR COLUNA DE SPLIT TREINO / VALIDAÇÃO
#    GEO GSE96058 contém coluna indicando "training" vs "validation".
#    Nome exato pode variar — identificar automaticamente e registrar.
# =============================================================================
.split_candidates <- grep(
  "train|valid|group|cohort|set",
  tolower(names(pheno)), value = TRUE
)
message("[02_SCANB] Candidatos a coluna de split: ",
        paste(names(pheno)[tolower(names(pheno)) %in% .split_candidates],
              collapse = ", "))

# Ajustar SPLIT_COL se o nome automático não for correto
SPLIT_COL <- names(pheno)[
  tolower(names(pheno)) %in% .split_candidates
][1]

if (is.na(SPLIT_COL) || length(SPLIT_COL) == 0) {
  stop(
    "[02_SCANB] Coluna de split nao detectada automaticamente.\n",
    "Inspecione names(pheno) e defina SPLIT_COL manualmente."
  )
}
message("[02_SCANB] Usando coluna de split: ", SPLIT_COL)
message("[02_SCANB] Valores unicos: ", paste(unique(pheno[[SPLIT_COL]]), collapse = " | "))

# =============================================================================
# 3) MAPEAMENTO DE COLUNAS (ajustar nomes conforme pData real)
#    Estes são os nomes ESPERADOS no pData GSE96058.
#    Se diferirem, ajuste o mapa abaixo.
# =============================================================================
COL_MAP <- list(
  sample_id   = "geo_accession",          # GSM ID
  patient_id  = "title",                  # ID do paciente (barcode SCAN-B)
  os_time_raw = "overall survival (months):ch1",   # tempo em meses (verificar unidade)
  os_event    = "overall survival event:ch1",       # 0/1 ou dead/alive
  age         = "age at diagnosis:ch1",
  er_status   = "er status:ch1",
  stage       = "tumor stage:ch1"
)

# =============================================================================
# 4) FUNÇÃO DE HARMONIZAÇÃO (aplicada a treino e validação)
# =============================================================================
harmonize_scanb <- function(df, cohort_label) {
  n_raw <- nrow(df)

  out <- tibble(
    sample_id  = normalize_id(df[[COL_MAP$sample_id]]),
    patient_id = normalize_id(df[[COL_MAP$patient_id]])
  )

  # Tempo de sobrevida — verificar se já está em meses ou está em dias
  raw_time <- suppressWarnings(as.numeric(df[[COL_MAP$os_time_raw]]))

  # Heurística: se mediana > 200 assume dias → converter
  if (median(raw_time, na.rm = TRUE) > 200) {
    message(sprintf(
      "[02_%s] Tempo em dias detectado (mediana=%.1f). Convertendo para meses.",
      cohort_label, median(raw_time, na.rm = TRUE)
    ))
    out$os_time_months <- raw_time / FREEZE$time_unit_divisor
  } else {
    message(sprintf(
      "[02_%s] Tempo ja em meses (mediana=%.1f).",
      cohort_label, median(raw_time, na.rm = TRUE)
    ))
    out$os_time_months <- raw_time
  }

  # Evento OS: padronizar para 0/1
  raw_event <- tolower(trimws(df[[COL_MAP$os_event]]))
  out$os_event <- dplyr::case_when(
    raw_event %in% c("1", "dead", "deceased", "yes", "true")  ~ 1L,
    raw_event %in% c("0", "alive", "living",  "no",  "false") ~ 0L,
    TRUE ~ NA_integer_
  )

  # Covariáveis CORE-A
  out$age      <- suppressWarnings(as.numeric(df[[COL_MAP$age]]))
  out$er_status <- dplyr::case_when(
    tolower(trimws(df[[COL_MAP$er_status]])) %in% c("positive", "pos", "1", "+") ~ "Positive",
    tolower(trimws(df[[COL_MAP$er_status]])) %in% c("negative", "neg", "0", "-") ~ "Negative",
    TRUE ~ NA_character_
  )
  out$stage <- trimws(df[[COL_MAP$stage]])

  # DROP: time <= 0
  n_le0 <- sum(out$os_time_months <= 0, na.rm = TRUE)
  out   <- out |> filter(os_time_months > 0)
  n_na  <- sum(is.na(out$os_time_months))
  out   <- out |> filter(!is.na(os_time_months))

  message(sprintf(
    "[02_%s] N raw=%d | Removidos time<=0: %d | Removidos time=NA: %d | N final=%d",
    cohort_label, n_raw, n_le0, n_na, nrow(out)
  ))
  message(sprintf(
    "[02_%s] Eventos OS: %d (%.1f%%) | ER+ %d | ER- %d | ER NA %d",
    cohort_label,
    sum(out$os_event, na.rm = TRUE),
    100 * mean(out$os_event, na.rm = TRUE),
    sum(out$er_status == "Positive", na.rm = TRUE),
    sum(out$er_status == "Negative", na.rm = TRUE),
    sum(is.na(out$er_status))
  ))

  out
}

# =============================================================================
# 5) SPLIT E HARMONIZAÇÃO
# =============================================================================
train_mask <- grepl("train", tolower(pheno[[SPLIT_COL]]))
valid_mask <- !train_mask

pheno_train <- pheno[train_mask, ]
pheno_valid <- pheno[valid_mask, ]

message(sprintf("[02_SCANB] Split: treino=%d | validacao=%d",
                nrow(pheno_train), nrow(pheno_valid)))

# Verificar disjunção (leakage-proof)
ids_train <- normalize_id(pheno_train[[COL_MAP$sample_id]])
ids_valid <- normalize_id(pheno_valid[[COL_MAP$sample_id]])
overlap   <- intersect(ids_train, ids_valid)
if (length(overlap) > 0) {
  stop(sprintf("[02_SCANB] LEAKAGE: %d IDs sobrepostos entre treino e validacao: %s",
               length(overlap), paste(overlap[1:min(5, length(overlap))], collapse = ",")))
}
message("[02_SCANB] Leakage check: PASS (overlap = 0)")

clin_train <- harmonize_scanb(pheno_train, "SCANB")
clin_valid <- harmonize_scanb(pheno_valid, "GSE96058")

# =============================================================================
# 6) SALVAR PARQUETS FINAIS
# =============================================================================
.save_clinical <- function(df, cohort) {
  dest_dir <- proc_cohort(cohort)
  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
  out_path <- file.path(dest_dir, "clinical_FINAL.parquet")
  arrow::write_parquet(df, out_path)

  h    <- sha256_file(out_path)
  size <- file.info(out_path)$size / 1024^2
  registry_append(cohort, "Clinical_FINAL", out_path, h,
                  "INTEGRO", SCRIPT_NAME, size)
  message(sprintf("[02_%s] Salvo: %s | %d amostras | %.2f MB | SHA256: %s",
                  cohort, out_path, nrow(df), size, h))
}

.save_clinical(clin_train, "SCANB")
.save_clinical(clin_valid, "GSE96058")

# =============================================================================
# 7) ENDPOINT MAPPING (dicionário de auditoria)
# =============================================================================
.save_endpoint_map <- function(cohort, n_samples, n_events,
                               time_col_origin, event_col_origin,
                               time_unit_origin) {
  out <- tibble(
    cohort           = cohort,
    endpoint         = "OS",
    col_time_origin  = time_col_origin,
    col_event_origin = event_col_origin,
    unit_origin      = time_unit_origin,
    unit_final       = "months",
    conversion       = ifelse(time_unit_origin == "days",
                              "/ 30.4375", "none"),
    event_1_means    = "death_any_cause",
    n_samples        = n_samples,
    n_events         = n_events,
    script           = SCRIPT_NAME,
    timestamp        = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  dest <- file.path(PATHS$docs, "endpoint_mapping_templates",
                    sprintf("endpoint_mapping_%s.csv", cohort))
  write_csv(out, dest)
  message(sprintf("[02_%s] Endpoint map salvo: %s", cohort, dest))
}

.save_endpoint_map("SCANB", nrow(clin_train),
                   sum(clin_train$os_event, na.rm = TRUE),
                   COL_MAP$os_time_raw, COL_MAP$os_event, "months_or_days")
.save_endpoint_map("GSE96058", nrow(clin_valid),
                   sum(clin_valid$os_event, na.rm = TRUE),
                   COL_MAP$os_time_raw, COL_MAP$os_event, "months_or_days")

message("\n[02_SCANB] Concluido: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message("Proximo passo: scripts/03_expression_preprocess_SCANB.R")
