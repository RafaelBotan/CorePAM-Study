# =============================================================================
# SCRIPT: 02_harmonize_clinical_TCGA_BRCA.R
# PURPOSE: Harmonização clínica do TCGA-BRCA (GDC harmonized).
# PROJETO: Core-PAM (Memorial v6.1)
#
# INPUTS:
#   01_Base_Pura_CorePAM/RAW/TCGA_BRCA/TCGA_BRCA_SE_counts_raw.rds  (SE)
#   01_Base_Pura_CorePAM/RAW/TCGA_BRCA/TCGA_BRCA_clinical_raw.rds   (GDCquery_clinic)
#
# OUTPUTS:
#   01_Base_Pura_CorePAM/PROCESSED/TCGA_BRCA/clinical_FINAL.parquet
#   01_docs/endpoint_mapping_templates/endpoint_mapping_TCGA_BRCA.csv
#
# REGRAS (Memorial v6.1):
#   - Endpoint: OS (GDC: days_to_death / days_to_last_follow_up).
#   - Tempo em meses (dias / 30.4375).
#   - time <= 0: DROP.
#   - Sensibilidade pré-declarada: horizonte 24m (executada em 07_*.R).
#   - Apenas Primary Tumor (excluir metastáticos, normais).
#   - Barcode TCGA: usar 12 caracteres (TCGA-XX-XXXX) como patient_id.
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "02_harmonize_clinical_TCGA_BRCA.R"

# =============================================================================
# 1) LEITURA ESTRITA
# =============================================================================
clin_path <- file.path(raw_cohort("TCGA_BRCA"), "TCGA_BRCA_clinical_raw.rds")
se_path   <- file.path(raw_cohort("TCGA_BRCA"), "TCGA_BRCA_SE_counts_raw.rds")

clin_raw <- strict_rds(clin_path)
se_obj   <- strict_rds(se_path)

message(sprintf("[02_TCGA] Clinica GDC: %d pacientes x %d colunas",
                nrow(clin_raw), ncol(clin_raw)))
message("[02_TCGA] Colunas clinicas: ", paste(names(clin_raw), collapse = ", "))

# Barcode das amostras no SE (Primary Tumor = tipo 01)
se_barcodes <- colnames(se_obj)
primary_barcodes <- se_barcodes[substr(se_barcodes, 14, 15) == "01"]
patient_ids_se   <- substr(primary_barcodes, 1, 12)

message(sprintf("[02_TCGA] SE: %d barcodes totais | %d Primary Tumor",
                length(se_barcodes), length(primary_barcodes)))

# =============================================================================
# 2) MAPEAMENTO DE COLUNAS GDC (nomes padrão GDCquery_clinic)
# =============================================================================
# GDC harmonized clinical usa nomes snake_case
PAT_COL           <- "submitter_id"          # TCGA-XX-XXXX
DAYS_DEATH_COL    <- "days_to_death"
DAYS_FOLLOW_COL   <- "days_to_last_follow_up"
VITAL_COL         <- "vital_status"          # "Dead" / "Alive"
AGE_COL           <- "age_at_index"
GENDER_COL        <- "gender"

# ER status vem do colData do SE (não da clínica GDC diretamente)
# Extrair do SE se disponível
if ("paper_BRCA_Subtype_PAM50" %in% names(SummarizedExperiment::colData(se_obj))) {
  se_coldata <- as.data.frame(SummarizedExperiment::colData(se_obj))
  se_coldata$patient_id <- substr(rownames(se_coldata), 1, 12)
  se_coldata <- se_coldata[substr(rownames(se_coldata), 14, 15) == "01", ]
  message(sprintf("[02_TCGA] SE colData com %d colunas encontradas.",
                  ncol(se_coldata)))
} else {
  se_coldata <- NULL
  message("[02_TCGA] AVISO: colData do SE sem colunas de subtipo — ER via clinica GDC.")
}

# =============================================================================
# 3) HARMONIZAÇÃO
# =============================================================================
n_raw <- nrow(clin_raw)

out <- tibble(
  patient_id = normalize_id(clin_raw[[PAT_COL]])
)

# Tempo OS: days_to_death (óbito) ou days_to_last_follow_up (censura)
days_death  <- suppressWarnings(as.numeric(clin_raw[[DAYS_DEATH_COL]]))
days_follow <- suppressWarnings(as.numeric(clin_raw[[DAYS_FOLLOW_COL]]))

out$os_event <- dplyr::case_when(
  tolower(trimws(clin_raw[[VITAL_COL]])) %in% c("dead", "deceased", "1") ~ 1L,
  tolower(trimws(clin_raw[[VITAL_COL]])) %in% c("alive", "living",   "0") ~ 0L,
  TRUE ~ NA_integer_
)

# Para óbitos: usar days_to_death; para vivos: days_to_last_follow_up
out$os_time_days <- dplyr::case_when(
  out$os_event == 1L & !is.na(days_death)  ~ days_death,
  out$os_event == 0L & !is.na(days_follow) ~ days_follow,
  !is.na(days_death)                        ~ days_death,
  TRUE                                       ~ days_follow
)
out$os_time_months <- out$os_time_days / FREEZE$time_unit_divisor

# Covariáveis
out$age    <- suppressWarnings(as.numeric(clin_raw[[AGE_COL]]))
out$gender <- tolower(trimws(clin_raw[[GENDER_COL]]))

# ER status via SE colData (se disponível) ou marcar como NA
if (!is.null(se_coldata) && "ER_Status_By_IHC" %in% names(se_coldata)) {
  er_lkp <- setNames(se_coldata$ER_Status_By_IHC,
                     normalize_id(se_coldata$patient_id))
  out$er_status <- dplyr::case_when(
    tolower(er_lkp[out$patient_id]) %in% c("positive", "pos") ~ "Positive",
    tolower(er_lkp[out$patient_id]) %in% c("negative", "neg") ~ "Negative",
    TRUE ~ NA_character_
  )
} else {
  out$er_status <- NA_character_
  message("[02_TCGA] AVISO: ER status nao disponivel — marcar como NA.")
}

# Filtrar apenas amostras com Primary Tumor no SE
out <- out |> filter(patient_id %in% normalize_id(patient_ids_se))

# DROP: time <= 0
n_le0 <- sum(out$os_time_months <= 0, na.rm = TRUE)
n_na  <- sum(is.na(out$os_time_months))
out   <- out |> filter(os_time_months > 0, !is.na(os_time_months))

message(sprintf(
  "[02_TCGA] N raw=%d | Primary Tumor filtrado | Removidos time<=0: %d | NA: %d | N final=%d",
  n_raw, n_le0, n_na, nrow(out)
))
message(sprintf(
  "[02_TCGA] OS eventos: %d (%.1f%%) | Follow-up mediana: %.1f meses",
  sum(out$os_event, na.rm = TRUE),
  100 * mean(out$os_event, na.rm = TRUE),
  median(out$os_time_months, na.rm = TRUE)
))
message("[02_TCGA] NOTA: sensibilidade horizonte 24m pre-declarada (executada em 07_*.R)")

# =============================================================================
# 4) SALVAR
# =============================================================================
dest_dir <- proc_cohort("TCGA_BRCA")
dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(dest_dir, "clinical_FINAL.parquet")

arrow::write_parquet(out, out_path)
h    <- sha256_file(out_path)
size <- file.info(out_path)$size / 1024^2
registry_append("TCGA_BRCA", "Clinical_FINAL", out_path, h,
                "INTEGRO", SCRIPT_NAME, size)
message(sprintf("[02_TCGA] Salvo: %s | %d amostras | %.2f MB | SHA256: %s",
                out_path, nrow(out), size, h))

# =============================================================================
# 5) ENDPOINT MAPPING
# =============================================================================
ep_map <- tibble(
  cohort              = "TCGA_BRCA",
  endpoint_primary    = "OS",
  col_days_death      = DAYS_DEATH_COL,
  col_days_follow_up  = DAYS_FOLLOW_COL,
  col_vital_status    = VITAL_COL,
  unit_origin         = "days",
  unit_final          = "months",
  conversion          = paste0("/ ", FREEZE$time_unit_divisor),
  event_1_means       = "death_any_cause",
  sensitivity_note    = "horizonte_24m_executado_em_07_survival_TCGA_BRCA.R",
  n_samples           = nrow(out),
  n_events            = sum(out$os_event, na.rm = TRUE),
  script              = SCRIPT_NAME,
  timestamp           = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)
ep_path <- file.path(PATHS$docs, "endpoint_mapping_templates",
                     "endpoint_mapping_TCGA_BRCA.csv")
write_csv(ep_map, ep_path)
message("[02_TCGA] Endpoint map: ", ep_path)

message("\n[02_TCGA] Concluido: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message("Proximo passo: scripts/03_expression_preprocess_TCGA_BRCA.R")
