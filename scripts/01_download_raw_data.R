# =============================================================================
# SCRIPT: 01_download_raw_data.R
# PURPOSE: Download automático de dados clínicos e de expressão de fontes
#          públicas (GEO, GDC/TCGA, cBioPortal) + verificação SHA-256.
# PROJETO: Core-PAM (Memorial v6.1 / Freeze Core-PAM)
# NÍVEL:   A1-Proof — rastreabilidade forense completa
#
# FONTES:
#   GSE96058  (SCAN-B treino + validação) — GEO via GEOquery
#   METABRIC  (validação microarray)      — cBioPortal (acesso aberto)
#   TCGA-BRCA (validação RNA-seq)         — GDC via TCGAbiolinks
#   GSE20685  (Taiwan, validação)         — GEO via GEOquery
#
# DEPENDÊNCIAS R:
#   CRAN:        tidyverse, digest, arrow, httr, curl
#   Bioconductor: GEOquery, TCGAbiolinks, SummarizedExperiment
#
# EXECUÇÃO:
#   Rscript scripts/01_download_raw_data.R
#
# SAÍDAS (todas em 01_Base_Pura_CorePAM/RAW/<COHORT>/):
#   Expressão + clínica raw de cada coorte
#   registry/study_registry.csv  — append de cada artefato
#   01_docs/registry/data_lake_audit_report.csv — relatório final
#   results/supp/leakage_check_scanb_vs_gse96058.csv — anti-leakage
#
# REGRAS (Memorial v6.1):
#   - RAW é somente leitura após ingestão (não sobrescrever).
#   - Warning de leitura = erro (strict I/O via 00_setup.R).
#   - Nenhum número de genes nos nomes de arquivo.
#   - SCAN-B (treino) e GSE96058 (validação) compartilham acesso GEO mas
#     têm amostras DISTINTAS — split feito em 02_harmonize_clinical_SCANB.R.
# =============================================================================

source("scripts/00_setup.R")

# Load packages with warn=0: suppressPackageStartupMessages() only suppresses
# message() calls, NOT warning() calls. Under warn=2 (strict mode set by 00_setup.R),
# internal package-load warnings (e.g. GEOquery libdeflate) become hard errors.
# Pattern: save warn level, set to 0 for library(), restore immediately after.
old_warn <- getOption("warn"); options(warn = 0)
suppressPackageStartupMessages({
  library(GEOquery)
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(httr)
  library(curl)
})
options(warn = old_warn)

SCRIPT_NAME <- "01_download_raw_data.R"

# =============================================================================
# HELPER: download_file_safe
#   Baixa um arquivo via URL com mode="wb" (Windows-safe), verifica tamanho,
#   calcula SHA-256 e registra no study_registry.csv.
#   Se arquivo já existe e size > 0, pula o download (idempotente).
# =============================================================================
download_file_safe <- function(url, destfile, cohort, file_type,
                               overwrite = FALSE) {
  dir.create(dirname(destfile), showWarnings = FALSE, recursive = TRUE)

  if (file.exists(destfile) && file.info(destfile)$size > 0 && !overwrite) {
    message(sprintf("  [SKIP] Ja existe: %s", basename(destfile)))
    h    <- sha256_file(destfile)
    size <- file.info(destfile)$size / 1024^2
    registry_append(cohort, file_type, destfile, h, "INTEGRO_CACHED",
                    SCRIPT_NAME, size)
    return(invisible(h))
  }

  message(sprintf("  [DOWN] %s  ->  %s", url, basename(destfile)))
  resp <- tryCatch(
    httr::GET(url,
              httr::write_disk(destfile, overwrite = TRUE),
              httr::progress(),
              httr::timeout(600)),
    error = function(e) {
      stop("Falha no download de ", url, "\n", conditionMessage(e))
    }
  )
  if (httr::http_error(resp)) {
    stop("HTTP ", httr::status_code(resp), " ao baixar: ", url)
  }

  info <- file.info(destfile)
  if (is.na(info$size) || info$size == 0L) {
    stop("Download resultou em arquivo vazio: ", destfile)
  }

  h    <- sha256_file(destfile)
  size <- info$size / 1024^2
  message(sprintf("  [OK]   %.1f MB | SHA256: %s", size, h))

  registry_append(cohort, file_type, destfile, h, "INTEGRO", SCRIPT_NAME, size)
  invisible(h)
}

# =============================================================================
# HELPER: geo_supp_download
#   Baixa arquivos suplementares de um acesso GEO para o diretório de destino.
#   Retorna vetor de caminhos dos arquivos baixados.
# =============================================================================
geo_supp_download <- function(geo_acc, dest_dir, cohort) {
  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
  message(sprintf("\n  [GEO-SUPP] Baixando suplementares de %s ...", geo_acc))

  # GEOquery salva em subpasta com nome do acesso; mover depois
  tmp_base <- tempdir()
  old_warn <- getOption("warn"); options(warn = 0)
  files <- tryCatch(
    GEOquery::getGEOSuppFiles(geo_acc, baseDir = tmp_base, fetch_files = TRUE),
    error = function(e) {
      options(warn = old_warn)
      stop("Falha ao baixar suplementares de ", geo_acc, ": ", conditionMessage(e))
    }
  )
  options(warn = old_warn)

  src_dir  <- file.path(tmp_base, geo_acc)
  src_files <- list.files(src_dir, full.names = TRUE)
  out_paths <- file.path(dest_dir, basename(src_files))

  for (i in seq_along(src_files)) {
    file.copy(src_files[i], out_paths[i], overwrite = TRUE)
    info <- file.info(out_paths[i])
    h    <- sha256_file(out_paths[i])
    registry_append(cohort, paste0("Expression_supp_", i),
                    out_paths[i], h, "INTEGRO", SCRIPT_NAME,
                    info$size / 1024^2)
    message(sprintf("  [OK] %s | %.1f MB", basename(out_paths[i]),
                    info$size / 1024^2))
  }
  invisible(out_paths)
}

# =============================================================================
# HELPER: geo_clinical_download
#   Baixa a series matrix de um acesso GEO e salva como RDS (pData).
# =============================================================================
geo_clinical_download <- function(geo_acc, dest_rds, cohort) {
  dir.create(dirname(dest_rds), showWarnings = FALSE, recursive = TRUE)

  if (file.exists(dest_rds) && file.info(dest_rds)$size > 0) {
    message(sprintf("  [SKIP] Clinica ja existe: %s", basename(dest_rds)))
    h <- sha256_file(dest_rds)
    registry_append(cohort, "Clinical_raw", dest_rds, h,
                    "INTEGRO_CACHED", SCRIPT_NAME,
                    file.info(dest_rds)$size / 1024^2)
    return(invisible(h))
  }

  message(sprintf("  [GEO-CLIN] Baixando series matrix de %s ...", geo_acc))
  old_warn <- getOption("warn"); options(warn = 0)
  gse <- tryCatch(
    GEOquery::getGEO(geo_acc, GSEMatrix = TRUE, getGPL = FALSE),
    error = function(e) {
      options(warn = old_warn)
      stop("Falha ao baixar GEO series matrix de ", geo_acc, ": ",
           conditionMessage(e))
    }
  )
  options(warn = old_warn)

  # Extrair pData (metadados clínicos) do primeiro elemento
  pheno <- Biobase::pData(gse[[1]])
  saveRDS(pheno, dest_rds)

  h    <- sha256_file(dest_rds)
  size <- file.info(dest_rds)$size / 1024^2
  message(sprintf("  [OK] %d amostras | %.1f MB | SHA256: %s",
                  nrow(pheno), size, h))
  registry_append(cohort, "Clinical_raw", dest_rds, h,
                  "INTEGRO", SCRIPT_NAME, size)
  invisible(h)
}

# =============================================================================
# SEÇÃO 1: GSE96058 — SCAN-B (treino + validação GSE96058)
#   Plataforma: RNA-seq (Illumina HiSeq 2000/2500)
#   Expressão:  matriz de contagens + FPKM (suplementares GEO)
#   Clínica:    series matrix pData
#   Nota: o split treino (SCAN-B) vs validação (GSE96058) é feito em
#         02_harmonize_clinical_SCANB.R com base nos IDs de paciente.
# =============================================================================
message("\n", strrep("=", 60))
message("[01_download] SECAO 1: GSE96058 (SCAN-B + validacao)")
message(strrep("=", 60))

GSE96058_RAW_DIR <- raw_cohort("GSE96058")

# 1a. Arquivos suplementares (expressão — conta genes, FPKM)
geo_supp_download("GSE96058", GSE96058_RAW_DIR, cohort = "GSE96058")

# 1b. Clínica via series matrix
geo_clinical_download(
  geo_acc  = "GSE96058",
  dest_rds = file.path(GSE96058_RAW_DIR, "GSE96058_clinical_raw.rds"),
  cohort   = "GSE96058"
)

# =============================================================================
# SEÇÃO 2: METABRIC — cBioPortal (acesso aberto)
#   Plataforma: Illumina HT-12 v3 microarray
#   Expressão:  data_mrna_illumina_microarray.txt  (log2-intensidade)
#   Clínica:    data_clinical_patient.txt + data_clinical_sample.txt
#   Fonte:      cBioPortal Data Hub (AWS S3, sem autenticação)
#   Referência: Curtis et al. Nature 2012; Pereira et al. Nat Commun 2016
# =============================================================================
message("\n", strrep("=", 60))
message("[01_download] SECAO 2: METABRIC (cBioPortal)")
message(strrep("=", 60))

METABRIC_RAW_DIR  <- raw_cohort("METABRIC")
METABRIC_TAR      <- file.path(METABRIC_RAW_DIR, "brca_metabric.tar.gz")
METABRIC_BASE_URL <- "https://cbioportal-datahub.s3.amazonaws.com/brca_metabric.tar.gz"

# 2a. Download do pacote completo
download_file_safe(
  url      = METABRIC_BASE_URL,
  destfile = METABRIC_TAR,
  cohort   = "METABRIC",
  file_type = "Package_tarball"
)

# 2b. Extração
message("  [EXTRACT] Extraindo METABRIC...")
old_warn <- getOption("warn"); options(warn = 0)
utils::untar(METABRIC_TAR, exdir = METABRIC_RAW_DIR)
options(warn = old_warn)

# 2c. Localizar e registrar arquivos extraídos
metabric_files <- list(
  expression = file.path(METABRIC_RAW_DIR, "brca_metabric",
                         "data_mrna_illumina_microarray.txt"),
  clin_patient = file.path(METABRIC_RAW_DIR, "brca_metabric",
                           "data_clinical_patient.txt"),
  clin_sample  = file.path(METABRIC_RAW_DIR, "brca_metabric",
                           "data_clinical_sample.txt")
)

for (ftype in names(metabric_files)) {
  fpath <- metabric_files[[ftype]]
  if (!file.exists(fpath)) {
    warning("Arquivo METABRIC nao encontrado apos extracao: ", fpath)
    next
  }
  h    <- sha256_file(fpath)
  size <- file.info(fpath)$size / 1024^2
  registry_append("METABRIC", paste0("Extracted_", ftype),
                  fpath, h, "INTEGRO", SCRIPT_NAME, size)
  message(sprintf("  [OK] %s | %.1f MB", basename(fpath), size))
}

# =============================================================================
# SEÇÃO 3: TCGA-BRCA — GDC / TCGAbiolinks
#   Plataforma: RNA-seq (Illumina HiSeq; STAR alinhador)
#   Expressão:  STAR raw counts (unstranded)
#   Clínica:    GDC clinical XML (harmonized)
#   Nota: download pesado (~5-10 GB); requer conexão estável.
#         TCGAbiolinks salva em subpasta GDCdata/ por default.
# =============================================================================
message("\n", strrep("=", 60))
message("[01_download] SECAO 3: TCGA-BRCA (GDC / TCGAbiolinks)")
message(strrep("=", 60))

TCGA_RAW_DIR <- raw_cohort("TCGA_BRCA")
dir.create(TCGA_RAW_DIR, showWarnings = FALSE, recursive = TRUE)

# 3a. Query expressão (STAR counts)
message("  [GDC] Construindo query de expressao TCGA-BRCA...")
old_warn <- getOption("warn"); options(warn = 0)
query_expr <- tryCatch(
  TCGAbiolinks::GDCquery(
    project            = "TCGA-BRCA",
    data.category      = "Transcriptome Profiling",
    data.type          = "Gene Expression Quantification",
    workflow.type      = "STAR - Counts",
    sample.type        = "Primary Tumor"
  ),
  error = function(e) {
    options(warn = old_warn)
    stop("Falha na GDCquery TCGA-BRCA expressao: ", conditionMessage(e))
  }
)
options(warn = old_warn)

message(sprintf("  [GDC] %d arquivos encontrados para download.",
                nrow(TCGAbiolinks::getResults(query_expr))))

# 3b. Download (idempotente — retoma se já parcialmente baixado)
message("  [GDC] Iniciando download (pode demorar >30 min)...")
old_warn <- getOption("warn"); options(warn = 0)
tryCatch(
  TCGAbiolinks::GDCdownload(
    query     = query_expr,
    directory = TCGA_RAW_DIR,
    files.per.chunk = 10
  ),
  error = function(e) {
    options(warn = old_warn)
    stop("Falha em GDCdownload TCGA-BRCA: ", conditionMessage(e))
  }
)
options(warn = old_warn)

# 3c. Preparar SummarizedExperiment e salvar como RDS
message("  [GDC] Preparando SummarizedExperiment e salvando como RDS...")
old_warn <- getOption("warn"); options(warn = 0)
tcga_se <- tryCatch(
  TCGAbiolinks::GDCprepare(query_expr, directory = TCGA_RAW_DIR),
  error = function(e) {
    options(warn = old_warn)
    stop("Falha em GDCprepare TCGA-BRCA: ", conditionMessage(e))
  }
)
options(warn = old_warn)

tcga_expr_rds <- file.path(TCGA_RAW_DIR, "TCGA_BRCA_SE_counts_raw.rds")
saveRDS(tcga_se, tcga_expr_rds)
h    <- sha256_file(tcga_expr_rds)
size <- file.info(tcga_expr_rds)$size / 1024^2
registry_append("TCGA_BRCA", "Expression_SE_counts", tcga_expr_rds,
                h, "INTEGRO", SCRIPT_NAME, size)
message(sprintf("  [OK] SE salvo: %d genes x %d amostras | %.1f MB",
                nrow(tcga_se), ncol(tcga_se), size))

# 3d. Clínica TCGA (harmonized clinical)
message("  [GDC] Baixando dados clinicos TCGA-BRCA...")
old_warn <- getOption("warn"); options(warn = 0)
tcga_clin <- tryCatch(
  TCGAbiolinks::GDCquery_clinic("TCGA-BRCA", type = "clinical"),
  error = function(e) {
    options(warn = old_warn)
    stop("Falha ao baixar clinica TCGA-BRCA: ", conditionMessage(e))
  }
)
options(warn = old_warn)

tcga_clin_rds <- file.path(TCGA_RAW_DIR, "TCGA_BRCA_clinical_raw.rds")
saveRDS(tcga_clin, tcga_clin_rds)
h    <- sha256_file(tcga_clin_rds)
size <- file.info(tcga_clin_rds)$size / 1024^2
registry_append("TCGA_BRCA", "Clinical_raw", tcga_clin_rds,
                h, "INTEGRO", SCRIPT_NAME, size)
message(sprintf("  [OK] Clinica: %d pacientes | %.1f MB | SHA256: %s",
                nrow(tcga_clin), size, h))

# =============================================================================
# SEÇÃO 4: GSE20685 — Taiwan (validação microarray)
#   Plataforma: Affymetrix Human Genome U133A (GPL96)
#   Expressão:  series matrix (intensidades log2)
#   Clínica:    series matrix pData
#   Referência: Lu et al. BMC Med Genomics 2012
# =============================================================================
message("\n", strrep("=", 60))
message("[01_download] SECAO 4: GSE20685 (Taiwan microarray)")
message(strrep("=", 60))

GSE20685_RAW_DIR <- raw_cohort("GSE20685")

# 4a. Suplementares (matriz de expressão)
geo_supp_download("GSE20685", GSE20685_RAW_DIR, cohort = "GSE20685")

# 4b. Clínica via series matrix
geo_clinical_download(
  geo_acc  = "GSE20685",
  dest_rds = file.path(GSE20685_RAW_DIR, "GSE20685_clinical_raw.rds"),
  cohort   = "GSE20685"
)

# =============================================================================
# SEÇÃO 5: RELATÓRIO FINAL E GO / NO-GO
# =============================================================================
message("\n", strrep("=", 60))
message("[01_download] SECAO 5: RELATORIO FINAL")
message(strrep("=", 60))

# Ler registry e filtrar apenas entradas deste script
old_warn <- getOption("warn"); options(warn = 0)
if (file.exists(PATHS$run_registry)) {
  df_reg <- readr::read_csv(PATHS$run_registry, show_col_types = FALSE)
  df_this <- df_reg |> filter(script == SCRIPT_NAME)
} else {
  df_this <- tibble()
}
options(warn = old_warn)

# Salvar relatório consolidado
audit_out <- file.path(PATHS$registry_docs, "data_lake_audit_report.csv")
write_csv(df_this, audit_out)
message("Relatorio salvo em: ", audit_out)

if (nrow(df_this) > 0) {
  n_ok  <- sum(grepl("INTEGRO", df_this$status))
  n_bad <- sum(!grepl("INTEGRO", df_this$status))
  message(sprintf("\nArtefatos integros : %d", n_ok))
  message(sprintf("Artefatos com falha: %d", n_bad))
  print(df_this |> select(cohort, file_type, status, size_mb))
}

# =============================================================================
# SEÇÃO 6: LEAKAGE-PROOF CHECK — SCAN-B (treino) vs GSE96058 (validação)
#   Memorial v6.1 §3.2 — interseção de IDs deve ser 0.
#   ATENÇÃO: este check opera sobre os pData (series matrix) do GEO.
#   O split definitivo (por barcode/patient_id) ocorre em 02_harmonize_*.R.
#   Aqui verificamos que os GSM IDs dos dois subconjuntos não se sobrepõem.
# =============================================================================
message("\n", strrep("=", 60))
message("[01_download] SECAO 6: LEAKAGE-PROOF CHECK")
message(strrep("=", 60))

scanb_clin_path <- file.path(raw_cohort("GSE96058"), "GSE96058_clinical_raw.rds")
leakage_out     <- file.path(PATHS$results$supp,
                             "leakage_check_scanb_vs_gse96058.csv")

# Nota: o dataset GSE96058 contém TODAS as amostras SCAN-B.
# O split SCAN-B (treino) vs validação é feito por coluna de status
# presente nos metadados (ex.: coluna "scan.b.training.validation").
# Aqui verificamos a coerência dos grupos dentro do dataset.
if (file.exists(scanb_clin_path)) {
  old_warn <- getOption("warn"); options(warn = 0)
  pheno <- readRDS(scanb_clin_path)
  options(warn = old_warn)

  # Detectar coluna de split (ajustar se nome diferir)
  split_col <- intersect(
    c("scan.b.training.validation", "group", "training_validation",
      "characteristics_ch1.6"),
    tolower(names(pheno))
  )

  if (length(split_col) > 0) {
    col_use   <- names(pheno)[tolower(names(pheno)) == split_col[1]][1]
    groups    <- table(pheno[[col_use]])
    train_ids <- rownames(pheno)[grepl("train",  tolower(pheno[[col_use]]))]
    valid_ids <- rownames(pheno)[grepl("valid|test", tolower(pheno[[col_use]]))]
    overlap   <- intersect(train_ids, valid_ids)

    leakage_df <- tibble(
      timestamp            = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      script               = SCRIPT_NAME,
      split_column_used    = col_use,
      n_train_gsm          = length(train_ids),
      n_validation_gsm     = length(valid_ids),
      n_overlap            = length(overlap),
      overlap_ids          = paste(overlap, collapse = ";"),
      sha256_train_ids     = digest::digest(sort(train_ids),  algo = "sha256"),
      sha256_valid_ids     = digest::digest(sort(valid_ids),  algo = "sha256"),
      result               = if (length(overlap) == 0) "PASS_NO_LEAKAGE"
                             else "FAIL_LEAKAGE_DETECTED"
    )
    write_csv(leakage_df, leakage_out)

    if (length(overlap) == 0) {
      message(sprintf("[leakage] PASS — treino: %d | validacao: %d | overlap: 0",
                      length(train_ids), length(valid_ids)))
      message("  Artefato: ", leakage_out)
    } else {
      stop(sprintf(
        "[leakage] FAIL — %d amostras sobrepostas entre treino e validacao.\n%s",
        length(overlap), paste(overlap, collapse = ", ")
      ))
    }
  } else {
    message("[leakage] AVISO: coluna de split nao detectada automaticamente.")
    message("  Inspecione pData(GSE96058) e ajuste split_col neste script.")
    message("  Colunas disponiveis: ", paste(names(pheno)[1:10], collapse = ", "))
  }
} else {
  message("[leakage] Adiado — arquivo clinico GSE96058 ainda ausente.")
}

message("\n[01_download] Concluido: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message("Proximo passo: scripts/02_harmonize_clinical_<COHORT>.R")
