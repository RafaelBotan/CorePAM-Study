# =============================================================================
# SCRIPT: 03_expression_preprocess_GSE20685.R
# PURPOSE: Pré-processamento de expressão GSE20685 (Taiwan, microarray).
#          Affymetrix HGU133A: log2 intensity as-is → HGNC map.
# PROJETO: Core-PAM (Memorial v6.1 §4.3)
#
# INPUTS:
#   01_Base_Pura_CorePAM/RAW/GSE20685/  (suplementares GEO — .gz)
#   01_Base_Pura_CorePAM/PROCESSED/GSE20685/clinical_FINAL.parquet
#
# OUTPUTS:
#   01_Base_Pura_CorePAM/PROCESSED/GSE20685/expression_genelevel_preZ.parquet
#   results/supp/gene_mapping_audit_GSE20685.csv
#
# REGRAS (Memorial v6.1 §4.3):
#   - Microarray: log2-intensity as-is (sem TMM).
#   - Probes Affymetrix → HGNC via hgu133a.db.
#   - Múltiplas probes → maior variância intra-coorte.
#   - Z-score intra-coorte em 06_zscore_and_score_GSE20685.R.
# =============================================================================

source("scripts/00_setup.R")
source("scripts/03_utils_gene_mapping.R")

SCRIPT_NAME <- "03_expression_preprocess_GSE20685.R"
COHORT      <- "GSE20685"

# =============================================================================
# 1) LOCALIZAR ARQUIVO DE EXPRESSÃO
#    GSE20685 suplementares podem ser:
#    (a) série matrix comprimida (series_matrix.txt.gz) — contém expr + clínica
#    (b) arquivo de expressão separado (.csv.gz ou .txt.gz)
# =============================================================================
raw_dir   <- raw_cohort(COHORT)
raw_files <- list.files(raw_dir, full.names = TRUE, recursive = FALSE)
message("[03_GSE20685] Arquivos em RAW/GSE20685:")
print(basename(raw_files))

# Prioridade: arquivo de matriz de expressão suplementar
expr_file <- raw_files[grep("matrix|expression|series_matrix",
                             basename(raw_files), ignore.case = TRUE)][1]

if (is.na(expr_file) || !file.exists(expr_file)) {
  stop("[03_GSE20685] Arquivo de expressão não encontrado em RAW/GSE20685/.\n",
       "Execute 01_download_raw_data.R primeiro.")
}
message("[03_GSE20685] Usando arquivo: ", basename(expr_file))

# =============================================================================
# 2) LEITURA DA SÉRIE MATRIX GEO
#    O series_matrix.txt.gz contém blocos de metadados (linhas "!") e
#    depois a tabela de expressão.
# =============================================================================
old_warn <- getOption("warn"); options(warn = 0)

if (grepl("series_matrix", basename(expr_file), ignore.case = TRUE)) {
  # Ler via GEOquery (já carregado indiretamente)
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop("GEOquery necessário: BiocManager::install('GEOquery')")
  }
  message("[03_GSE20685] Lendo series_matrix via GEOquery...")
  gse       <- GEOquery::getGEO(filename = expr_file)
  expr_mat  <- Biobase::exprs(gse)
  probe_ids <- rownames(expr_mat)
  platform  <- annotation(gse)
  message(sprintf("[03_GSE20685] Platform: %s | %d probes x %d amostras",
                  platform, nrow(expr_mat), ncol(expr_mat)))

} else {
  # Arquivo tabular suplementar
  message("[03_GSE20685] Lendo arquivo tabular suplementar...")
  expr_df  <- read_tsv(expr_file, show_col_types = FALSE, name_repair = "unique")
  gene_col <- names(expr_df)[1]
  probe_ids <- expr_df[[gene_col]]
  expr_mat  <- as.matrix(expr_df[, -1])
  rownames(expr_mat) <- probe_ids
  storage.mode(expr_mat) <- "numeric"
  platform  <- "hgu133a"   # GSE20685 usa HGU133A
  message(sprintf("[03_GSE20685] %d probes x %d amostras", nrow(expr_mat), ncol(expr_mat)))
}
options(warn = old_warn)

# =============================================================================
# 3) FILTRAR AMOSTRAS COM CLÍNICA
# =============================================================================
clin_path <- file.path(proc_cohort(COHORT), "clinical_FINAL.parquet")
clin      <- strict_parquet(clin_path)
clin_ids  <- normalize_id(clin$sample_id)

col_ids_norm <- normalize_id(colnames(expr_mat))
keep_cols    <- col_ids_norm %in% clin_ids
n_matched    <- sum(keep_cols)

message(sprintf("[03_GSE20685] Match amostras: %d de %d IDs clínicos.",
                n_matched, length(clin_ids)))

if (n_matched == 0) {
  # Fallback: usar todos (série matrix inclui todas as amostras do estudo)
  message("[03_GSE20685] AVISO: IDs não cruzaram. Usando todas as amostras.")
  message("  Verifique se sample_id em clinical_FINAL corresponde a geo_accession.")
  keep_cols <- rep(TRUE, ncol(expr_mat))
}

expr_mat <- expr_mat[, keep_cols, drop = FALSE]
colnames(expr_mat) <- col_ids_norm[keep_cols]

# =============================================================================
# 4) VERIFICAR ESCALA LOG2
# =============================================================================
p95 <- quantile(expr_mat, 0.95, na.rm = TRUE)
if (p95 > 20) {
  message(sprintf("[03_GSE20685] p95=%.1f > 20 → aplicando log2(x + 1).", p95))
  expr_mat <- log2(expr_mat + 1)
} else {
  message(sprintf("[03_GSE20685] p95=%.1f ≤ 20 → dados já em log2 (as-is).", p95))
}

# =============================================================================
# 5) MAPEAMENTO PROBES AFFYMETRIX → HGNC
#    GSE20685: Affymetrix Human Genome U133A (GPL96) ou U133 Plus 2.0 (GPL570)
# =============================================================================
platform_pkg <- dplyr::case_when(
  grepl("GPL96|hgu133a|133a$",  tolower(platform)) ~ "hgu133a",
  grepl("GPL570|hgu133plus2",   tolower(platform)) ~ "hgu133plus2",
  grepl("GPL6947|hgu133b",      tolower(platform)) ~ "hgu133b",
  TRUE ~ "hgu133a"   # default para GSE20685
)
message(sprintf("[03_GSE20685] Platform detectado: %s → usando %s.db",
                platform, platform_pkg))

probe_map    <- build_hgnc_map_affymetrix(rownames(expr_mat), platform_pkg)
gene_symbols <- probe_map$SYMBOL[match(rownames(expr_mat), probe_map$PROBEID)]

# Validar symbols via HGNChelper
if (any(!is.na(gene_symbols))) {
  valid_symbols <- unique(gene_symbols[!is.na(gene_symbols)])
  hgnc_check   <- build_hgnc_map_symbol(valid_symbols)
  symbol_map   <- setNames(hgnc_check$hgnc_symbol, hgnc_check$original)
  gene_symbols <- ifelse(!is.na(gene_symbols),
                         symbol_map[gene_symbols],
                         NA_character_)
}

keep_genes   <- !is.na(gene_symbols) & gene_symbols != ""
n_unmapped   <- sum(!keep_genes)
expr_mat     <- expr_mat[keep_genes, , drop = FALSE]
gene_symbols <- gene_symbols[keep_genes]
message(sprintf("[03_GSE20685] Mapeados: %d probes | Descartados: %d",
                sum(keep_genes), n_unmapped))

# =============================================================================
# 6) COLAPSO DE PROBES (maior variância intra-coorte)
# =============================================================================
n_before <- nrow(expr_mat)
expr_mat  <- collapse_probes_by_variance(expr_mat, gene_symbols)
message(sprintf("[03_GSE20685] Colapso probes: %d → %d genes únicos HGNC",
                n_before, nrow(expr_mat)))

# =============================================================================
# 7) QC: remover genes com missingness > 20%
# =============================================================================
miss_frac <- rowMeans(is.na(expr_mat))
n_high_miss <- sum(miss_frac > 0.20)
if (n_high_miss > 0) {
  message(sprintf("[03_GSE20685] Removendo %d genes com >20%% missingness.",
                  n_high_miss))
  expr_mat <- expr_mat[miss_frac <= 0.20, , drop = FALSE]
}

# =============================================================================
# 8) AUDIT PAM50
# =============================================================================
pam50_genes <- c(
  "ACTR3B","ANLN","BAG1","BCL2","BIRC5","BLVRA","CCNB1","CCNE1","CDC20",
  "CDC6","CDH3","CENPF","CEP55","CXXC5","EGFR","ERBB2","ESR1","EXO1",
  "FGFR4","FOXA1","FOXC1","GPR160","GRB7","KIF2C","KRT14","KRT17","KRT5",
  "MAPT","MDM2","MELK","MIA","MKI67","MLPH","MMP11","MYBL2","MYC","NAT1",
  "NDC80","NUF2","ORC6","PGR","PHGDH","PTTG1","RRM2","SFRP1","SLC39A6",
  "TMEM45B","TYMS","UBE2C","UBE2T"
)
pam50_present <- intersect(pam50_genes, rownames(expr_mat))
pam50_missing <- setdiff(pam50_genes, rownames(expr_mat))
message(sprintf("[03_GSE20685] PAM50: %d/50 presentes | ausentes: %s",
                length(pam50_present),
                if (length(pam50_missing) == 0) "nenhum"
                else paste(pam50_missing, collapse = ",")))

audit_df <- tibble(
  cohort   = COHORT,
  gene     = rownames(expr_mat),
  present  = TRUE,
  in_pam50 = rownames(expr_mat) %in% pam50_genes
) |>
  bind_rows(tibble(cohort = COHORT, gene = pam50_missing,
                   present = FALSE, in_pam50 = TRUE))
save_gene_mapping_audit(audit_df, COHORT, SCRIPT_NAME)

# =============================================================================
# 9) SALVAR
# =============================================================================
dest_dir <- proc_cohort(COHORT)
dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(dest_dir, "expression_genelevel_preZ.parquet")

expr_df  <- tibble::rownames_to_column(as.data.frame(expr_mat), var = "gene")
arrow::write_parquet(expr_df, out_path)
h    <- sha256_file(out_path)
size <- file.info(out_path)$size / 1024^2
registry_append(COHORT, "Expression_preZ", out_path, h,
                "INTEGRO", SCRIPT_NAME, size)
message(sprintf("[03_GSE20685] Salvo: %d genes x %d amostras | %.2f MB | SHA256: %s",
                nrow(expr_mat), ncol(expr_mat), size, h))

message("\n[03_GSE20685] Concluido: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message("Proximo passo: scripts/04_gene_audit_freeze.R")
