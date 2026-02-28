# =============================================================================
# SCRIPT: 03_expression_preprocess_SCANB.R
# PURPOSE: Pré-processamento de expressão SCAN-B (subconjunto treino).
#          RNA-seq: raw counts → edgeR TMM → logCPM; mapeamento HGNC.
# PROJETO: Core-PAM (Memorial v6.1 §4.2)
#
# INPUTS:
#   01_Base_Pura_CorePAM/RAW/GSE96058/  (arquivos suplementares GEO)
#   01_Base_Pura_CorePAM/PROCESSED/SCANB/clinical_FINAL.parquet  (IDs treino)
#
# OUTPUTS:
#   01_Base_Pura_CorePAM/PROCESSED/SCANB/expression_genelevel_preZ.parquet
#   results/supp/gene_mapping_audit_SCANB.csv
#
# PIPELINE (congelado):
#   raw counts → filtrar amostras SCANB → edgeR TMM → logCPM → HGNC map →
#   colapso por variância → expression_genelevel_preZ (genes × amostras)
#
# NOTA: Z-score intra-coorte é feito em 06_zscore_and_score_SCANB.R.
# =============================================================================

source("scripts/00_setup.R")
source("scripts/03_utils_gene_mapping.R")

suppressPackageStartupMessages(library(edgeR))

SCRIPT_NAME <- "03_expression_preprocess_SCANB.R"
COHORT      <- "SCANB"

# =============================================================================
# 1) CARREGAR IDs SCANB (do clinical_FINAL para filtrar amostras de treino)
# =============================================================================
clin_path <- file.path(proc_cohort(COHORT), "clinical_FINAL.parquet")
clin      <- strict_parquet(clin_path)
scanb_ids <- normalize_id(clin$sample_id)
message(sprintf("[03_SCANB] IDs treino carregados: %d amostras", length(scanb_ids)))

# =============================================================================
# 2) LOCALIZAR ARQUIVO DE EXPRESSÃO RAW NO DATA LAKE
#    GSE96058 disponibiliza vários formatos suplementares.
#    Prioridade: (1) arquivo de counts brutos; (2) matriz de genes filtrada.
# =============================================================================
raw_dir   <- raw_cohort("GSE96058")
raw_files <- list.files(raw_dir, full.names = TRUE, recursive = FALSE)
message("[03_SCANB] Arquivos disponíveis em RAW/GSE96058:")
print(basename(raw_files))

# Detectar arquivo de counts (prioridade) vs gene expression matrix
counts_file <- raw_files[grep("count|raw|RAW\\.tar", basename(raw_files),
                               ignore.case = TRUE)][1]
expr_file   <- raw_files[grep("gene_expression|genematrix|expression",
                               basename(raw_files), ignore.case = TRUE)][1]

# =============================================================================
# 3A) LEITURA DE COUNTS BRUTOS (se disponível — preferido pelo Memorial)
# =============================================================================
if (!is.na(counts_file) && file.exists(counts_file)) {
  message("[03_SCANB] Usando arquivo de counts brutos: ", basename(counts_file))

  if (grepl("\\.tar$", counts_file)) {
    # GSE96058_RAW.tar → extrair para temp e combinar counts individuais
    tmp_dir <- file.path(raw_dir, "RAW_extracted")
    dir.create(tmp_dir, showWarnings = FALSE)
    old_warn <- getOption("warn"); options(warn = 0)
    utils::untar(counts_file, exdir = tmp_dir)
    options(warn = old_warn)

    count_files <- list.files(tmp_dir, pattern = "\\.txt(\\.gz)?$",
                               full.names = TRUE, recursive = TRUE)
    message(sprintf("[03_SCANB] %d arquivos de counts individuais extraídos.",
                    length(count_files)))

    # Ler e combinar: cada arquivo = uma amostra (formato: GeneID \t Count)
    counts_list <- lapply(count_files, function(f) {
      old_warn <- getOption("warn"); options(warn = 0)
      d <- read_tsv(f, col_names = c("gene_id", "count"),
                    show_col_types = FALSE, skip = 0)
      options(warn = old_warn)
      sample_name <- sub("(_count|_counts|\\.txt.*)", "",
                         basename(f), ignore.case = TRUE)
      setNames(d$count, d$gene_id)
    })
    sample_names <- sub("(_count|_counts|\\.txt.*)", "",
                        basename(count_files), ignore.case = TRUE)

    all_genes  <- Reduce(intersect, lapply(counts_list, names))
    count_mat  <- do.call(cbind, lapply(counts_list, function(x) x[all_genes]))
    colnames(count_mat) <- sample_names
    count_mat  <- count_mat[!grepl("^__", rownames(count_mat)), ] # remover linhas de QC HTSeq

    input_type <- "raw_counts"

  } else if (grepl("\\.gz$|\\.csv$|\\.tsv$", counts_file)) {
    old_warn <- getOption("warn"); options(warn = 0)
    count_df <- read_csv(counts_file, show_col_types = FALSE)
    options(warn = old_warn)
    gene_col  <- names(count_df)[1]
    count_mat <- as.matrix(count_df[, -1])
    rownames(count_mat) <- count_df[[gene_col]]
    input_type <- "raw_counts"
  }

# =============================================================================
# 3B) FALLBACK: gene expression matrix (se counts não disponíveis)
# =============================================================================
} else if (!is.na(expr_file) && file.exists(expr_file)) {
  message("[03_SCANB] FALLBACK: usando gene expression matrix (não são counts brutos).")
  message("  Nota: TMM não será aplicado; assumir pré-normalizado.")

  old_warn <- getOption("warn"); options(warn = 0)
  expr_df <- read_csv(expr_file, show_col_types = FALSE)
  options(warn = old_warn)
  gene_col  <- names(expr_df)[1]
  count_mat <- as.matrix(expr_df[, -1])
  rownames(count_mat) <- expr_df[[gene_col]]
  storage.mode(count_mat) <- "numeric"
  input_type <- "preprocessed_matrix"

} else {
  stop("[03_SCANB] Nenhum arquivo de expressão encontrado em RAW/GSE96058/.\n",
       "Execute 01_download_raw_data.R primeiro.")
}

message(sprintf("[03_SCANB] Matriz raw: %d genes x %d amostras | tipo: %s",
                nrow(count_mat), ncol(count_mat), input_type))

# =============================================================================
# 4) FILTRAR AMOSTRAS SCANB (treino)
# =============================================================================
# Normalizar nomes de colunas para cruzar com sample_ids clínicos
col_ids_norm <- normalize_id(colnames(count_mat))
keep_cols    <- col_ids_norm %in% scanb_ids
n_matched    <- sum(keep_cols)

message(sprintf("[03_SCANB] Match amostras: %d de %d IDs clínicos encontrados na matriz.",
                n_matched, length(scanb_ids)))

if (n_matched == 0) {
  stop("[03_SCANB] Nenhuma amostra SCANB encontrada na matriz de expressão.\n",
       "Verifique se os IDs em clinical_FINAL.parquet correspondem aos nomes de coluna.")
}
if (n_matched < 0.9 * length(scanb_ids)) {
  warning(sprintf("[03_SCANB] Apenas %.1f%% das amostras clínicas encontradas na matriz.",
                  100 * n_matched / length(scanb_ids)))
}

count_mat <- count_mat[, keep_cols, drop = FALSE]
colnames(count_mat) <- col_ids_norm[keep_cols]

message(sprintf("[03_SCANB] Matriz filtrada: %d genes x %d amostras SCANB",
                nrow(count_mat), ncol(count_mat)))

# =============================================================================
# 5) NORMALIZAÇÃO RNA-seq: edgeR TMM → logCPM (intra-coorte)
# =============================================================================
if (input_type == "raw_counts") {
  message("[03_SCANB] Aplicando edgeR TMM → logCPM...")

  # Filtro de genes com expressão muito baixa (recomendação edgeR)
  dge      <- edgeR::DGEList(counts = count_mat)
  keep     <- edgeR::filterByExpr(dge)
  dge      <- dge[keep, , keep.lib.sizes = FALSE]
  message(sprintf("[03_SCANB] Genes apos filtro edgeR: %d (removidos: %d)",
                  sum(keep), sum(!keep)))

  dge      <- edgeR::calcNormFactors(dge, method = "TMM")
  expr_mat <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

} else {
  # Já normalizado — usar como está (assumir log2)
  expr_mat <- count_mat
  message("[03_SCANB] Matriz pré-normalizada usada diretamente (sem TMM).")
}

message(sprintf("[03_SCANB] Matriz logCPM: %d genes x %d amostras",
                nrow(expr_mat), ncol(expr_mat)))

# =============================================================================
# 6) MAPEAMENTO HGNC
#    Detectar tipo de ID (Ensembl vs Symbol) e mapear
# =============================================================================
gene_ids_raw <- rownames(expr_mat)
is_ensembl   <- grepl("^ENSG", gene_ids_raw[1])

if (is_ensembl) {
  message("[03_SCANB] IDs detectados: Ensembl → mapeando para HGNC...")
  hgnc_map <- build_hgnc_map_ensembl(gene_ids_raw)

  # Resolver duplicatas: manter mapeamento 1:1 (remover Ensembl sem symbol)
  hgnc_map <- hgnc_map[hgnc_map$hgnc_symbol != "" & !is.na(hgnc_map$hgnc_symbol), ]
  hgnc_map <- hgnc_map[!duplicated(hgnc_map$ensembl_gene_id), ]

  gene_symbols <- hgnc_map$hgnc_symbol[match(gene_ids_raw, hgnc_map$ensembl_gene_id)]

} else {
  message("[03_SCANB] IDs detectados: Gene symbols → padronizando via HGNChelper...")
  hgnc_res     <- build_hgnc_map_symbol(gene_ids_raw)
  gene_symbols <- hgnc_res$hgnc_symbol
}

n_mapped   <- sum(!is.na(gene_symbols) & gene_symbols != "")
n_unmapped <- sum(is.na(gene_symbols)  | gene_symbols == "")
message(sprintf("[03_SCANB] Mapeamento: %d mapeados | %d nao mapeados (descartados)",
                n_mapped, n_unmapped))

# Manter apenas genes mapeados
keep_genes   <- !is.na(gene_symbols) & gene_symbols != ""
expr_mat     <- expr_mat[keep_genes, , drop = FALSE]
gene_symbols <- gene_symbols[keep_genes]

# =============================================================================
# 7) COLAPSO DE PROBES: probe de maior variância por gene (§4.4)
# =============================================================================
n_before <- nrow(expr_mat)
expr_mat  <- collapse_probes_by_variance(expr_mat, gene_symbols)
n_after   <- nrow(expr_mat)
message(sprintf("[03_SCANB] Colapso probes: %d → %d genes únicos (HGNC)",
                n_before, n_after))

# =============================================================================
# 8) AUDIT DE MAPEAMENTO
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

message(sprintf("[03_SCANB] PAM50: %d/%d genes presentes | ausentes: %s",
                length(pam50_present), length(pam50_genes),
                if (length(pam50_missing) == 0) "nenhum"
                else paste(pam50_missing, collapse = ",")))

audit_df <- tibble(
  cohort          = COHORT,
  gene            = rownames(expr_mat),
  present         = TRUE,
  in_pam50        = rownames(expr_mat) %in% pam50_genes
) |>
  bind_rows(tibble(
    cohort   = COHORT,
    gene     = pam50_missing,
    present  = FALSE,
    in_pam50 = TRUE
  ))

save_gene_mapping_audit(audit_df, COHORT, SCRIPT_NAME)

# =============================================================================
# 9) SALVAR expression_genelevel_preZ.parquet
#    Formato: linhas = genes (rownames), colunas = amostras + coluna "gene"
# =============================================================================
dest_dir <- proc_cohort(COHORT)
dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(dest_dir, "expression_genelevel_preZ.parquet")

expr_df <- as.data.frame(expr_mat)
expr_df <- tibble::rownames_to_column(expr_df, var = "gene")

arrow::write_parquet(expr_df, out_path)
h    <- sha256_file(out_path)
size <- file.info(out_path)$size / 1024^2
registry_append(COHORT, "Expression_preZ", out_path, h,
                "INTEGRO", SCRIPT_NAME, size)
message(sprintf("[03_SCANB] Salvo: %s | %d genes x %d amostras | %.2f MB | SHA256: %s",
                out_path, nrow(expr_mat), ncol(expr_mat), size, h))

message("\n[03_SCANB] Concluido: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message("Proximo passo: scripts/03_expression_preprocess_GSE96058.R")
