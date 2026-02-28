# =============================================================================
# UTILITÁRIO: 03_utils_gene_mapping.R
# PURPOSE: Funções compartilhadas de mapeamento gênico (HGNC) e colapso de
#          probes para todos os scripts 03_expression_preprocess_<COHORT>.R
# PROJETO: Core-PAM (Memorial v6.1 §4.4)
#
# REGRAS:
#   - IDs (Symbol/Entrez/RefSeq/Ensembl) → HGNC official symbol.
#   - Múltiplas probes → gene: selecionar probe de maior variância intra-coorte.
#   - Salvar gene_mapping_audit_<COHORT>.csv por coorte.
#   - NÃO fazer Z-score aqui — output é pre-Z (logCPM ou log2 intensity).
# =============================================================================

suppressPackageStartupMessages({
  library(biomaRt)
  library(HGNChelper)
})

# --------------------------------------------------------------------------
# build_hgnc_map_ensembl
#   Converte Ensembl IDs para HGNC official symbols via biomaRt.
#   Cache: salva em registry para evitar re-download.
# --------------------------------------------------------------------------
build_hgnc_map_ensembl <- function(ensembl_ids,
                                   cache_path = file.path(
                                     PATHS$registry_docs,
                                     "hgnc_ensembl_map_cache.rds")) {
  if (file.exists(cache_path)) {
    message("  [gene_map] Usando cache HGNC Ensembl: ", cache_path)
    old_warn <- getOption("warn"); options(warn = 0)
    map <- readRDS(cache_path)
    options(warn = old_warn)
  } else {
    message("  [gene_map] Baixando mapa Ensembl→HGNC via biomaRt...")
    old_warn <- getOption("warn"); options(warn = 0)
    mart <- tryCatch(
      biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl"),
      error = function(e) {
        options(warn = old_warn)
        stop("biomaRt indisponivel: ", conditionMessage(e),
             "\nVerifique conexao ou use cache manual.")
      }
    )
    map <- biomaRt::getBM(
      attributes = c("ensembl_gene_id", "hgnc_symbol",
                     "entrezgene_id", "gene_biotype"),
      filters    = "ensembl_gene_id",
      values     = unique(ensembl_ids),
      mart       = mart
    )
    options(warn = old_warn)
    saveRDS(map, cache_path)
    message("  [gene_map] Cache salvo: ", cache_path)
  }
  map
}

# --------------------------------------------------------------------------
# build_hgnc_map_symbol
#   Padroniza gene symbols existentes usando HGNChelper (offline).
#   Retorna data.frame com original → hgnc_symbol.
# --------------------------------------------------------------------------
build_hgnc_map_symbol <- function(symbols) {
  message("  [gene_map] Padronizando symbols via HGNChelper...")
  old_warn <- getOption("warn"); options(warn = 0)
  res <- HGNChelper::checkGeneSymbols(
    symbols,
    unmapped.as.na = TRUE,
    species = "human"
  )
  options(warn = old_warn)
  data.frame(
    original    = res$x,
    hgnc_symbol = res$Suggested.Symbol,
    approved    = res$Approved,
    stringsAsFactors = FALSE
  )
}

# --------------------------------------------------------------------------
# build_hgnc_map_entrez
#   Converte Entrez IDs para HGNC symbols via org.Hs.eg.db.
# --------------------------------------------------------------------------
build_hgnc_map_entrez <- function(entrez_ids) {
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("Pacote org.Hs.eg.db necessario: ",
         "BiocManager::install('org.Hs.eg.db')")
  }
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  message("  [gene_map] Mapeando Entrez→HGNC via org.Hs.eg.db...")
  old_warn <- getOption("warn"); options(warn = 0)
  map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys    = as.character(unique(entrez_ids[!is.na(entrez_ids)])),
    columns = c("ENTREZID", "SYMBOL"),
    keytype = "ENTREZID"
  )
  options(warn = old_warn)
  map
}

# --------------------------------------------------------------------------
# build_hgnc_map_affymetrix
#   Converte probe IDs Affymetrix HGU133A/Plus2 para HGNC symbols.
# --------------------------------------------------------------------------
build_hgnc_map_affymetrix <- function(probe_ids, platform = "hgu133a") {
  pkg <- paste0(platform, ".db")
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Pacote %s necessario: BiocManager::install('%s')", pkg, pkg))
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  message(sprintf("  [gene_map] Mapeando probes %s→SYMBOL via %s...",
                  platform, pkg))
  db_obj <- get(pkg)
  old_warn <- getOption("warn"); options(warn = 0)
  map <- AnnotationDbi::select(
    db_obj,
    keys    = unique(probe_ids),
    columns = c("PROBEID", "SYMBOL", "ENTREZID"),
    keytype = "PROBEID"
  )
  options(warn = old_warn)
  map
}

# --------------------------------------------------------------------------
# collapse_probes_by_variance
#   Para matrizes com múltiplas linhas por gene:
#   seleciona a linha (probe) com maior variância intra-coorte.
#   Input:  mat     = matrix/data.frame (linhas=probes, colunas=amostras)
#           gene_id = vetor de gene symbols (mesmo comprimento que nrow(mat))
#   Output: matrix com uma linha por gene (probe de maior variância)
# --------------------------------------------------------------------------
collapse_probes_by_variance <- function(mat, gene_id) {
  stopifnot(length(gene_id) == nrow(mat))

  df <- data.frame(
    gene     = gene_id,
    row_idx  = seq_len(nrow(mat)),
    variance = apply(mat, 1, var, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  # Remover genes NA ou vazios
  df <- df[!is.na(df$gene) & nchar(trimws(df$gene)) > 0, ]

  # Selecionar probe de maior variância por gene
  best <- df[order(df$gene, -df$variance), ]
  best <- best[!duplicated(best$gene), ]

  mat_collapsed <- mat[best$row_idx, , drop = FALSE]
  rownames(mat_collapsed) <- best$gene
  mat_collapsed
}

# --------------------------------------------------------------------------
# save_gene_mapping_audit
#   Salva gene_mapping_audit_<cohort>.csv e registra hash.
# --------------------------------------------------------------------------
save_gene_mapping_audit <- function(audit_df, cohort, script_name) {
  dest_dir <- file.path(PATHS$results$supp)
  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
  out_path <- file.path(dest_dir,
                        sprintf("gene_mapping_audit_%s.csv", cohort))
  write_csv(audit_df, out_path)
  h    <- sha256_file(out_path)
  size <- file.info(out_path)$size / 1024^2
  registry_append(cohort, "Gene_Mapping_Audit", out_path, h,
                  "INTEGRO", script_name, size)
  message(sprintf("  [gene_map] Audit salvo: %s | %.2f MB", out_path, size))
  invisible(out_path)
}

message("[03_utils] Utilitarios de mapeamento genico carregados.")
