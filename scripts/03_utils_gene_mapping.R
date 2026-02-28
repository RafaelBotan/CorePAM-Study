# =============================================================================
# UTILITY: 03_utils_gene_mapping.R
# PURPOSE: Shared gene mapping (HGNC) and probe collapse functions for all
#          03_expression_preprocess_<COHORT>.R scripts.
# PROJETO: Core-PAM (Memorial v6.1 §4.4)
#
# RULES:
#   - IDs (Symbol/Entrez/RefSeq/Ensembl) → HGNC official symbol.
#   - Multiple probes → gene: select probe with highest intra-cohort variance.
#   - Save gene_mapping_audit_<COHORT>.csv per cohort.
#   - Do NOT Z-score here — output is pre-Z (logCPM or log2 intensity).
# =============================================================================

suppressPackageStartupMessages({
  library(biomaRt)
  library(HGNChelper)
})

# --------------------------------------------------------------------------
# build_hgnc_map_ensembl
#   Converts Ensembl IDs to HGNC official symbols via biomaRt.
#   Cache: saves in registry to avoid re-download.
# --------------------------------------------------------------------------
build_hgnc_map_ensembl <- function(ensembl_ids,
                                   cache_path = file.path(
                                     PATHS$registry_docs,
                                     "hgnc_ensembl_map_cache.rds")) {
  if (file.exists(cache_path)) {
    message("  [gene_map] Using HGNC Ensembl cache: ", cache_path)
    old_warn <- getOption("warn"); options(warn = 0)
    map <- readRDS(cache_path)
    options(warn = old_warn)
  } else {
    message("  [gene_map] Downloading Ensembl→HGNC map via biomaRt...")
    old_warn <- getOption("warn"); options(warn = 0)
    mart <- tryCatch(
      biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl"),
      error = function(e) {
        options(warn = old_warn)
        stop("biomaRt unavailable: ", conditionMessage(e),
             "\nCheck connection or use manual cache.")
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
    message("  [gene_map] Cache saved: ", cache_path)
  }
  map
}

# --------------------------------------------------------------------------
# build_hgnc_map_symbol
#   Standardizes existing gene symbols using HGNChelper (offline).
#   Returns data.frame with original → hgnc_symbol.
# --------------------------------------------------------------------------
build_hgnc_map_symbol <- function(symbols) {
  message("  [gene_map] Standardizing symbols via HGNChelper...")
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
#   Converts Entrez IDs to HGNC symbols via org.Hs.eg.db.
# --------------------------------------------------------------------------
build_hgnc_map_entrez <- function(entrez_ids) {
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("Package org.Hs.eg.db required: ",
         "BiocManager::install('org.Hs.eg.db')")
  }
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  message("  [gene_map] Mapping Entrez→HGNC via org.Hs.eg.db...")
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
#   Converts Affymetrix HGU133A/Plus2 probe IDs to HGNC symbols.
# --------------------------------------------------------------------------
build_hgnc_map_affymetrix <- function(probe_ids, platform = "hgu133a") {
  pkg <- paste0(platform, ".db")
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package %s required: BiocManager::install('%s')", pkg, pkg))
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  message(sprintf("  [gene_map] Mapping probes %s→SYMBOL via %s...",
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
#   For matrices with multiple rows per gene:
#   selects the row (probe) with highest intra-cohort variance.
#   Input:  mat     = matrix/data.frame (rows=probes, columns=samples)
#           gene_id = vector of gene symbols (same length as nrow(mat))
#   Output: matrix with one row per gene (highest variance probe)
# --------------------------------------------------------------------------
collapse_probes_by_variance <- function(mat, gene_id) {
  stopifnot(length(gene_id) == nrow(mat))

  df <- data.frame(
    gene     = gene_id,
    row_idx  = seq_len(nrow(mat)),
    variance = apply(mat, 1, var, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  # Remove NA or empty genes
  df <- df[!is.na(df$gene) & nchar(trimws(df$gene)) > 0, ]

  # Select highest variance probe per gene
  best <- df[order(df$gene, -df$variance), ]
  best <- best[!duplicated(best$gene), ]

  mat_collapsed <- mat[best$row_idx, , drop = FALSE]
  rownames(mat_collapsed) <- best$gene
  mat_collapsed
}

# --------------------------------------------------------------------------
# save_gene_mapping_audit
#   Saves gene_mapping_audit_<cohort>.csv and registers hash.
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
  message(sprintf("  [gene_map] Audit saved: %s | %.2f MB", out_path, size))
  invisible(out_path)
}

message("[03_utils] Gene mapping utilities loaded.")
