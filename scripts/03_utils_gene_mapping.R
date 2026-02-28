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

old_warn <- getOption("warn"); options(warn = 0)
suppressPackageStartupMessages({
  library(HGNChelper)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})
options(warn = old_warn)

# --------------------------------------------------------------------------
# build_hgnc_map_ensembl
#   Converts Ensembl IDs to HGNC official symbols via org.Hs.eg.db (offline).
#   Replaces biomaRt (which had load failures on this R installation).
#   Cache: saves in registry to avoid re-query.
# --------------------------------------------------------------------------
build_hgnc_map_ensembl <- function(ensembl_ids,
                                   cache_path = file.path(
                                     PATHS$registry_docs,
                                     "hgnc_ensembl_map_cache.rds")) {
  if (file.exists(cache_path)) {
    message("  [gene_map] Using HGNC Ensembl cache: ", cache_path)
    old_warn <- getOption("warn"); on.exit(options(warn = old_warn))
    options(warn = 0)
    map <- readRDS(cache_path)
    return(map)
  }

  message("  [gene_map] Building Ensembl→HGNC map via org.Hs.eg.db (offline)...")
  old_warn <- getOption("warn"); on.exit(options(warn = old_warn))
  options(warn = 0)

  clean_ids <- unique(sub("\\..*$", "", ensembl_ids))   # strip version suffix

  map_raw <- tryCatch(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys    = clean_ids[!is.na(clean_ids) & nchar(clean_ids) > 0],
      columns = c("ENSEMBL", "SYMBOL", "ENTREZID", "GENETYPE"),
      keytype = "ENSEMBL"
    ),
    error = function(e) {
      message("  [gene_map] org.Hs.eg.db mapping failed: ", conditionMessage(e))
      data.frame(ENSEMBL = character(), SYMBOL = character(),
                 ENTREZID = character(), GENETYPE = character(),
                 stringsAsFactors = FALSE)
    }
  )

  # Rename columns to match old biomaRt output format used downstream
  map <- data.frame(
    ensembl_gene_id = map_raw$ENSEMBL,
    hgnc_symbol     = map_raw$SYMBOL,
    entrezgene_id   = suppressWarnings(as.integer(map_raw$ENTREZID)),
    gene_biotype    = map_raw$GENETYPE,
    stringsAsFactors = FALSE
  )

  dir.create(dirname(cache_path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(map, cache_path)
  message(sprintf("  [gene_map] Cache saved: %s (%d mappings)", cache_path, nrow(map)))
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
