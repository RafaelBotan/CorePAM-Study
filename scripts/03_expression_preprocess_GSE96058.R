# =============================================================================
# SCRIPT: 03_expression_preprocess_GSE96058.R
# PURPOSE: Expression preprocessing for GSE96058 (validation subset).
#          Same raw source as SCANB, but different samples — separate TMM.
# PROJETO: Core-PAM (Memorial v6.1 §4.2)
#
# INPUTS:
#   01_Base_Pura_CorePAM/RAW/GSE96058/  (same raw as SCANB)
#   01_Base_Pura_CorePAM/PROCESSED/GSE96058/clinical_FINAL.parquet (validation IDs)
#
# OUTPUTS:
#   01_Base_Pura_CorePAM/PROCESSED/GSE96058/expression_genelevel_preZ.parquet
#   results/supp/gene_mapping_audit_GSE96058.csv
#
# CRITICAL RULE (Memorial §1.2 / §3.2):
#   TMM calculated ONLY on GSE96058 samples (validation).
#   NEVER combine with SCANB for normalization.
#   Intra-cohort Z-score done in 06_zscore_and_score_GSE96058.R.
# =============================================================================

source("scripts/00_setup.R")
source("scripts/03_utils_gene_mapping.R")

suppressPackageStartupMessages(library(edgeR))

SCRIPT_NAME <- "03_expression_preprocess_GSE96058.R"
COHORT      <- "GSE96058"

# =============================================================================
# 1) LOAD GSE96058 IDs (validation)
# =============================================================================
clin_path <- file.path(proc_cohort(COHORT), "clinical_FINAL.parquet")
clin      <- strict_parquet(clin_path)
valid_ids  <- normalize_id(clin$sample_id)
message(sprintf("[03_GSE96058] Validation IDs loaded: %d samples", length(valid_ids)))

# =============================================================================
# 2) LOCATE EXPRESSION FILE (same directory as SCANB)
# =============================================================================
raw_dir   <- raw_cohort("GSE96058")
raw_files <- list.files(raw_dir, full.names = TRUE, recursive = FALSE)

counts_file <- raw_files[grep("count|raw|RAW\\.tar", basename(raw_files),
                               ignore.case = TRUE)][1]
expr_file   <- raw_files[grep("gene_expression|genematrix|expression",
                               basename(raw_files), ignore.case = TRUE)][1]

# =============================================================================
# 3) READ (identical to SCANB — reuses same parsing logic)
# =============================================================================
if (!is.na(counts_file) && file.exists(counts_file)) {
  message("[03_GSE96058] Using raw counts: ", basename(counts_file))

  if (grepl("\\.tar$", counts_file)) {
    tmp_dir     <- file.path(raw_dir, "RAW_extracted")
    count_files <- list.files(tmp_dir, pattern = "\\.txt(\\.gz)?$",
                               full.names = TRUE, recursive = TRUE)
    if (length(count_files) == 0) {
      dir.create(tmp_dir, showWarnings = FALSE)
      old_warn <- getOption("warn"); options(warn = 0)
      utils::untar(counts_file, exdir = tmp_dir)
      options(warn = old_warn)
      count_files <- list.files(tmp_dir, pattern = "\\.txt(\\.gz)?$",
                                 full.names = TRUE, recursive = TRUE)
    }
    counts_list  <- lapply(count_files, function(f) {
      old_warn <- getOption("warn"); options(warn = 0)
      d <- read_tsv(f, col_names = c("gene_id", "count"),
                    show_col_types = FALSE)
      options(warn = old_warn)
      setNames(d$count, d$gene_id)
    })
    sample_names <- sub("(_count|_counts|\\.txt.*)", "",
                        basename(count_files), ignore.case = TRUE)
    all_genes    <- Reduce(intersect, lapply(counts_list, names))
    count_mat    <- do.call(cbind, lapply(counts_list, function(x) x[all_genes]))
    colnames(count_mat) <- sample_names
    count_mat    <- count_mat[!grepl("^__", rownames(count_mat)), ]
    input_type   <- "raw_counts"

  } else {
    old_warn <- getOption("warn"); options(warn = 0)
    count_df  <- read_csv(counts_file, show_col_types = FALSE)
    options(warn = old_warn)
    gene_col  <- names(count_df)[1]
    count_mat <- as.matrix(count_df[, -1])
    rownames(count_mat) <- count_df[[gene_col]]
    input_type <- "raw_counts"
  }

} else if (!is.na(expr_file) && file.exists(expr_file)) {
  message("[03_GSE96058] FALLBACK: pre-normalized gene expression matrix.")
  old_warn <- getOption("warn"); options(warn = 0)
  expr_df  <- read_csv(expr_file, show_col_types = FALSE)
  options(warn = old_warn)
  gene_col  <- names(expr_df)[1]
  count_mat <- as.matrix(expr_df[, -1])
  rownames(count_mat) <- expr_df[[gene_col]]
  storage.mode(count_mat) <- "numeric"
  input_type <- "preprocessed_matrix"

} else {
  stop("[03_GSE96058] No expression file found in RAW/GSE96058/.")
}

# =============================================================================
# 4) FILTER ONLY GSE96058 SAMPLES (validation) — SEPARATE TMM FROM SCANB
# =============================================================================
col_ids_norm <- normalize_id(colnames(count_mat))
keep_cols    <- col_ids_norm %in% valid_ids
n_matched    <- sum(keep_cols)

message(sprintf("[03_GSE96058] Match: %d of %d validation IDs found.",
                n_matched, length(valid_ids)))

if (n_matched == 0) {
  stop("[03_GSE96058] No GSE96058 samples found. ",
       "Check IDs in clinical_FINAL.parquet vs column names.")
}

count_mat <- count_mat[, keep_cols, drop = FALSE]
colnames(count_mat) <- col_ids_norm[keep_cols]

message(sprintf("[03_GSE96058] Filtered matrix: %d genes x %d GSE96058 samples",
                nrow(count_mat), ncol(count_mat)))

# =============================================================================
# 5) NORMALIZATION (intra-cohort GSE96058 ONLY)
# =============================================================================
if (input_type == "raw_counts") {
  message("[03_GSE96058] Applying edgeR TMM → logCPM (intra-GSE96058)...")
  dge  <- edgeR::DGEList(counts = count_mat)
  keep <- edgeR::filterByExpr(dge)
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  message(sprintf("[03_GSE96058] Genes after edgeR filter: %d", sum(keep)))
  dge      <- edgeR::calcNormFactors(dge, method = "TMM")
  expr_mat <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
} else {
  expr_mat <- count_mat
  message("[03_GSE96058] Pre-normalized matrix used directly.")
}

# =============================================================================
# 6) HGNC MAPPING
# =============================================================================
gene_ids_raw <- rownames(expr_mat)
is_ensembl   <- grepl("^ENSG", gene_ids_raw[1])

if (is_ensembl) {
  hgnc_map     <- build_hgnc_map_ensembl(gene_ids_raw)
  hgnc_map     <- hgnc_map[hgnc_map$hgnc_symbol != "" & !is.na(hgnc_map$hgnc_symbol), ]
  hgnc_map     <- hgnc_map[!duplicated(hgnc_map$ensembl_gene_id), ]
  gene_symbols <- hgnc_map$hgnc_symbol[match(gene_ids_raw, hgnc_map$ensembl_gene_id)]
} else {
  hgnc_res     <- build_hgnc_map_symbol(gene_ids_raw)
  gene_symbols <- hgnc_res$hgnc_symbol
}

keep_genes   <- !is.na(gene_symbols) & gene_symbols != ""
expr_mat     <- expr_mat[keep_genes, , drop = FALSE]
gene_symbols <- gene_symbols[keep_genes]
message(sprintf("[03_GSE96058] Mapped: %d | Discarded: %d",
                sum(keep_genes), sum(!keep_genes)))

# =============================================================================
# 7) PROBE COLLAPSE
# =============================================================================
n_before <- nrow(expr_mat)
expr_mat  <- collapse_probes_by_variance(expr_mat, gene_symbols)
message(sprintf("[03_GSE96058] Collapse: %d → %d unique genes", n_before, nrow(expr_mat)))

# =============================================================================
# 8) PAM50 AUDIT
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
message(sprintf("[03_GSE96058] PAM50: %d/50 present", length(pam50_present)))

audit_df <- tibble(cohort = COHORT, gene = rownames(expr_mat),
                   present = TRUE, in_pam50 = rownames(expr_mat) %in% pam50_genes) |>
  bind_rows(tibble(cohort = COHORT, gene = pam50_missing,
                   present = FALSE, in_pam50 = TRUE))
save_gene_mapping_audit(audit_df, COHORT, SCRIPT_NAME)

# =============================================================================
# 9) SAVE
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
message(sprintf("[03_GSE96058] Saved: %d genes x %d samples | %.2f MB | SHA256: %s",
                nrow(expr_mat), ncol(expr_mat), size, h))

message("\n[03_GSE96058] Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message("Next step: scripts/03_expression_preprocess_TCGA_BRCA.R")
