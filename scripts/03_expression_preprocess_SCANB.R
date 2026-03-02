# =============================================================================
# SCRIPT: 03_expression_preprocess_SCANB.R
# PURPOSE: Expression preprocessing for SCAN-B (training cohort, all 3069 samples).
#          RNA-seq: raw counts → edgeR TMM → logCPM; HGNC mapping.
# PROJETO: Core-PAM (Memorial v6.1 §4.2)
#
# INPUTS:
#   01_Base_Pura_CorePAM/RAW/GSE96058/  (GEO supplementary files)
#   01_Base_Pura_CorePAM/PROCESSED/SCANB/clinical_FINAL.parquet  (training IDs)
#
# OUTPUTS:
#   01_Base_Pura_CorePAM/PROCESSED/SCANB/expression_genelevel_preZ.parquet
#   results/supp/gene_mapping_audit_SCANB.csv
#
# PIPELINE (frozen):
#   raw counts → filter SCANB samples → edgeR TMM → logCPM → HGNC map →
#   collapse by variance → expression_genelevel_preZ (genes x samples)
#
# NOTE: Intra-cohort Z-score is done in 06_zscore_and_score_SCANB.R.
# =============================================================================

source("scripts/00_setup.R")
source("scripts/03_utils_gene_mapping.R")

# NOTE: edgeR is loaded lazily inside the normalization block (section 5)
# to avoid conflicts with Bioconductor method group bindings when
# org.Hs.eg.db/AnnotationDbi is also loaded in the same session.

SCRIPT_NAME <- "03_expression_preprocess_SCANB.R"
COHORT      <- "SCANB"

# =============================================================================
# 1) LOAD SCANB IDs (from clinical_FINAL to filter training samples)
# =============================================================================
clin_path  <- file.path(proc_cohort(COHORT), "clinical_FINAL.parquet")
clin       <- strict_parquet(clin_path)
# GSE96058 expression matrix uses SCAN-B internal IDs (patient_id = "title", e.g. "F1")
# NOT the GEO accession (sample_id = "geo_accession", e.g. "GSM2528079")
scanb_ids  <- normalize_id(clin$patient_id)
sample_lookup <- setNames(normalize_id(clin$sample_id), scanb_ids)   # patientId → sampleId
message(sprintf("[03_SCANB] Training IDs loaded: %d samples (matching by patient_id)", length(scanb_ids)))

# =============================================================================
# 2) LOCATE RAW EXPRESSION FILE IN DATA LAKE
#    GSE96058 provides several supplementary formats.
#    Priority: (1) raw counts file; (2) filtered gene matrix.
# =============================================================================
raw_dir   <- raw_cohort("GSE96058")
raw_files <- list.files(raw_dir, full.names = TRUE, recursive = FALSE)
message("[03_SCANB] Files available in RAW/GSE96058:")
print(basename(raw_files))

# Detect counts file (priority) vs gene expression matrix
# Exclude .rds files (clinical metadata) from counts detection
non_rds <- raw_files[!grepl("\\.rds$", raw_files, ignore.case = TRUE)]
counts_file <- non_rds[grep("count|RAW\\.tar", basename(non_rds),
                             ignore.case = TRUE)][1]
expr_file   <- non_rds[grep("gene_expression|genematrix", basename(non_rds),
                             ignore.case = TRUE)][1]

# =============================================================================
# 3A) READ RAW COUNTS (if available — preferred by Memorial)
# =============================================================================
if (!is.na(counts_file) && file.exists(counts_file)) {
  message("[03_SCANB] Usando arquivo de counts brutos: ", basename(counts_file))

  if (grepl("\\.tar$", counts_file)) {
    # GSE96058_RAW.tar → extract to temp and combine individual counts
    tmp_dir <- file.path(raw_dir, "RAW_extracted")
    dir.create(tmp_dir, showWarnings = FALSE)
    old_warn <- getOption("warn"); options(warn = 0)
    utils::untar(counts_file, exdir = tmp_dir)
    options(warn = old_warn)

    count_files <- list.files(tmp_dir, pattern = "\\.txt(\\.gz)?$",
                               full.names = TRUE, recursive = TRUE)
    message(sprintf("[03_SCANB] %d individual count files extracted.",
                    length(count_files)))

    # Read and combine: each file = one sample (format: GeneID \t Count)
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
    count_mat  <- count_mat[!grepl("^__", rownames(count_mat)), ] # remove HTSeq QC lines

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
# 3B) FALLBACK: gene expression matrix (if counts not available)
# =============================================================================
} else if (!is.na(expr_file) && file.exists(expr_file)) {
  message("[03_SCANB] FALLBACK: using gene expression matrix (not raw counts).")
  message("  Note: TMM will not be applied; assume pre-normalized.")

  old_warn <- getOption("warn"); options(warn = 0)
  expr_df <- read_csv(expr_file, show_col_types = FALSE)
  options(warn = old_warn)
  gene_col  <- names(expr_df)[1]
  count_mat <- as.matrix(expr_df[, -1])
  rownames(count_mat) <- expr_df[[gene_col]]
  storage.mode(count_mat) <- "numeric"
  input_type <- "preprocessed_matrix"

} else {
  stop("[03_SCANB] No expression file found in RAW/GSE96058/.\n",
       "Run 01_download_raw_data.R first.")
}

message(sprintf("[03_SCANB] Raw matrix: %d genes x %d samples | type: %s",
                nrow(count_mat), ncol(count_mat), input_type))

# =============================================================================
# 4) FILTER SCANB SAMPLES (training)
# =============================================================================
# Normalize column names to match clinical sample_ids
col_ids_norm <- normalize_id(colnames(count_mat))
keep_cols    <- col_ids_norm %in% scanb_ids
n_matched    <- sum(keep_cols)

message(sprintf("[03_SCANB] Sample match: %d of %d clinical IDs found in matrix.",
                n_matched, length(scanb_ids)))

if (n_matched == 0) {
  stop("[03_SCANB] No SCANB samples found in expression matrix.\n",
       "Check whether IDs in clinical_FINAL.parquet match column names.")
}
if (n_matched < 0.9 * length(scanb_ids)) {
  warning(sprintf("[03_SCANB] Only %.1f%% of clinical samples found in matrix.",
                  100 * n_matched / length(scanb_ids)))
}

count_mat <- count_mat[, keep_cols, drop = FALSE]
colnames(count_mat) <- col_ids_norm[keep_cols]

message(sprintf("[03_SCANB] Filtered matrix: %d genes x %d SCANB samples",
                nrow(count_mat), ncol(count_mat)))

# =============================================================================
# 5) RNA-seq NORMALIZATION: edgeR TMM → logCPM (intra-cohort)
# =============================================================================
if (input_type == "raw_counts") {
  message("[03_SCANB] Applying edgeR TMM → logCPM...")
  old_warn_pkg <- getOption("warn"); options(warn = 0)
  suppressPackageStartupMessages(library(edgeR))
  options(warn = old_warn_pkg)

  # Filter genes with very low expression (edgeR recommendation)
  dge      <- edgeR::DGEList(counts = count_mat)
  keep     <- edgeR::filterByExpr(dge)
  dge      <- dge[keep, , keep.lib.sizes = FALSE]
  message(sprintf("[03_SCANB] Genes after edgeR filter: %d (removed: %d)",
                  sum(keep), sum(!keep)))

  dge      <- edgeR::calcNormFactors(dge, method = "TMM")
  expr_mat <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

} else {
  # Already normalized — use as-is (assume log2)
  expr_mat <- count_mat
  message("[03_SCANB] Pre-normalized matrix used directly (no TMM).")
}

message(sprintf("[03_SCANB] logCPM matrix: %d genes x %d samples",
                nrow(expr_mat), ncol(expr_mat)))

# =============================================================================
# 6) HGNC MAPPING
#    Detect ID type (Ensembl vs Symbol) and map
# =============================================================================
gene_ids_raw <- rownames(expr_mat)
is_ensembl   <- grepl("^ENSG", gene_ids_raw[1])

if (is_ensembl) {
  message("[03_SCANB] IDs detected: Ensembl → mapping to HGNC...")
  hgnc_map <- build_hgnc_map_ensembl(gene_ids_raw)

  # Resolve duplicates: keep 1:1 mapping (remove Ensembl without symbol)
  hgnc_map <- hgnc_map[hgnc_map$hgnc_symbol != "" & !is.na(hgnc_map$hgnc_symbol), ]
  hgnc_map <- hgnc_map[!duplicated(hgnc_map$ensembl_gene_id), ]

  gene_symbols <- hgnc_map$hgnc_symbol[match(gene_ids_raw, hgnc_map$ensembl_gene_id)]

} else {
  message("[03_SCANB] IDs detected: Gene symbols → standardizing via HGNChelper...")
  hgnc_res     <- build_hgnc_map_symbol(gene_ids_raw)
  gene_symbols <- hgnc_res$hgnc_symbol
}

n_mapped   <- sum(!is.na(gene_symbols) & gene_symbols != "")
n_unmapped <- sum(is.na(gene_symbols)  | gene_symbols == "")
message(sprintf("[03_SCANB] Mapping: %d mapped | %d unmapped (discarded)",
                n_mapped, n_unmapped))

# Keep only mapped genes
keep_genes   <- !is.na(gene_symbols) & gene_symbols != ""
expr_mat     <- expr_mat[keep_genes, , drop = FALSE]
gene_symbols <- gene_symbols[keep_genes]

# =============================================================================
# 7) PROBE COLLAPSE: highest variance probe per gene (§4.4)
# =============================================================================
n_before <- nrow(expr_mat)
expr_mat  <- collapse_probes_by_variance(expr_mat, gene_symbols)
n_after   <- nrow(expr_mat)
message(sprintf("[03_SCANB] Probe collapse: %d → %d unique genes (HGNC)",
                n_before, n_after))

# =============================================================================
# 8) MAPPING AUDIT
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

message(sprintf("[03_SCANB] PAM50: %d/%d genes present | missing: %s",
                length(pam50_present), length(pam50_genes),
                if (length(pam50_missing) == 0) "none"
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
# 9) SAVE expression_genelevel_preZ.parquet
#    Format: rows = genes (rownames), columns = samples + "gene" column
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
message(sprintf("[03_SCANB] Saved: %s | %d genes x %d samples | %.2f MB | SHA256: %s",
                out_path, nrow(expr_mat), ncol(expr_mat), size, h))

message("\n[03_SCANB] Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message("Next step: scripts/03_expression_preprocess_TCGA_BRCA.R")
