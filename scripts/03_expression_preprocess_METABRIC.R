# =============================================================================
# SCRIPT: 03_expression_preprocess_METABRIC.R
# PURPOSE: Expression preprocessing for METABRIC.
#          Illumina HT-12 v3 microarray: log2 intensity as-is → HGNC map.
# PROJETO: Core-PAM (Memorial v6.1 §4.3)
#
# INPUT:
#   01_Base_Pura_CorePAM/RAW/METABRIC/brca_metabric/data_mrna_illumina_microarray.txt
#   01_Base_Pura_CorePAM/PROCESSED/METABRIC/clinical_FINAL.parquet
#
# OUTPUTS:
#   01_Base_Pura_CorePAM/PROCESSED/METABRIC/expression_genelevel_preZ.parquet
#   results/supp/gene_mapping_audit_METABRIC.csv
#
# RULES (Memorial v6.1 §4.3):
#   - Microarray: log2-intensity "as-is" (no TMM).
#   - Intra-cohort Z-score standardization in 06_zscore_and_score_METABRIC.R.
#   - Multiple probes per gene → highest intra-cohort variance.
#   - HGNC mapping: cBioPortal file already has Hugo_Symbol → validate via HGNChelper.
# =============================================================================

source("scripts/00_setup.R")
source("scripts/03_utils_gene_mapping.R")

SCRIPT_NAME <- "03_expression_preprocess_METABRIC.R"
COHORT      <- "METABRIC"

# =============================================================================
# 1) READ MATRIX (cBioPortal format)
#    Format: Hugo_Symbol | Entrez_Gene_Id | SAMPLE1 | SAMPLE2 | ...
#    Lines starting with "#" are metadata → skip
# =============================================================================
expr_path <- file.path(raw_cohort(COHORT), "brca_metabric",
                       "data_mrna_illumina_microarray.txt")
if (!file.exists(expr_path)) {
  stop("[03_METABRIC] File not found: ", expr_path,
       "\nRun 01_download_raw_data.R first.")
}

message("[03_METABRIC] Reading cBioPortal expression matrix...")
lines  <- readLines(expr_path, warn = FALSE)
skip_n <- sum(startsWith(lines, "#"))

old_warn <- getOption("warn"); options(warn = 0)
expr_raw <- read_tsv(expr_path, skip = skip_n, show_col_types = FALSE,
                     name_repair = "unique")
options(warn = old_warn)

message(sprintf("[03_METABRIC] Read: %d rows x %d columns (including Hugo_Symbol, Entrez)",
                nrow(expr_raw), ncol(expr_raw)))

# Separate gene and sample columns
GENE_COL   <- "Hugo_Symbol"
ENTREZ_COL <- "Entrez_Gene_Id"
sample_cols <- setdiff(names(expr_raw), c(GENE_COL, ENTREZ_COL))

gene_symbols_raw <- expr_raw[[GENE_COL]]
entrez_ids_raw   <- expr_raw[[ENTREZ_COL]]

expr_mat <- as.matrix(expr_raw[, sample_cols])
rownames(expr_mat) <- gene_symbols_raw
storage.mode(expr_mat) <- "numeric"

message(sprintf("[03_METABRIC] Matrix: %d genes x %d samples",
                nrow(expr_mat), ncol(expr_mat)))

# =============================================================================
# 2) FILTER SAMPLES WITH CLINICAL DATA
# =============================================================================
clin_path <- file.path(proc_cohort(COHORT), "clinical_FINAL.parquet")
clin      <- strict_parquet(clin_path)
clin_ids  <- normalize_id(clin$sample_id)

col_ids_norm <- normalize_id(colnames(expr_mat))
keep_cols    <- col_ids_norm %in% clin_ids
n_matched    <- sum(keep_cols)

message(sprintf("[03_METABRIC] Sample match: %d of %d clinical IDs found.",
                n_matched, length(clin_ids)))

if (n_matched == 0) {
  stop("[03_METABRIC] No clinical samples found in matrix.\n",
       "Check IDs in clinical_FINAL.parquet vs sample columns.")
}

expr_mat <- expr_mat[, keep_cols, drop = FALSE]
colnames(expr_mat) <- col_ids_norm[keep_cols]

message(sprintf("[03_METABRIC] Filtered matrix: %d genes x %d samples",
                nrow(expr_mat), ncol(expr_mat)))

# =============================================================================
# 3) CHECK LOG2 SCALE
#    Illumina HT-12 cBioPortal: typically already in log2.
#    Heuristic: if 95th percentile < 20 → assume log2; else → apply log2.
# =============================================================================
p95 <- quantile(expr_mat, 0.95, na.rm = TRUE)
if (p95 > 20) {
  message(sprintf("[03_METABRIC] p95=%.1f > 20 → applying log2(x + 1).", p95))
  expr_mat <- log2(expr_mat + 1)
} else {
  message(sprintf("[03_METABRIC] p95=%.1f <= 20 → data already in log2 (as-is).", p95))
}

# =============================================================================
# 4) HGNC MAPPING (Hugo_Symbol already available — validate via HGNChelper)
# =============================================================================
message("[03_METABRIC] Validating Hugo Symbols via HGNChelper...")
hgnc_res     <- build_hgnc_map_symbol(rownames(expr_mat))
gene_symbols <- hgnc_res$hgnc_symbol

# Remove genes without valid HGNC mapping
keep_genes   <- !is.na(gene_symbols) & gene_symbols != ""
n_unmapped   <- sum(!keep_genes)
expr_mat     <- expr_mat[keep_genes, , drop = FALSE]
gene_symbols <- gene_symbols[keep_genes]

message(sprintf("[03_METABRIC] HGNC: %d valid | %d discarded",
                sum(keep_genes), n_unmapped))

# =============================================================================
# 5) PROBE COLLAPSE (highest intra-cohort variance)
# =============================================================================
n_before <- nrow(expr_mat)
expr_mat  <- collapse_probes_by_variance(expr_mat, gene_symbols)
n_after   <- nrow(expr_mat)
message(sprintf("[03_METABRIC] Probe collapse: %d → %d unique HGNC genes",
                n_before, n_after))

# =============================================================================
# 6) QC: remove genes with missingness > 20%
# =============================================================================
miss_frac <- rowMeans(is.na(expr_mat))
n_high_miss <- sum(miss_frac > 0.20)
if (n_high_miss > 0) {
  message(sprintf("[03_METABRIC] Removing %d genes with >20%% missingness.",
                  n_high_miss))
  expr_mat <- expr_mat[miss_frac <= 0.20, , drop = FALSE]
}

# =============================================================================
# 7) AUDIT PAM50
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
message(sprintf("[03_METABRIC] PAM50: %d/50 present | missing: %s",
                length(pam50_present),
                if (length(pam50_missing) == 0) "none"
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
# 8) SAVE
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
message(sprintf("[03_METABRIC] Saved: %d genes x %d samples | %.2f MB | SHA256: %s",
                nrow(expr_mat), ncol(expr_mat), size, h))

message("\n[03_METABRIC] Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message("Next step: scripts/03_expression_preprocess_GSE20685.R")
