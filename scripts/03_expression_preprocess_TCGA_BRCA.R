# =============================================================================
# SCRIPT: 03_expression_preprocess_TCGA_BRCA.R
# PURPOSE: Expression preprocessing for TCGA-BRCA.
#          RNA-seq STAR counts → edgeR TMM → logCPM; HGNC mapping.
# PROJETO: Core-PAM (Memorial v6.1 §4.2)
#
# INPUTS:
#   01_Base_Pura_CorePAM/RAW/TCGA_BRCA/TCGA_BRCA_SE_counts_raw.rds
#   01_Base_Pura_CorePAM/PROCESSED/TCGA_BRCA/clinical_FINAL.parquet
#
# OUTPUTS:
#   01_Base_Pura_CorePAM/PROCESSED/TCGA_BRCA/expression_genelevel_preZ.parquet
#   results/supp/gene_mapping_audit_TCGA_BRCA.csv
#
# NOTE: TCGA IDs are Ensembl (ENSG...) with version (.X) — remove version.
#       Intra-cohort Z-score done in 06_zscore_and_score_TCGA_BRCA.R.
# =============================================================================

source("scripts/00_setup.R")
# Load edgeR and SummarizedExperiment BEFORE 03_utils_gene_mapping.R
# to avoid Bioconductor locked binding conflicts (edgeR must load before AnnotationDbi)
old_warn_pkg <- getOption("warn"); options(warn = 0)
suppressPackageStartupMessages({
  library(edgeR)
  library(SummarizedExperiment)
})
options(warn = old_warn_pkg)

source("scripts/03_utils_gene_mapping.R")

SCRIPT_NAME <- "03_expression_preprocess_TCGA_BRCA.R"
COHORT      <- "TCGA_BRCA"

# =============================================================================
# 1) LOAD SummarizedExperiment
# =============================================================================
se_path <- file.path(raw_cohort(COHORT), "TCGA_BRCA_SE_counts_raw.rds")
se_obj  <- strict_rds(se_path)
message(sprintf("[03_TCGA] SE loaded: %d genes x %d samples",
                nrow(se_obj), ncol(se_obj)))
message("[03_TCGA] Available assays: ", paste(assayNames(se_obj), collapse = ", "))

# =============================================================================
# 2) EXTRACT RAW COUNTS (unstranded = GDC STAR standard)
# =============================================================================
count_assay <- if ("unstranded" %in% assayNames(se_obj)) "unstranded" else assayNames(se_obj)[1]
message("[03_TCGA] Usando assay: ", count_assay)
count_mat <- assay(se_obj, count_assay)

# =============================================================================
# 3) FILTER Primary Tumor only + samples with clinical data
# =============================================================================
clin_path <- file.path(proc_cohort(COHORT), "clinical_FINAL.parquet")
clin      <- strict_parquet(clin_path)
pat_ids   <- normalize_id(clin$patient_id)

# TCGA barcode: first 12 characters = patient_id; positions 14-15 = type
barcodes      <- colnames(count_mat)
sample_type   <- substr(barcodes, 14, 15)
patient_from_bc <- normalize_id(substr(barcodes, 1, 12))

keep_cols <- sample_type == "01" & patient_from_bc %in% pat_ids
message(sprintf("[03_TCGA] Primary Tumor with clinical data: %d of %d samples",
                sum(keep_cols), length(barcodes)))

count_mat <- count_mat[, keep_cols, drop = FALSE]
# Rename columns to patient_id (12 chars) — one per patient
colnames(count_mat) <- normalize_id(substr(colnames(count_mat), 1, 12))

# Remove patient duplicates (multiple tubes) — keep first occurrence
dup_mask  <- duplicated(colnames(count_mat))
if (any(dup_mask)) {
  message(sprintf("[03_TCGA] Removing %d duplicate barcodes (same patient_id).",
                  sum(dup_mask)))
  count_mat <- count_mat[, !dup_mask, drop = FALSE]
}
message(sprintf("[03_TCGA] Filtered matrix: %d genes x %d patients",
                nrow(count_mat), ncol(count_mat)))

# =============================================================================
# 4) REMOVE ENSEMBL ID VERSION (ENSG00000XXX.Y → ENSG00000XXX)
# =============================================================================
gene_ids_raw <- rownames(count_mat)
gene_ids_clean <- sub("\\.\\d+$", "", gene_ids_raw)
rownames(count_mat) <- gene_ids_clean

# Remove genes with non-Ensembl IDs (PAR_Y, etc.)
keep_genes <- grepl("^ENSG", gene_ids_clean)
count_mat  <- count_mat[keep_genes, , drop = FALSE]
message(sprintf("[03_TCGA] After removing non-Ensembl: %d genes", nrow(count_mat)))

# =============================================================================
# 5) NORMALIZATION: edgeR TMM → logCPM (intra-TCGA)
# =============================================================================
message("[03_TCGA] Applying edgeR TMM → logCPM...")
dge  <- edgeR::DGEList(counts = count_mat)
old_warn_filter <- getOption("warn"); options(warn = 0)
keep <- edgeR::filterByExpr(dge)
options(warn = old_warn_filter)
dge  <- dge[keep, , keep.lib.sizes = FALSE]
message(sprintf("[03_TCGA] Genes after edgeR filter: %d (removed: %d)",
                sum(keep), sum(!keep)))
dge      <- edgeR::calcNormFactors(dge, method = "TMM")
expr_mat <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
message(sprintf("[03_TCGA] logCPM: %d genes x %d samples", nrow(expr_mat), ncol(expr_mat)))

# =============================================================================
# 6) ENSEMBL → HGNC MAPPING
# =============================================================================
message("[03_TCGA] Mapping Ensembl → HGNC official symbols...")
hgnc_map     <- build_hgnc_map_ensembl(rownames(expr_mat))
hgnc_map     <- hgnc_map[hgnc_map$hgnc_symbol != "" & !is.na(hgnc_map$hgnc_symbol), ]
hgnc_map     <- hgnc_map[!duplicated(hgnc_map$ensembl_gene_id), ]
gene_symbols <- hgnc_map$hgnc_symbol[match(rownames(expr_mat), hgnc_map$ensembl_gene_id)]

n_mapped   <- sum(!is.na(gene_symbols) & gene_symbols != "")
n_unmapped <- sum(is.na(gene_symbols)  | gene_symbols == "")
message(sprintf("[03_TCGA] Mapped: %d | Unmapped: %d", n_mapped, n_unmapped))

keep_g       <- !is.na(gene_symbols) & gene_symbols != ""
expr_mat     <- expr_mat[keep_g, , drop = FALSE]
gene_symbols <- gene_symbols[keep_g]

# =============================================================================
# 7) COLLAPSE (one gene per HGNC symbol)
# =============================================================================
n_before <- nrow(expr_mat)
expr_mat  <- collapse_probes_by_variance(expr_mat, gene_symbols)
message(sprintf("[03_TCGA] Collapse: %d → %d unique HGNC genes",
                n_before, nrow(expr_mat)))

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
message(sprintf("[03_TCGA] PAM50: %d/50 present | missing: %s",
                length(pam50_present),
                if (length(pam50_missing) == 0) "none"
                else paste(pam50_missing, collapse = ",")))

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
message(sprintf("[03_TCGA] Saved: %d genes x %d samples | %.2f MB | SHA256: %s",
                nrow(expr_mat), ncol(expr_mat), size, h))

message("\n[03_TCGA] Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message("Next step: scripts/03_expression_preprocess_METABRIC.R")
