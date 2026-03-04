# =============================================================================
# SCRIPT: 20_prepare_pCR_ISPY2.R
# PURPOSE: Harmonize clinical (pCR + covariates), load batch-corrected
#          gene-level expression, compute CorePAM score for I-SPY2 (GSE194040).
#          Produces analysis_ready.parquet for 21_ispy2_analysis.R.
#
# COHORT:  ISPY2 | GEO: GSE194040 | Endpoint: pCR
# PLATFORM: Agilent 44K (GPL20078 + GPL30493) — batch-corrected by authors
#           Gene-level file already batch-corrected (ComBat by Agendia/UCSF).
#           CorePAM coverage: 24/24 genes confirmed present (19,134 genes total).
# N expected: 988 pre-treatment biopsies (both platforms combined)
# CLINICAL:  pCR (0/1), HR, HER2, treatment arm (10 arms), MammaPrint (mp1/mp2)
#            NOTE: age NOT available in GEO metadata for I-SPY2.
# ROLE:      Exploratory External Cohort (supplementary analysis only)
# REF:       Hirst GL et al. NPJ Breast Cancer 2020; NACT multi-arm trial
#
# ANALYSIS PLAN:
#   (i)   Univariate:  pCR ~ score_z
#   (ii)  Adjusted:    pCR ~ score_z + hr + her2 + arm
#   (iii) Control-only sensitivity: subset arm == "Paclitaxel"
#
# PROJECT: Core-PAM (Memorial v6.1 / Freeze Core-PAM — pCR block)
# =============================================================================

source("scripts/00_setup.R")
source("scripts/20_utils_pcr_extract.R")

SCRIPT_NAME <- "20_prepare_pCR_ISPY2.R"
COHORT      <- "ISPY2"
GEO_ACC     <- "GSE194040"
PLATFORM    <- "agilent_44k_batch_corrected"

FORCE <- as.logical(Sys.getenv("FORCE_RERUN", "FALSE"))

out_path <- file.path(proc_pcr_cohort(COHORT), "analysis_ready.parquet")
if (!FORCE && file.exists(out_path)) {
  message(sprintf("[%s] Output exists — skipping. Set FORCE_RERUN=TRUE to rerun.",
                  SCRIPT_NAME))
  quit(save = "no", status = 0)
}

message(sprintf("[%s] Starting pCR prepare for %s (%s)", SCRIPT_NAME, COHORT, GEO_ACC))

old_warn <- getOption("warn"); options(warn = 0)
suppressPackageStartupMessages(library(GEOquery))
suppressPackageStartupMessages(library(Biobase))
options(warn = old_warn)

# --------------------------------------------------------------------------
# 1) Load clinical from series matrix files (two platforms combined)
#    GPL20078 = 654 samples, GPL30493 = 334 samples → 988 total
# --------------------------------------------------------------------------
raw_dir <- raw_pcr_cohort(COHORT)
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)

# The two series matrix files should already be present from 19b download.
# If missing, download them now.
sm_files_expected <- c(
  file.path(raw_dir, sprintf("%s-GPL20078_series_matrix.txt.gz", GEO_ACC)),
  file.path(raw_dir, sprintf("%s-GPL30493_series_matrix.txt.gz", GEO_ACC))
)

if (!all(file.exists(sm_files_expected))) {
  message(sprintf("[%s] Downloading series matrix files for %s...", SCRIPT_NAME, GEO_ACC))
  old_warn2 <- getOption("warn"); options(warn = 0)
  gse_list <- tryCatch(
    GEOquery::getGEO(GEO = GEO_ACC, destdir = raw_dir,
                     GSEMatrix = TRUE, getGPL = FALSE, AnnotGPL = FALSE),
    error = function(e) { message("ERROR: ", e$message); NULL }
  )
  options(warn = old_warn2)
  if (is.null(gse_list)) stop(sprintf("[%s] Could not download %s", SCRIPT_NAME, GEO_ACC))
}
stopifnot(all(file.exists(sm_files_expected)))

message(sprintf("[%s] Reading two series matrix files...", SCRIPT_NAME))

read_pdata_only <- function(sm_path, cohort_tag) {
  old_warn <- getOption("warn"); options(warn = 0)
  gse_sm <- GEOquery::getGEO(filename = sm_path)
  options(warn = old_warn)
  if (is.list(gse_sm)) gse_sm <- gse_sm[[1]]
  pd <- Biobase::pData(gse_sm)
  message(sprintf("[%s]   %s: %d samples, %d columns",
                  SCRIPT_NAME, basename(sm_path), nrow(pd), ncol(pd)))
  pd
}

pd_gpl20078 <- read_pdata_only(sm_files_expected[1], "GPL20078")
pd_gpl30493 <- read_pdata_only(sm_files_expected[2], "GPL30493")

# Standardize column names to lowercase for both
names(pd_gpl20078) <- tolower(names(pd_gpl20078))
names(pd_gpl30493) <- tolower(names(pd_gpl30493))

# Combine using common columns (bind_rows fills missing with NA)
old_warn <- getOption("warn"); options(warn = 0)
pdata_combined <- tryCatch(
  dplyr::bind_rows(pd_gpl20078, pd_gpl30493),
  error = function(e) {
    common_cols <- intersect(names(pd_gpl20078), names(pd_gpl30493))
    rbind(pd_gpl20078[, common_cols, drop = FALSE],
          pd_gpl30493[, common_cols, drop = FALSE])
  }
)
options(warn = old_warn)
message(sprintf("[%s] Combined pData: %d samples | %d columns",
                SCRIPT_NAME, nrow(pdata_combined), ncol(pdata_combined)))

# --------------------------------------------------------------------------
# 2) Extract patient IDs from pData (used to match gene matrix columns)
#    Column: "patient id:ch1" → numeric IDs like 1001, 1002, ...
# --------------------------------------------------------------------------
# After tolower() rename, the column is named "patient id:ch1"
patient_id_col <- grep("patient.id", names(pdata_combined), ignore.case = TRUE, value = TRUE)
if (length(patient_id_col) == 0) {
  # Try GSM-based (geo_accession) as fallback
  patient_id_col <- "geo_accession"
  warning(sprintf("[%s] 'patient id:ch1' not found; using geo_accession as ID", SCRIPT_NAME))
}
patient_id_col <- patient_id_col[1]
message(sprintf("[%s] Patient ID column: '%s'", SCRIPT_NAME, patient_id_col))

raw_ids <- trimws(sub("^[^:]+:", "", as.character(pdata_combined[[patient_id_col]])))
pdata_combined$patient_id_clean <- raw_ids
message(sprintf("[%s] Sample patient IDs (first 5): %s",
                SCRIPT_NAME, paste(head(raw_ids, 5), collapse = ", ")))

# --------------------------------------------------------------------------
# 3) Extract pCR (0/1) from pData via utility
# --------------------------------------------------------------------------
pcr_result <- extract_pcr_column(pdata_combined, cohort = COHORT)
pdata_combined$pcr <- pcr_result$pcr
message(sprintf("[%s] pCR column: '%s' | pCR+=%d | pCR-=%d | NA=%d",
                SCRIPT_NAME, pcr_result$col_used,
                sum(pdata_combined$pcr == 1L, na.rm = TRUE),
                sum(pdata_combined$pcr == 0L, na.rm = TRUE),
                sum(is.na(pdata_combined$pcr))))

# --------------------------------------------------------------------------
# 4) Extract I-SPY2 specific clinical covariates
#    hr:ch1   — hormone receptor (0=neg, 1=pos)
#    her2:ch1 — HER2 status   (0=neg, 1=pos)
#    arm:ch1  — treatment arm (10 levels; control = "Paclitaxel")
#    mp:ch1   — MammaPrint risk (mp1/mp2 or similar)
# --------------------------------------------------------------------------
.geo_val <- function(pdata, pattern) {
  col <- grep(pattern, names(pdata), ignore.case = TRUE, value = TRUE)
  if (length(col) == 0) return(rep(NA_character_, nrow(pdata)))
  col <- col[1]
  trimws(sub("^[^:]+:", "", as.character(pdata[[col]])))
}

# HR status → 0/1 integer
hr_raw  <- .geo_val(pdata_combined, "^hr")
her2_raw <- .geo_val(pdata_combined, "her2")
arm_raw  <- .geo_val(pdata_combined, "^arm")
mp_raw   <- .geo_val(pdata_combined, "^mp|mammaprint")

pdata_combined$hr   <- suppressWarnings(as.integer(hr_raw))
pdata_combined$her2 <- suppressWarnings(as.integer(her2_raw))
pdata_combined$arm  <- arm_raw
pdata_combined$mp   <- mp_raw

message(sprintf("[%s] HR valid: %d | HER2 valid: %d | arm unique: %d | mp unique: %d",
                SCRIPT_NAME,
                sum(!is.na(pdata_combined$hr)),
                sum(!is.na(pdata_combined$her2)),
                length(unique(pdata_combined$arm[!is.na(pdata_combined$arm)])),
                length(unique(pdata_combined$mp[!is.na(pdata_combined$mp)]))))
message(sprintf("[%s] Arm levels: %s",
                SCRIPT_NAME,
                paste(sort(unique(pdata_combined$arm[!is.na(pdata_combined$arm)])),
                      collapse = " | ")))

# --------------------------------------------------------------------------
# 5) Load batch-corrected gene-level expression matrix
#    Format: rows = genes (rownames), columns = patient IDs (numeric)
#    Already batch-corrected across both platforms by Agendia/UCSF (ComBat).
# --------------------------------------------------------------------------
gene_file <- file.path(raw_dir, sprintf("%s_geneLevel_n988.txt.gz", GEO_ACC))
stopifnot(file.exists(gene_file))

message(sprintf("[%s] Loading gene-level expression matrix: %s (%.1f MB)",
                SCRIPT_NAME, basename(gene_file),
                file.info(gene_file)$size / 1e6))

old_warn <- getOption("warn"); options(warn = 0)
expr_mat <- read.table(gene_file, sep = "\t", header = TRUE,
                        row.names = 1, check.names = FALSE, comment.char = "")
options(warn = old_warn)

# Ensure matrix type
expr_mat <- as.matrix(expr_mat)
message(sprintf("[%s] Expression matrix: %d genes × %d samples",
                SCRIPT_NAME, nrow(expr_mat), ncol(expr_mat)))
message(sprintf("[%s] Sample IDs from matrix (first 5): %s",
                SCRIPT_NAME, paste(head(colnames(expr_mat), 5), collapse = ", ")))

# --------------------------------------------------------------------------
# 6) Match patient IDs between pData and expression matrix columns
#    Gene matrix column names = numeric patient IDs (as character)
#    pData$patient_id_clean = same numeric IDs extracted from "patient id:ch1"
# --------------------------------------------------------------------------
mat_ids    <- colnames(expr_mat)
common_ids <- intersect(pdata_combined$patient_id_clean, mat_ids)

n_pdata  <- nrow(pdata_combined)
n_matrix <- length(mat_ids)
n_common <- length(common_ids)

# IDs present in pData but absent from expression matrix
ids_pdata_only  <- setdiff(pdata_combined$patient_id_clean, mat_ids)
# IDs present in matrix but absent from pData
ids_matrix_only <- setdiff(mat_ids, pdata_combined$patient_id_clean)

message(sprintf("[%s] pData rows: %d | Matrix cols: %d | Common: %d",
                SCRIPT_NAME, n_pdata, n_matrix, n_common))
if (length(ids_pdata_only) > 0)
  message(sprintf("[%s]   pData-only IDs (no expression): %s",
                  SCRIPT_NAME, paste(ids_pdata_only, collapse = ", ")))
if (length(ids_matrix_only) > 0)
  message(sprintf("[%s]   Matrix-only IDs (no clinical): %s",
                  SCRIPT_NAME, paste(ids_matrix_only, collapse = ", ")))

if (n_common < 900) {
  stop(sprintf("[%s] ID matching failure: only %d common IDs (expected ~988). Check ID format.",
               SCRIPT_NAME, n_common))
}

# Subset and align
expr_matched  <- expr_mat[, common_ids, drop = FALSE]
pdata_matched <- pdata_combined[match(common_ids, pdata_combined$patient_id_clean), ]
rownames(pdata_matched) <- common_ids

stopifnot(identical(colnames(expr_matched), pdata_matched$patient_id_clean))
message(sprintf("[%s] After alignment: %d samples matched", SCRIPT_NAME, ncol(expr_matched)))

# --------------------------------------------------------------------------
# 7) Intra-cohort Z-score (applied on top of batch-corrected values)
#    Z-score per gene across all 988 matched samples (intra-cohort rule)
# --------------------------------------------------------------------------
z_t <- t(expr_matched)    # samples × genes
old_warn <- getOption("warn"); options(warn = 0)
gene_means <- colMeans(z_t, na.rm = TRUE)
gene_sds   <- apply(z_t, 2, sd, na.rm = TRUE)
options(warn = old_warn)

zero_sd <- gene_sds == 0 | is.na(gene_sds)
if (any(zero_sd))
  message(sprintf("[%s] %d genes with zero/NA SD — set to 0 after z-score", SCRIPT_NAME, sum(zero_sd)))
gene_sds[zero_sd] <- 1
z_t[, zero_sd]    <- 0

z_mat <- sweep(sweep(z_t, 2, gene_means, "-"), 2, gene_sds, "/")
# z_mat: samples (rows) × genes (cols)

# --------------------------------------------------------------------------
# 8) CorePAM weights + score computation
# --------------------------------------------------------------------------
weights_df  <- strict_csv(file.path(PATHS$results$corepam, "CorePAM_weights.csv"))
weights_df  <- weights_df[weights_df$weight != 0, ]
panel_genes <- weights_df$gene
panel_w     <- setNames(weights_df$weight, weights_df$gene)

genes_present <- intersect(panel_genes, colnames(z_mat))
genes_missing <- setdiff(panel_genes, colnames(z_mat))
n_present     <- length(genes_present)
frac_present  <- n_present / length(panel_genes)

message(sprintf("[%s] CorePAM genes present: %d / %d (%.1f%%)",
                SCRIPT_NAME, n_present, length(panel_genes), frac_present * 100))
if (length(genes_missing) > 0)
  message(sprintf("[%s] Missing genes: %s", SCRIPT_NAME, paste(genes_missing, collapse = ", ")))

if (frac_present < FREEZE$min_genes_fraction)
  stop(sprintf("[%s] COVERAGE FAILURE: %.1f%% < %.0f%% minimum",
               SCRIPT_NAME, frac_present * 100, FREEZE$min_genes_fraction * 100))

w_present      <- panel_w[genes_present]
denom_sum_absw <- sum(abs(w_present))
score_raw      <- as.vector(z_mat[, genes_present] %*% w_present) / denom_sum_absw
names(score_raw) <- rownames(z_mat)

old_warn <- getOption("warn"); options(warn = 0)
score_z_raw <- as.vector(scale(score_raw))
options(warn = old_warn)
names(score_z_raw) <- rownames(z_mat)

# --------------------------------------------------------------------------
# 9) Assemble final data frame and save
# --------------------------------------------------------------------------
final_df <- tibble(
  sample_id      = common_ids,
  pcr            = pdata_matched$pcr,
  hr             = pdata_matched$hr,
  her2           = pdata_matched$her2,
  arm            = pdata_matched$arm,
  mp             = pdata_matched$mp,
  er_status      = dplyr::case_when(
                     pdata_matched$hr == 1L ~ "positive",
                     pdata_matched$hr == 0L ~ "negative",
                     TRUE                   ~ NA_character_
                   ),
  score          = score_raw[common_ids],
  score_z        = score_z_raw[common_ids],
  genes_present  = n_present,
  denom_sum_absw = denom_sum_absw,
  cohort         = COHORT
)

# Drop samples without pCR
n_before <- nrow(final_df)
final_df <- final_df[!is.na(final_df$pcr), ]
n_dropped_na_pcr <- n_before - nrow(final_df)
if (n_dropped_na_pcr > 0)
  message(sprintf("[%s] Dropped %d samples with NA pCR", SCRIPT_NAME, n_dropped_na_pcr))

# --------------------------------------------------------------------------
# Inclusion log — documents every exclusion step for transparent reporting
# --------------------------------------------------------------------------
inclusion_log <- data.frame(
  step         = c("pData from series matrix",
                   "Expression matrix columns",
                   "After ID matching (intersect)",
                   "pData-only IDs (missing expression)",
                   "Matrix-only IDs (missing clinical)",
                   "After dropping NA pCR",
                   "FINAL ANALYZED N"),
  n            = c(n_pdata,
                   n_matrix,
                   n_common,
                   length(ids_pdata_only),
                   length(ids_matrix_only),
                   n_dropped_na_pcr,
                   nrow(final_df)),
  note         = c("GPL20078 (654) + GPL30493 (334)",
                   "Gene-level batch-corrected file GSE194040_geneLevel_n988.txt.gz",
                   "Samples with both clinical and expression data",
                   paste0(if (length(ids_pdata_only) > 0) paste(ids_pdata_only, collapse=";") else "none",
                          " — excluded (no expression data found)"),
                   paste0(if (length(ids_matrix_only) > 0) paste(ids_matrix_only, collapse=";") else "none",
                          " — excluded (no clinical data found)"),
                   "NA in pcr:ch1 field",
                   "Used in all downstream analyses"),
  stringsAsFactors = FALSE
)

incl_log_path <- file.path(PATHS$results$pcr, "ISPY2_inclusion_log.csv")
write.csv(inclusion_log, incl_log_path, row.names = FALSE)
message(sprintf("[%s] Inclusion log saved: %s", SCRIPT_NAME, incl_log_path))

message(sprintf("[%s] Final N=%d | pCR+=%d (%.1f%%) | HR+=%d | HER2+=%d",
                SCRIPT_NAME,
                nrow(final_df),
                sum(final_df$pcr == 1L, na.rm = TRUE),
                100 * mean(final_df$pcr == 1L, na.rm = TRUE),
                sum(final_df$hr   == 1L, na.rm = TRUE),
                sum(final_df$her2 == 1L, na.rm = TRUE)))

out_dir <- proc_pcr_cohort(COHORT)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

old_warn <- getOption("warn"); options(warn = 0)
arrow::write_parquet(final_df, out_path)
options(warn = old_warn)

h  <- sha256_file(out_path)
sz <- file.info(out_path)$size / 1e6
registry_append(COHORT, "pcr_analysis_ready", out_path, h, "ok", SCRIPT_NAME, sz,
                list(n_samples      = nrow(final_df),
                     pcr1           = sum(final_df$pcr == 1L, na.rm = TRUE),
                     pcr_rate       = round(mean(final_df$pcr == 1L, na.rm = TRUE), 3),
                     genes_present  = n_present,
                     n_arms         = length(unique(final_df$arm[!is.na(final_df$arm)])),
                     platform       = PLATFORM))

message(sprintf("[%s] COMPLETED — %s | N=%d | pCR=%.1f%% | genes=%d/%d",
                SCRIPT_NAME, COHORT, nrow(final_df),
                100 * mean(final_df$pcr == 1L, na.rm = TRUE),
                n_present, length(panel_genes)))
message(sprintf("[%s] Output: %s (%.2f MB)", SCRIPT_NAME, out_path, sz))
