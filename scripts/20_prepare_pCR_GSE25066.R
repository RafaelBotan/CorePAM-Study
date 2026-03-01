# =============================================================================
# SCRIPT: 20_prepare_pCR_GSE25066.R
# PURPOSE: Harmonize clinical (pCR), preprocess expression (Affymetrix
#          HGU133Plus2), compute CorePAM score for GSE25066 (Hatzis 2011).
#          Produces analysis_ready.parquet for 21_pCR_logistic_analysis.R.
#
# COHORT:  GSE25066 | Endpoint: pCR (pathologic complete response to NACT)
# PLATFORM: HGU133Plus2 (GPL570)
# N expected: ~508 pre-treatment biopsies
# REF: Hatzis C et al. JAMA 2011;305:1873-1881
#
# PROJETO: Core-PAM (Memorial v6.1 / Freeze Core-PAM — pCR block)
# =============================================================================

source("scripts/00_setup.R")
source("scripts/03_utils_gene_mapping.R")

SCRIPT_NAME <- "20_prepare_pCR_GSE25066.R"
COHORT      <- "GSE25066"
GEO_ACC     <- "GSE25066"
PLATFORM    <- "hgu133plus2"

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
# 1) Load GEO series matrix
# --------------------------------------------------------------------------
raw_dir <- raw_pcr_cohort(COHORT)
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)

sm_path <- file.path(raw_dir, sprintf("%s_series_matrix.txt.gz", GEO_ACC))
if (!file.exists(sm_path)) {
  message(sprintf("[%s] Series matrix not found; downloading %s...", SCRIPT_NAME, GEO_ACC))
  old_warn2 <- getOption("warn"); options(warn = 0)
  GEOquery::getGEO(GEO = GEO_ACC, destdir = raw_dir, GSEMatrix = TRUE,
                   getGPL = FALSE, AnnotGPL = FALSE)
  options(warn = old_warn2)
  sm_files <- list.files(raw_dir, pattern = "series_matrix.txt.gz", full.names = TRUE)
  if (length(sm_files) > 0 && sm_files[1] != sm_path) file.rename(sm_files[1], sm_path)
}

stopifnot(file.exists(sm_path))
message(sprintf("[%s] Reading series matrix: %s", SCRIPT_NAME, basename(sm_path)))

old_warn <- getOption("warn"); options(warn = 0)
gse <- GEOquery::getGEO(filename = sm_path)
options(warn = old_warn)

expr_mat  <- Biobase::exprs(gse)        # probes × samples (log2 intensity)
pdata     <- Biobase::pData(gse)
probe_ids <- rownames(expr_mat)
sample_ids <- colnames(expr_mat)

message(sprintf("[%s] GEO loaded: %d probes × %d samples", SCRIPT_NAME,
                nrow(expr_mat), ncol(expr_mat)))
message(sprintf("[%s] pData columns: %s", SCRIPT_NAME,
                paste(names(pdata), collapse = ", ")))

# --------------------------------------------------------------------------
# 2) Extract pCR from pData (auto-detect column)
# --------------------------------------------------------------------------
# pCR values in GSE25066 are encoded in "characteristics_ch1.X" columns.
# The column containing "pCR" or "pathological complete response" is searched.
source("scripts/20_utils_pcr_extract.R")
pcr_result <- extract_pcr_column(pdata, cohort = COHORT)
pdata$pcr  <- pcr_result$pcr
message(sprintf("[%s] pCR column used: %s | pCR=1: %d | pCR=0: %d | NA: %d",
                SCRIPT_NAME, pcr_result$col_used,
                sum(pdata$pcr == 1, na.rm = TRUE),
                sum(pdata$pcr == 0, na.rm = TRUE),
                sum(is.na(pdata$pcr))))

# --------------------------------------------------------------------------
# 3) Extract optional covariates (age, er_status if available)
# --------------------------------------------------------------------------
pdata$age       <- extract_numeric_covariate(pdata, pattern = "age")
pdata$er_status <- extract_binary_covariate(pdata,
                     pattern = "er.status|estrogen.receptor",
                     pos_vals = c("positive", "pos", "er+", "1"),
                     neg_vals = c("negative", "neg", "er-", "0"))

# --------------------------------------------------------------------------
# 4) Probe → gene mapping (HGU133Plus2 → HGNC symbol)
# --------------------------------------------------------------------------
probe_map <- build_hgnc_map_affymetrix(probe_ids, platform = PLATFORM)
probe_map  <- probe_map[!is.na(probe_map$SYMBOL) & nchar(probe_map$SYMBOL) > 0, ]
message(sprintf("[%s] Probe map: %d probe→symbol pairs", SCRIPT_NAME, nrow(probe_map)))

# Align expression matrix to mapped probes
mapped_probes <- intersect(rownames(expr_mat), probe_map$PROBEID)
expr_sub      <- expr_mat[mapped_probes, , drop = FALSE]
gene_labels   <- probe_map$SYMBOL[match(mapped_probes, probe_map$PROBEID)]

# Collapse: highest-variance probe per gene
expr_collapsed <- collapse_probes_by_variance(expr_sub, gene_labels)
message(sprintf("[%s] Expression after collapse: %d genes × %d samples",
                SCRIPT_NAME, nrow(expr_collapsed), ncol(expr_collapsed)))

# Gene mapping audit
audit_df <- data.frame(
  probe_id  = mapped_probes,
  gene      = gene_labels,
  stringsAsFactors = FALSE
)
save_gene_mapping_audit(
  data.frame(cohort = COHORT, n_probes_mapped = nrow(audit_df),
             n_genes_after_collapse = nrow(expr_collapsed)),
  cohort = COHORT, script_name = SCRIPT_NAME
)

# --------------------------------------------------------------------------
# 5) Load CorePAM frozen weights
# --------------------------------------------------------------------------
weights_path   <- file.path(PATHS$results$corepam, "CorePAM_weights.csv")
weights_df     <- strict_csv(weights_path)
weights_df     <- weights_df[weights_df$weight != 0, ]
panel_genes    <- weights_df$gene
panel_weights  <- setNames(weights_df$weight, weights_df$gene)
message(sprintf("[%s] CorePAM panel: %d genes", SCRIPT_NAME, length(panel_genes)))

# --------------------------------------------------------------------------
# 6) Intra-cohort Z-score per gene
# --------------------------------------------------------------------------
# expr_collapsed is gene × sample; transpose for Z-score
z_t <- t(expr_collapsed)   # samples × genes

old_warn <- getOption("warn"); options(warn = 0)
gene_means <- colMeans(z_t, na.rm = TRUE)
gene_sds   <- apply(z_t, 2, sd, na.rm = TRUE)
options(warn = old_warn)

zero_sd <- gene_sds == 0 | is.na(gene_sds)
gene_sds[zero_sd] <- 1
z_t[, zero_sd] <- 0
z_mat <- sweep(sweep(z_t, 2, gene_means, "-"), 2, gene_sds, "/")

# --------------------------------------------------------------------------
# 7) CorePAM score = Σ(wᵢ·zᵢ) / Σ|wᵢ|  (present genes only)
# --------------------------------------------------------------------------
genes_present <- intersect(panel_genes, colnames(z_mat))
genes_missing <- setdiff(panel_genes, colnames(z_mat))
n_present     <- length(genes_present)
n_panel       <- length(panel_genes)
frac_present  <- n_present / n_panel

message(sprintf("[%s] Genes present: %d / %d (%.1f%%)",
                SCRIPT_NAME, n_present, n_panel, frac_present * 100))
if (length(genes_missing) > 0)
  message(sprintf("[%s] Missing: %s", SCRIPT_NAME, paste(genes_missing, collapse = ", ")))

if (frac_present < FREEZE$min_genes_fraction) {
  stop(sprintf("[%s] COVERAGE FAILURE: %.1f%% < %.0f%% minimum",
               SCRIPT_NAME, frac_present * 100, FREEZE$min_genes_fraction * 100))
}

w_present      <- panel_weights[genes_present]
denom_sum_absw <- sum(abs(w_present))
score_raw      <- as.vector(z_mat[, genes_present, drop = FALSE] %*% w_present) /
                  denom_sum_absw

old_warn <- getOption("warn"); options(warn = 0)
score_z_raw <- as.vector(scale(score_raw))
options(warn = old_warn)

# --------------------------------------------------------------------------
# 8) Assemble final data frame and save analysis_ready.parquet
# --------------------------------------------------------------------------
# Match samples: GEO sample IDs vs pData rownames
geo_ids  <- colnames(expr_collapsed)   # same as sample_ids after collapse
pdata_id <- rownames(pdata)

# pData rownames = GSM IDs; z_mat rownames = same order as expr_mat colnames
# Build score data.frame indexed by geo_ids
score_df <- data.frame(
  geo_id         = geo_ids,
  score          = score_raw,
  score_z        = score_z_raw,
  genes_present  = n_present,
  denom_sum_absw = denom_sum_absw,
  stringsAsFactors = FALSE
)

pdata$geo_id <- rownames(pdata)
df_joined <- merge(pdata, score_df, by = "geo_id", sort = FALSE)
message(sprintf("[%s] Final join: %d samples (pCR=1: %d | pCR=0: %d | NA: %d)",
                SCRIPT_NAME, nrow(df_joined),
                sum(df_joined$pcr == 1, na.rm = TRUE),
                sum(df_joined$pcr == 0, na.rm = TRUE),
                sum(is.na(df_joined$pcr))))

# Build final tibble with mandatory columns
final_df <- tibble(
  sample_id      = df_joined$geo_id,
  pcr            = df_joined$pcr,
  age            = if ("age" %in% names(df_joined)) df_joined$age else NA_real_,
  er_status      = if ("er_status" %in% names(df_joined)) df_joined$er_status else NA_character_,
  score          = df_joined$score,
  score_z        = df_joined$score_z,
  genes_present  = df_joined$genes_present,
  denom_sum_absw = df_joined$denom_sum_absw,
  cohort         = COHORT
)

# Remove samples with no pCR information (cannot be analysed)
n_before_drop <- nrow(final_df)
final_df <- final_df[!is.na(final_df$pcr), ]
message(sprintf("[%s] After removing NA pCR: %d samples (dropped %d)",
                SCRIPT_NAME, nrow(final_df), n_before_drop - nrow(final_df)))

out_dir <- proc_pcr_cohort(COHORT)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

old_warn <- getOption("warn"); options(warn = 0)
arrow::write_parquet(final_df, out_path)
options(warn = old_warn)

h  <- sha256_file(out_path)
sz <- file.info(out_path)$size / 1e6
message(sprintf("[%s] Saved: %s (%.2f MB | SHA256: %s)", SCRIPT_NAME, out_path, sz, h))
registry_append(COHORT, "pcr_analysis_ready", out_path, h, "ok", SCRIPT_NAME, sz,
                list(n_samples = nrow(final_df), pcr1 = sum(final_df$pcr),
                     genes_present = n_present))

message(sprintf("[%s] COMPLETED for %s | N=%d | pCR rate=%.1f%%",
                SCRIPT_NAME, COHORT, nrow(final_df),
                100 * mean(final_df$pcr, na.rm = TRUE)))
