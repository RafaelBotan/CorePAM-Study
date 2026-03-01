# =============================================================================
# SCRIPT: 20_prepare_pCR_ISPY1.R
# PURPOSE: Harmonize clinical (pCR), preprocess expression (Affymetrix
#          HGU133A), compute CorePAM score for I-SPY1 (GSE22226).
#          Produces analysis_ready.parquet for 21_pCR_logistic_analysis.R.
#
# COHORT:  ISPY1 | GEO: GSE22226 | Endpoint: pCR
# PLATFORM: HGU133A (GPL96) — NOTE: same restricted gene space as GSE20685.
#           CorePAM coverage may be < 100% due to platform limitations.
#           Minimum gene coverage: 80% (FREEZE$min_genes_fraction).
# N expected: ~149 pre-treatment biopsies
# REF: Esserman LJ et al. Sci Transl Med 2012; NACT (AC → T)
#
# PROJETO: Core-PAM (Memorial v6.1 / Freeze Core-PAM — pCR block)
# =============================================================================

source("scripts/00_setup.R")
source("scripts/03_utils_gene_mapping.R")

SCRIPT_NAME <- "20_prepare_pCR_ISPY1.R"
COHORT      <- "ISPY1"
GEO_ACC     <- "GSE22226"
PLATFORM    <- "hgu133a"    # GPL96 — limited to ~12970 genes

FORCE <- as.logical(Sys.getenv("FORCE_RERUN", "FALSE"))

out_path <- file.path(proc_pcr_cohort(COHORT), "analysis_ready.parquet")
if (!FORCE && file.exists(out_path)) {
  message(sprintf("[%s] Output exists — skipping. Set FORCE_RERUN=TRUE to rerun.",
                  SCRIPT_NAME))
  quit(save = "no", status = 0)
}

message(sprintf("[%s] Starting pCR prepare for %s (%s) [HGU133A]",
                SCRIPT_NAME, COHORT, GEO_ACC))

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
  message(sprintf("[%s] Downloading %s...", SCRIPT_NAME, GEO_ACC))
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

# GSE22226 may have multiple ExpressionSets (multiple platforms or sub-series)
# Pick the one with the most samples on HGU133A
if (is.list(gse)) {
  sizes     <- sapply(gse, function(x) ncol(Biobase::exprs(x)))
  platforms <- sapply(gse, function(x) annotation(x))
  message(sprintf("[%s] Multiple ExpressionSets: %s",
                  SCRIPT_NAME, paste(names(gse), platforms, sep = "=", collapse = ", ")))
  # Prefer hgu133a or the largest
  hgu_idx <- which(grepl("hgu133a", platforms, ignore.case = TRUE))
  idx     <- if (length(hgu_idx) > 0) hgu_idx[which.max(sizes[hgu_idx])] else which.max(sizes)
  gse     <- gse[[idx]]
  message(sprintf("[%s] Selected ExpressionSet index %d (%s, n=%d)",
                  SCRIPT_NAME, idx, annotation(gse), ncol(Biobase::exprs(gse))))
}

expr_mat  <- Biobase::exprs(gse)
pdata     <- Biobase::pData(gse)
probe_ids <- rownames(expr_mat)

message(sprintf("[%s] GEO loaded: %d probes × %d samples", SCRIPT_NAME,
                nrow(expr_mat), ncol(expr_mat)))
message(sprintf("[%s] pData columns: %s", SCRIPT_NAME,
                paste(names(pdata), collapse = ", ")))

# --------------------------------------------------------------------------
# 2) Extract pCR and covariates
# --------------------------------------------------------------------------
source("scripts/20_utils_pcr_extract.R")
pcr_result <- extract_pcr_column(pdata, cohort = COHORT)
pdata$pcr  <- pcr_result$pcr
message(sprintf("[%s] pCR column: '%s' | pCR+=%d | pCR-=%d | NA=%d",
                SCRIPT_NAME, pcr_result$col_used,
                sum(pdata$pcr == 1, na.rm = TRUE),
                sum(pdata$pcr == 0, na.rm = TRUE),
                sum(is.na(pdata$pcr))))

pdata$age       <- extract_numeric_covariate(pdata, pattern = "age")
pdata$er_status <- extract_binary_covariate(pdata, pattern = "er.status|estrogen")

# --------------------------------------------------------------------------
# 3) Probe → gene mapping (HGU133A — GPL96)
# --------------------------------------------------------------------------
probe_map <- build_hgnc_map_affymetrix(probe_ids, platform = PLATFORM)
probe_map  <- probe_map[!is.na(probe_map$SYMBOL) & nchar(probe_map$SYMBOL) > 0, ]

mapped_probes  <- intersect(rownames(expr_mat), probe_map$PROBEID)
expr_sub       <- expr_mat[mapped_probes, , drop = FALSE]
gene_labels    <- probe_map$SYMBOL[match(mapped_probes, probe_map$PROBEID)]
expr_collapsed <- collapse_probes_by_variance(expr_sub, gene_labels)

message(sprintf("[%s] Genes after collapse: %d (HGU133A has ~12970 genes)",
                SCRIPT_NAME, nrow(expr_collapsed)))
save_gene_mapping_audit(
  data.frame(cohort = COHORT, n_probes_mapped = length(mapped_probes),
             n_genes_after_collapse = nrow(expr_collapsed)),
  cohort = COHORT, script_name = SCRIPT_NAME
)

# --------------------------------------------------------------------------
# 4) CorePAM weights + Z-score + score (coverage check critical for HGU133A)
# --------------------------------------------------------------------------
weights_df     <- strict_csv(file.path(PATHS$results$corepam, "CorePAM_weights.csv"))
weights_df     <- weights_df[weights_df$weight != 0, ]
panel_genes    <- weights_df$gene
panel_weights  <- setNames(weights_df$weight, weights_df$gene)

z_t <- t(expr_collapsed)
old_warn <- getOption("warn"); options(warn = 0)
gene_means <- colMeans(z_t, na.rm = TRUE)
gene_sds   <- apply(z_t, 2, sd, na.rm = TRUE)
options(warn = old_warn)
zero_sd <- gene_sds == 0 | is.na(gene_sds)
gene_sds[zero_sd] <- 1; z_t[, zero_sd] <- 0
z_mat <- sweep(sweep(z_t, 2, gene_means, "-"), 2, gene_sds, "/")

genes_present <- intersect(panel_genes, colnames(z_mat))
genes_missing <- setdiff(panel_genes, colnames(z_mat))
n_present     <- length(genes_present)
frac_present  <- n_present / length(panel_genes)

message(sprintf("[%s] CorePAM genes present: %d / %d (%.1f%%)",
                SCRIPT_NAME, n_present, length(panel_genes), frac_present * 100))
if (length(genes_missing) > 0)
  message(sprintf("[%s] Missing (HGU133A platform ceiling): %s",
                  SCRIPT_NAME, paste(genes_missing, collapse = ", ")))

if (frac_present < FREEZE$min_genes_fraction)
  stop(sprintf("[%s] COVERAGE FAILURE: %.1f%% < %.0f%% minimum — HGU133A insufficient",
               SCRIPT_NAME, frac_present * 100, FREEZE$min_genes_fraction * 100))

w_present      <- panel_weights[genes_present]
denom_sum_absw <- sum(abs(w_present))
score_raw      <- as.vector(z_mat[, genes_present] %*% w_present) / denom_sum_absw
old_warn <- getOption("warn"); options(warn = 0)
score_z_raw <- as.vector(scale(score_raw))
options(warn = old_warn)

# --------------------------------------------------------------------------
# 5) Assemble and save
# --------------------------------------------------------------------------
geo_ids  <- colnames(expr_collapsed)
score_df <- data.frame(geo_id = geo_ids, score = score_raw, score_z = score_z_raw,
                       genes_present = n_present, denom_sum_absw = denom_sum_absw,
                       stringsAsFactors = FALSE)
pdata$geo_id <- rownames(pdata)
df_joined <- merge(pdata, score_df, by = "geo_id", sort = FALSE)

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
final_df <- final_df[!is.na(final_df$pcr), ]

out_dir <- proc_pcr_cohort(COHORT)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
old_warn <- getOption("warn"); options(warn = 0)
arrow::write_parquet(final_df, out_path)
options(warn = old_warn)

h  <- sha256_file(out_path); sz <- file.info(out_path)$size / 1e6
registry_append(COHORT, "pcr_analysis_ready", out_path, h, "ok", SCRIPT_NAME, sz,
                list(n_samples = nrow(final_df), pcr1 = sum(final_df$pcr),
                     genes_present = n_present, platform = PLATFORM))

message(sprintf("[%s] COMPLETED for %s | N=%d | pCR rate=%.1f%% | genes=%.0f%%",
                SCRIPT_NAME, COHORT, nrow(final_df),
                100 * mean(final_df$pcr, na.rm = TRUE),
                100 * frac_present))
