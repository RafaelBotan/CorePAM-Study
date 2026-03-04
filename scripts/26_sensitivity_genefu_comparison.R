# =============================================================================
# SCRIPT: 26_sensitivity_genefu_comparison.R
# PURPOSE: Sensitivity analysis — head-to-head C-index comparison of CorePAM
#          vs research-based PAM50-ROR-S and OncotypeDX-RS (genefu package).
#          These are heuristic estimates of commercial signatures, not exact
#          replications of proprietary algorithms.
# PROJECT: Core-PAM (Memorial CorePAM)
#
# OUTPUT:
#   results/supp/EXP_genefu_comparison.csv
# =============================================================================
source("scripts/00_setup.R")
SCRIPT_NAME <- "26_sensitivity_genefu_comparison.R"

suppressPackageStartupMessages({
  library(survival)
  library(arrow)
})

OUT_DIR <- file.path(PATHS$results$supp)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

message("[26_SENS] ========== GENEFU COMPARISON ==========")

# Install genefu if needed (only for centroid data)
if (!requireNamespace("genefu", quietly = TRUE)) {
  message("[26_SENS] Installing genefu from Bioconductor...")
  old_warn <- getOption("warn"); options(warn = 0)
  BiocManager::install("genefu", update = FALSE, ask = FALSE)
  options(warn = old_warn)
}

# Load genefu data (suppress namespace conflict warning)
old_warn <- getOption("warn"); options(warn = 0)
suppressPackageStartupMessages(library(genefu))
options(warn = old_warn)

# Load PAM50 centroids from genefu
data(pam50, package = "genefu")
pam50_centroids <- pam50$centroids       # 50 genes x 5 subtypes
pam50_map <- pam50$centroids.map         # gene annotation

# Load OncotypeDX signature from genefu
data(sig.oncotypedx, package = "genefu")
odx_sig <- sig.oncotypedx  # 21 genes with group + weight

COHORTS <- c("SCANB", "TCGA_BRCA", "METABRIC", "GSE20685")

results <- list()

for (coh in COHORTS) {
  message(sprintf("[26_SENS] Processing %s...", coh))

  ar <- as.data.frame(read_parquet(file.path(PATHS$processed, coh, "analysis_ready.parquet")))
  expr <- as.data.frame(read_parquet(file.path(PATHS$processed, coh, "expression_genelevel_preZ.parquet")))

  # Expression: genes as rows, samples as columns
  gene_col <- names(expr)[1]
  genes <- expr[[gene_col]]
  expr_mat <- as.matrix(expr[, -1])
  rownames(expr_mat) <- genes

  # ---- PAM50 subtyping (correlation to centroids) ----
  pam50_genes <- rownames(pam50_centroids)
  common_pam50 <- intersect(genes, pam50_genes)
  message(sprintf("[26_SENS] %s: PAM50 gene overlap = %d/%d", coh, length(common_pam50), length(pam50_genes)))

  ror_scores <- NULL
  pam50_subtypes_str <- "FAIL"

  if (length(common_pam50) >= 30) {
    # Subset to common genes
    expr_pam50 <- expr_mat[common_pam50, ]
    cent_pam50 <- pam50_centroids[common_pam50, ]

    # Compute Spearman correlation of each sample to each centroid
    n_samples <- ncol(expr_pam50)
    subtypes <- c("Basal", "Her2", "LumA", "LumB", "Normal")
    cor_mat <- matrix(NA, nrow = n_samples, ncol = 5)
    colnames(cor_mat) <- subtypes

    old_warn <- getOption("warn"); options(warn = 0)
    for (i in seq_len(n_samples)) {
      sample_expr <- expr_pam50[, i]
      for (j in seq_along(subtypes)) {
        cor_mat[i, j] <- cor(sample_expr, cent_pam50[, subtypes[j]], method = "spearman", use = "complete.obs")
      }
    }
    options(warn = old_warn)

    # Assign subtype by highest correlation
    assigned_subtypes <- subtypes[apply(cor_mat, 1, which.max)]
    names(assigned_subtypes) <- colnames(expr_mat)

    pam50_subtypes_str <- paste(names(table(assigned_subtypes)), collapse = ",")
    message(sprintf("[26_SENS] %s: PAM50 subtypes: %s", coh, paste(table(assigned_subtypes), collapse=", ")))

    # Compute ROR-S from centroid correlations
    # Convert correlations to "probabilities" via softmax
    exp_cor <- exp(cor_mat)
    probs <- exp_cor / rowSums(exp_cor)
    colnames(probs) <- subtypes

    # ROR-S = -0.34*LumA + 0.23*LumB + 0.05*Basal + 0.12*Her2 + 0*Normal
    # (Wallden et al. 2015 / Parker et al. 2009)
    ror_coeffs <- c("Basal" = 0.05, "Her2" = 0.12, "LumA" = -0.34, "LumB" = 0.23, "Normal" = 0)
    ror_raw <- rep(0, n_samples)
    for (st in names(ror_coeffs)) {
      ror_raw <- ror_raw + ror_coeffs[st] * probs[, st]
    }
    # Scale to 0-100
    ror_range <- range(ror_raw, na.rm = TRUE)
    if (diff(ror_range) > 0) {
      ror_scores <- (ror_raw - ror_range[1]) / diff(ror_range) * 100
      names(ror_scores) <- colnames(expr_mat)
    }
  }

  # ---- OncotypeDX Recurrence Score (research-based) ----
  # ODX uses 16 cancer-related genes + 5 reference genes, grouped into:
  # proliferation, invasion, her2, estrogen groups
  odx_genes <- odx_sig$symbol
  common_odx <- intersect(genes, odx_genes)
  message(sprintf("[26_SENS] %s: ODX gene overlap = %d/%d", coh, length(common_odx), length(odx_genes)))

  odx_scores <- NULL
  if (length(common_odx) >= 14) {  # Need at least ~2/3 of 21 genes
    # Compute group scores
    odx_sub <- odx_sig[odx_sig$symbol %in% common_odx, ]
    # Remove rows with NA weight
    odx_sub <- odx_sub[!is.na(odx_sub$weight), ]

    old_warn <- getOption("warn"); options(warn = 0)
    # Simple weighted sum approach (research-based approximation)
    n_samples <- ncol(expr_mat)
    odx_raw <- rep(0, n_samples)
    total_weight <- 0

    for (g in seq_len(nrow(odx_sub))) {
      gene_name <- odx_sub$symbol[g]
      w <- as.numeric(odx_sub$weight[g])
      if (!is.na(w) && gene_name %in% rownames(expr_mat)) {
        gene_vals <- as.numeric(expr_mat[gene_name, ])
        gene_vals[is.na(gene_vals)] <- 0  # impute NAs with 0
        odx_raw <- odx_raw + w * gene_vals
        total_weight <- total_weight + abs(w)
      }
    }
    options(warn = old_warn)

    if (!is.na(total_weight) && total_weight > 0) {
      odx_raw <- odx_raw / total_weight
      # Scale to 0-100
      odx_range <- range(odx_raw, na.rm = TRUE)
      if (diff(odx_range) > 0) {
        odx_scores <- (odx_raw - odx_range[1]) / diff(odx_range) * 100
        names(odx_scores) <- colnames(expr_mat)
      }
    }
  }

  # ---- Match samples and merge ----
  expr_sample_ids <- colnames(expr_mat)

  # Determine join key
  if ("sample_id" %in% names(ar)) {
    overlap_pid <- length(intersect(expr_sample_ids, ar$patient_id))
    overlap_sid <- length(intersect(expr_sample_ids, ar$sample_id))
    if (overlap_pid >= overlap_sid) {
      ar$join_key <- ar$patient_id
    } else {
      ar$join_key <- ar$sample_id
    }
  } else {
    ar$join_key <- ar$patient_id
  }

  # Build score data.frame
  score_df <- data.frame(
    join_key = expr_sample_ids,
    stringsAsFactors = FALSE
  )
  if (!is.null(ror_scores)) score_df$ror_s <- ror_scores[expr_sample_ids]
  if (!is.null(odx_scores)) score_df$odx_rs <- odx_scores[expr_sample_ids]

  merged <- merge(ar, score_df, by = "join_key")
  message(sprintf("[26_SENS] %s: merged N=%d, ROR=%s, ODX=%s",
                  coh, nrow(merged),
                  if(!is.null(ror_scores)) sprintf("%d valid", sum(!is.na(merged$ror_s))) else "FAIL",
                  if(!is.null(odx_scores)) sprintf("%d valid", sum(!is.na(merged$odx_rs))) else "FAIL"))

  # ---- Survival analysis: head-to-head C-index ----
  time_col  <- if (coh == "METABRIC") "dss_time_months" else "os_time_months"
  event_col <- if (coh == "METABRIC") "dss_event" else "os_event"

  merged <- merged[!is.na(merged[[time_col]]) & merged[[time_col]] > 0 &
                   !is.na(merged[[event_col]]), ]

  old_warn <- getOption("warn"); options(warn = 0)
  surv_obj <- Surv(merged[[time_col]], merged[[event_col]])

  # CorePAM C-index
  c_corepam <- tryCatch({
    fit <- coxph(surv_obj ~ score_z, data = merged)
    concordance(fit)$concordance
  }, error = function(e) NA_real_)

  # CorePAM HR
  hr_corepam <- NA_real_; p_corepam <- NA_real_
  tryCatch({
    fit_cp <- coxph(surv_obj ~ score_z, data = merged)
    hr_corepam <- exp(coef(fit_cp))
    p_corepam <- summary(fit_cp)$coefficients[,"Pr(>|z|)"]
  }, error = function(e) NULL)

  # ROR-S C-index
  c_ror <- NA_real_; hr_ror <- NA_real_; p_ror <- NA_real_
  if ("ror_s" %in% names(merged) && sum(!is.na(merged$ror_s)) > 50) {
    merged$ror_z <- scale(merged$ror_s)
    tryCatch({
      fit <- coxph(surv_obj ~ ror_z, data = merged)
      c_ror <- concordance(fit)$concordance
      hr_ror <- exp(coef(fit))
      p_ror <- summary(fit)$coefficients[,"Pr(>|z|)"]
    }, error = function(e) NULL)
  }

  # ODX RS C-index
  c_odx <- NA_real_; hr_odx <- NA_real_; p_odx <- NA_real_
  if ("odx_rs" %in% names(merged) && sum(!is.na(merged$odx_rs)) > 50) {
    merged$odx_z <- scale(merged$odx_rs)
    tryCatch({
      fit <- coxph(surv_obj ~ odx_z, data = merged)
      c_odx <- concordance(fit)$concordance
      hr_odx <- exp(coef(fit))
      p_odx <- summary(fit)$coefficients[,"Pr(>|z|)"]
    }, error = function(e) NULL)
  }

  options(warn = old_warn)

  # Number of genes used
  n_pam50_genes <- length(common_pam50)
  n_odx_genes <- length(common_odx)

  results[[coh]] <- data.frame(
    cohort = coh,
    n = nrow(merged),
    n_pam50_genes = n_pam50_genes,
    n_odx_genes = n_odx_genes,
    c_corepam = round(c_corepam, 4),
    hr_corepam = round(as.numeric(hr_corepam), 4),
    p_corepam = signif(as.numeric(p_corepam), 4),
    c_ror_s = round(c_ror, 4),
    hr_ror_s = round(as.numeric(hr_ror), 4),
    p_ror_s = signif(as.numeric(p_ror), 4),
    c_odx_rs = round(c_odx, 4),
    hr_odx_rs = round(as.numeric(hr_odx), 4),
    p_odx_rs = signif(as.numeric(p_odx), 4),
    pam50_subtypes = pam50_subtypes_str,
    stringsAsFactors = FALSE
  )

  message(sprintf("[26_SENS] %s: C_CorePAM=%.4f | C_ROR=%.4f | C_ODX=%.4f",
                  coh, c_corepam, c_ror, c_odx))
}

# --------------------------------------------------------------------------
# Save
# --------------------------------------------------------------------------
results_df <- do.call(rbind, results)
out_path <- file.path(OUT_DIR, "EXP_genefu_comparison.csv")
write.csv(results_df, out_path, row.names = FALSE)
message(sprintf("[26_SENS] Results saved to %s", out_path))
message("[26_SENS] ========== DONE ==========")
