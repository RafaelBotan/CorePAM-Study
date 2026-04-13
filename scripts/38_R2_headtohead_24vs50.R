# =============================================================================
# SCRIPT: 38_R2_headtohead_24vs50.R
# PURPOSE: Head-to-head external comparison of CorePAM (24 genes) vs
#          PAM50-full elastic-net (50 genes) — both framework-matched
#          (intra-cohort z-score) and deployment-matched (frozen z-score).
#          Addresses Reviewer 2 / GPT cross-check: demonstrate that
#          non-inferiority holds externally, not just in derivation.
# PROJECT: Core-PAM — Breast Cancer Research R2
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "38_R2_headtohead_24vs50.R"

suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
  library(dplyr)
})

set.seed(FREEZE$seed_folds)
message(sprintf("[%s] Starting head-to-head 24 vs 50 comparison", SCRIPT_NAME))

# --------------------------------------------------------------------------
# 1. Load frozen coefficients for both models
# --------------------------------------------------------------------------
w24_df <- read.csv(file.path(PATHS$results$corepam, "CorePAM_weights.csv"),
                   stringsAsFactors = FALSE)
w50_df <- read.csv(file.path(PATHS$results$corepam, "PAM50full_weights.csv"),
                   stringsAsFactors = FALSE)
w24 <- setNames(w24_df$weight, w24_df$gene)
w50 <- setNames(w50_df$weight, w50_df$gene)

message(sprintf("[%s] Loaded weights: 24-gene (%d), 50-gene (%d)",
                SCRIPT_NAME, length(w24), length(w50)))

# --------------------------------------------------------------------------
# 2. Load or build SCAN-B reference mean/SD for ALL 50 genes (frozen scoring)
# --------------------------------------------------------------------------
ref_out <- file.path(PATHS$results$corepam, "SCANB_reference_meanSD_50genes.csv")

if (file.exists(ref_out)) {
  scanb_ref <- read.csv(ref_out, stringsAsFactors = FALSE)
  message(sprintf("[%s] SCAN-B reference loaded from cache: %s (%d genes)",
                  SCRIPT_NAME, ref_out, nrow(scanb_ref)))
} else {
  scanb_expr_path <- file.path(proc_cohort("SCANB"),
                               "expression_genelevel_preZ.parquet")
  scanb_expr <- strict_parquet(scanb_expr_path)
  scanb_mat  <- as.matrix(scanb_expr[, -1])
  rownames(scanb_mat) <- scanb_expr$gene

  pam50_genes <- names(w50)
  scanb_ref <- data.frame(
    gene     = pam50_genes,
    mean_scanb = sapply(pam50_genes, function(g) {
      if (g %in% rownames(scanb_mat)) mean(scanb_mat[g, ], na.rm = TRUE) else NA
    }),
    sd_scanb   = sapply(pam50_genes, function(g) {
      if (g %in% rownames(scanb_mat)) sd(scanb_mat[g, ], na.rm = TRUE) else NA
    }),
    stringsAsFactors = FALSE
  )
  write.csv(scanb_ref, ref_out, row.names = FALSE)
  message(sprintf("[%s] SCAN-B reference computed and saved: %s (%d genes, %d valid)",
                  SCRIPT_NAME, ref_out, nrow(scanb_ref), sum(!is.na(scanb_ref$mean_scanb))))
  rm(scanb_expr, scanb_mat); gc(verbose = FALSE)
}

# --------------------------------------------------------------------------
# 3. Score computation function
# --------------------------------------------------------------------------
compute_lp <- function(expr_mat, weights, ref_meansd = NULL) {
  # expr_mat: genes x samples matrix
  # weights: named vector of gene weights
  # ref_meansd: NULL for intra-cohort z-score, or data.frame(gene, mean_scanb, sd_scanb)

  genes_avail <- intersect(names(weights), rownames(expr_mat))
  w <- weights[genes_avail]

  if (is.null(ref_meansd)) {
    # Intra-cohort z-score
    z <- t(scale(t(expr_mat[genes_avail, , drop = FALSE])))
  } else {
    # Frozen z-score using SCAN-B reference
    ref <- ref_meansd[match(genes_avail, ref_meansd$gene), ]
    z <- (expr_mat[genes_avail, , drop = FALSE] - ref$mean_scanb) / ref$sd_scanb
  }

  # Handle NAs in z
  z[is.na(z)] <- 0

  # Weighted sum / sum(|w|) — denominator renormalisation
  lp <- as.vector(t(z) %*% w / sum(abs(w)))
  names(lp) <- colnames(expr_mat)

  list(lp = lp, n_genes = length(genes_avail), genes_missing = setdiff(names(weights), genes_avail))
}

# --------------------------------------------------------------------------
# 4. Bootstrap paired delta-C
# --------------------------------------------------------------------------
boot_delta_c <- function(time, event, lp1, lp2, B = 2000) {
  n <- length(time)
  delta_boot <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample.int(n, replace = TRUE)
    c1 <- tryCatch(survival::concordance(Surv(time[idx], event[idx]) ~ I(-lp1[idx]))$concordance,
                   error = function(e) NA)
    c2 <- tryCatch(survival::concordance(Surv(time[idx], event[idx]) ~ I(-lp2[idx]))$concordance,
                   error = function(e) NA)
    delta_boot[b] <- c1 - c2
  }
  delta_boot <- delta_boot[!is.na(delta_boot)]
  c(mean = mean(delta_boot), lo = quantile(delta_boot, 0.025),
    hi = quantile(delta_boot, 0.975))
}

# --------------------------------------------------------------------------
# 5. Cohort definitions
# --------------------------------------------------------------------------
cohorts <- list(
  list(id = "TCGA_BRCA", time_col = "os_time_months", event_col = "os_event",
       endpoint = "OS"),
  list(id = "METABRIC", time_col = "dss_time_months", event_col = "dss_event",
       endpoint = "DSS"),
  list(id = "GSE20685", time_col = "os_time_months", event_col = "os_event",
       endpoint = "OS"),
  list(id = "GSE1456", time_col = "os_time_months", event_col = "os_event",
       endpoint = "OS")
)

# --------------------------------------------------------------------------
# 6. Main loop
# --------------------------------------------------------------------------
results <- list()

for (coh in cohorts) {
  message(sprintf("\n[%s] === %s (%s) ===", SCRIPT_NAME, coh$id, coh$endpoint))

  # Load expression
  expr_path <- file.path(proc_cohort(coh$id), "expression_genelevel_preZ.parquet")
  expr_df   <- strict_parquet(expr_path)
  expr_mat  <- as.matrix(expr_df[, -1])
  rownames(expr_mat) <- expr_df$gene

  # Load clinical/analysis-ready
  ar_path <- file.path(proc_cohort(coh$id), "analysis_ready.parquet")
  ar      <- strict_parquet(ar_path)

  # Filter valid
  ar <- ar[!is.na(ar[[coh$time_col]]) & ar[[coh$time_col]] > 0 &
             !is.na(ar[[coh$event_col]]), ]

  # Detect ID column (TCGA uses patient_id, others use sample_id)
  id_col <- if ("sample_id" %in% names(ar)) "sample_id" else "patient_id"

  # Match samples
  common_samples <- intersect(ar[[id_col]], colnames(expr_mat))
  ar       <- ar[match(common_samples, ar[[id_col]]), ]
  expr_sub <- expr_mat[, common_samples, drop = FALSE]

  time_vec  <- ar[[coh$time_col]]
  event_vec <- ar[[coh$event_col]]
  N <- length(time_vec)
  n_events <- sum(event_vec)

  message(sprintf("[%s] %s: N=%d, Events=%d", SCRIPT_NAME, coh$id, N, n_events))

  # --- A. Framework-matched (intra-cohort z-score) ---
  lp24_intra <- compute_lp(expr_sub, w24, ref_meansd = NULL)
  lp50_intra <- compute_lp(expr_sub, w50, ref_meansd = NULL)

  c24_intra <- concordance(Surv(time_vec, event_vec) ~ I(-lp24_intra$lp))$concordance
  c50_intra <- concordance(Surv(time_vec, event_vec) ~ I(-lp50_intra$lp))$concordance
  dc_intra  <- c24_intra - c50_intra

  # Bootstrap paired ΔC
  boot_intra <- boot_delta_c(time_vec, event_vec, lp24_intra$lp, lp50_intra$lp, B = 2000)

  # Correlation
  cor_intra <- cor(lp24_intra$lp, lp50_intra$lp, method = "pearson")

  message(sprintf("[%s] %s INTRA: C24=%.4f C50=%.4f ΔC=%.4f (%.4f,%.4f) r=%.3f | genes24=%d genes50=%d",
                  SCRIPT_NAME, coh$id, c24_intra, c50_intra, dc_intra,
                  boot_intra["lo.2.5%"], boot_intra["hi.97.5%"],
                  cor_intra, lp24_intra$n_genes, lp50_intra$n_genes))

  # --- B. Deployment-matched (frozen z-score) ---
  lp24_frozen <- compute_lp(expr_sub, w24, ref_meansd = scanb_ref)
  lp50_frozen <- compute_lp(expr_sub, w50, ref_meansd = scanb_ref)

  c24_frozen <- concordance(Surv(time_vec, event_vec) ~ I(-lp24_frozen$lp))$concordance
  c50_frozen <- concordance(Surv(time_vec, event_vec) ~ I(-lp50_frozen$lp))$concordance
  dc_frozen  <- c24_frozen - c50_frozen

  boot_frozen <- boot_delta_c(time_vec, event_vec, lp24_frozen$lp, lp50_frozen$lp, B = 2000)
  cor_frozen  <- cor(lp24_frozen$lp, lp50_frozen$lp, method = "pearson")

  message(sprintf("[%s] %s FROZEN: C24=%.4f C50=%.4f ΔC=%.4f (%.4f,%.4f) r=%.3f",
                  SCRIPT_NAME, coh$id, c24_frozen, c50_frozen, dc_frozen,
                  boot_frozen["lo.2.5%"], boot_frozen["hi.97.5%"], cor_frozen))

  # Collect results
  results[[length(results) + 1]] <- data.frame(
    cohort = coh$id, endpoint = coh$endpoint, N = N, events = n_events,
    mode = "intra",
    c24 = c24_intra, c50 = c50_intra, delta_c = dc_intra,
    dc_lo = boot_intra["lo.2.5%"], dc_hi = boot_intra["hi.97.5%"],
    cor_lp = cor_intra,
    genes24 = lp24_intra$n_genes, genes50 = lp50_intra$n_genes,
    stringsAsFactors = FALSE
  )
  results[[length(results) + 1]] <- data.frame(
    cohort = coh$id, endpoint = coh$endpoint, N = N, events = n_events,
    mode = "frozen",
    c24 = c24_frozen, c50 = c50_frozen, delta_c = dc_frozen,
    dc_lo = boot_frozen["lo.2.5%"], dc_hi = boot_frozen["hi.97.5%"],
    cor_lp = cor_frozen,
    genes24 = lp24_frozen$n_genes, genes50 = lp50_frozen$n_genes,
    stringsAsFactors = FALSE
  )

  rm(expr_df, expr_mat, expr_sub, ar); gc(verbose = FALSE)
}

# --------------------------------------------------------------------------
# 7. Combine and save results
# --------------------------------------------------------------------------
res_df <- bind_rows(results)
rownames(res_df) <- NULL

out_csv <- file.path(PATHS$results$supp, "headtohead_24vs50_external.csv")
write.csv(res_df, out_csv, row.names = FALSE)
message(sprintf("\n[%s] Results saved: %s", SCRIPT_NAME, out_csv))

# Print summary
message("\n========== HEAD-TO-HEAD SUMMARY ==========")
for (i in seq_len(nrow(res_df))) {
  r <- res_df[i, ]
  ni_flag <- ifelse(r$dc_lo > -0.010, "NI-OK", "")
  message(sprintf("  %s [%s] C24=%.4f C50=%.4f ΔC=%+.4f (%.4f,%.4f) r=%.3f %s",
                  r$cohort, r$mode, r$c24, r$c50, r$delta_c,
                  r$dc_lo, r$dc_hi, r$cor_lp, ni_flag))
}

# --------------------------------------------------------------------------
# 8. Forest plot of ΔC = C24 - C50
# --------------------------------------------------------------------------
plot_df <- res_df %>%
  mutate(
    label = paste0(cohort, " (", mode, ")"),
    label = factor(label, levels = rev(unique(label)))
  )

p <- ggplot(plot_df, aes(x = delta_c, y = label)) +
  geom_vline(xintercept = 0, linetype = "solid", colour = "grey50") +
  geom_vline(xintercept = -0.010, linetype = "dashed", colour = "red", linewidth = 0.5) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = dc_lo, xmax = dc_hi), width = 0.2, orientation = "y") +
  labs(
    x = expression(Delta*C ~ "(C"[24] ~ "-" ~ "C"[50]*")"),
    y = NULL,
    title = "External head-to-head: CorePAM (24) vs PAM50-full (50)",
    subtitle = "Dashed red line: non-inferiority margin (−0.010)"
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12))

out_pdf <- file.path(PATHS$figures$supp_en_pdf, "FigS_HeadToHead_24vs50.pdf")
out_png <- file.path(PATHS$figures$supp_en_png, "FigS_HeadToHead_24vs50.png")

cairo_pdf(out_pdf, width = 8, height = 5)
print(p)
dev.off()

png(out_png, width = 8, height = 5, units = "in", res = 300)
print(p)
dev.off()

message(sprintf("[%s] Figure saved: %s", SCRIPT_NAME, out_png))
message(sprintf("[%s] COMPLETED", SCRIPT_NAME))
