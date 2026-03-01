# =============================================================================
# SCRIPT: 07y_pam50full_comparison.R
# PURPOSE: Extract PAM50-full weights (max df lambda) from CorePAM_model.rds,
#          compute PAM50full scores in all 4 OS cohorts, compare with CorePAM
#          scores, generate scatter + forest figures.
# OUTPUT:
#   results/corepam/PAM50full_weights.csv
#   results/supp/pam50full_comparison.csv
#   figures/supp/FigS_CorePAM_vs_PAM50full_EN.pdf/.png  (scatter, overwrite)
#   figures/supp/FigS_CorePAM_vs_PAM50full_HR_EN.pdf/.png (forest, new)
# PROJETO: Core-PAM (Memorial v6.1 / Freeze Core-PAM)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(survival)
  library(ggplot2)
  library(patchwork)
})

setwd("Y:/Phd-Genomic-claude")
source("scripts/00_setup.R")
SCRIPT_NAME <- "07y_pam50full_comparison.R"

message(sprintf("[%s] Starting PAM50full extraction + comparison", SCRIPT_NAME))

# --------------------------------------------------------------------------
# 0) Paths
# --------------------------------------------------------------------------
out_weights  <- file.path(PATHS$results$corepam, "PAM50full_weights.csv")
out_comp_csv <- file.path(PATHS$results$supp,    "pam50full_comparison.csv")

# --------------------------------------------------------------------------
# Helper: save PDF + PNG — EN only (supp/en/)
# --------------------------------------------------------------------------
save_fig <- function(p, name, w = 10, h = 8, lang = "en", dpi = 300) {
  pdf_path <- file.path(PATHS$figures[[paste0("supp_", lang, "_pdf")]], paste0(name, ".pdf"))
  png_path <- file.path(PATHS$figures[[paste0("supp_", lang, "_png")]], paste0(name, ".png"))
  old_warn <- getOption("warn"); options(warn = 0)
  tryCatch({
    cairo_pdf(pdf_path, width = w, height = h); print(p); dev.off()
    png(png_path, width = w, height = h, units = "in", res = dpi); print(p); dev.off()
    registry_append("ALL", name, pdf_path, sha256_file(pdf_path), "ok",
                    SCRIPT_NAME, file.info(pdf_path)$size / 1e6)
    message(sprintf("[%s] Saved: %s | %s", SCRIPT_NAME, pdf_path, png_path))
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    message(sprintf("[%s] ERROR saving %s: %s", SCRIPT_NAME, name, e$message))
  })
  options(warn = old_warn)
  invisible(list(pdf = pdf_path, png = png_path))
}

# --------------------------------------------------------------------------
# Helper: compute score from z-matrix and named weight vector
# --------------------------------------------------------------------------
compute_score <- function(z_mat, named_weights) {
  genes_present <- intersect(names(named_weights), colnames(z_mat))
  if (length(genes_present) == 0) stop("No matching genes between weights and z-matrix")
  w_present      <- named_weights[genes_present]
  denom_sum_absw <- sum(abs(w_present))
  z_sub          <- z_mat[, genes_present, drop = FALSE]
  score_raw      <- as.vector(z_sub %*% w_present) / denom_sum_absw
  list(
    score         = score_raw,
    genes_present = genes_present,
    n_present     = length(genes_present),
    frac_present  = length(genes_present) / length(named_weights)
  )
}

# --------------------------------------------------------------------------
# Helper: bootstrap C-index (Harrell, adjusted so C >= 0.5)
# --------------------------------------------------------------------------
bootstrap_cindex <- function(time, event, score_z, n_boot = 500, seed = 42) {
  set.seed(seed)
  n     <- length(time)
  cvals <- numeric(n_boot)
  old_warn <- getOption("warn"); options(warn = 0)
  for (i in seq_len(n_boot)) {
    idx    <- sample(n, n, replace = TRUE)
    cx_raw <- tryCatch(
      concordance(Surv(time[idx], event[idx]) ~ score_z[idx])$concordance,
      error = function(e) NA_real_
    )
    cvals[i] <- if (is.na(cx_raw)) NA_real_ else max(cx_raw, 1 - cx_raw)
  }
  options(warn = old_warn)
  old_warn2 <- getOption("warn"); options(warn = 0)
  c_raw <- concordance(Surv(time, event) ~ score_z)$concordance
  options(warn = old_warn2)
  list(
    c_index  = max(c_raw, 1 - c_raw),
    ci_low   = quantile(cvals, 0.025, na.rm = TRUE),
    ci_high  = quantile(cvals, 0.975, na.rm = TRUE)
  )
}

# =============================================================================
# STEP 1: Extract PAM50-full weights from glmnet model (max df lambda)
# =============================================================================
message(sprintf("[%s] === Step 1: Extract PAM50full weights ===", SCRIPT_NAME))

model_path <- file.path(PATHS$results$corepam, "CorePAM_model.rds")
model_obj  <- strict_rds(model_path)

# The model is a list; glmnet_fit contains the coxnet fit
glm_fit <- model_obj$glmnet_fit

# Find the lambda index with maximum non-zero genes (max df)
df_vals    <- glm_fit$df          # number of non-zero coefficients per lambda
lambda_idx <- which.max(df_vals)  # index with max df

message(sprintf("[%s] Max df = %d at lambda index %d (lambda = %.6f)",
                SCRIPT_NAME, max(df_vals), lambda_idx, glm_fit$lambda[lambda_idx]))

# Access beta matrix directly (sparse matrix: rows=genes, cols=lambdas)
# glmnet$beta has rownames = gene names; use beta[, lambda_idx] for that lambda
beta_mat <- glm_fit$beta
b_vec    <- beta_mat[, lambda_idx]            # named numeric vector
nonzero  <- b_vec[b_vec != 0]                 # keep only non-zero

pam50full_df <- tibble(
  gene   = names(nonzero),
  weight = as.numeric(nonzero)
) %>%
  arrange(desc(abs(weight)))

message(sprintf("[%s] PAM50full genes (non-zero): %d", SCRIPT_NAME, nrow(pam50full_df)))
message(sprintf("[%s] Genes: %s", SCRIPT_NAME, paste(pam50full_df$gene, collapse = ", ")))

# Save PAM50full weights
write_csv(pam50full_df, out_weights)
h_w  <- sha256_file(out_weights)
sz_w <- file.info(out_weights)$size / 1e6
registry_append("ALL", "PAM50full_weights", out_weights, h_w, "ok",
                SCRIPT_NAME, sz_w)
message(sprintf("[%s] PAM50full weights saved: %s (%.3f MB | SHA256: %s)",
                SCRIPT_NAME, out_weights, sz_w, h_w))

pam50full_weights <- setNames(pam50full_df$weight, pam50full_df$gene)

# =============================================================================
# STEP 2: Load CorePAM weights
# =============================================================================
corepam_df      <- strict_csv(file.path(PATHS$results$corepam, "CorePAM_weights.csv"))
corepam_weights <- setNames(corepam_df$weight, corepam_df$gene)
message(sprintf("[%s] CorePAM genes: %d", SCRIPT_NAME, length(corepam_weights)))

# =============================================================================
# STEP 3: Compute scores and compare per cohort
# =============================================================================
message(sprintf("[%s] === Step 3: Compute scores per cohort ===", SCRIPT_NAME))

cohort_config <- list(
  SCANB     = list(time_col = "os_time_months",  event_col = "os_event",  label = "SCAN-B",    endpoint = "OS"),
  TCGA_BRCA = list(time_col = "os_time_months",  event_col = "os_event",  label = "TCGA-BRCA", endpoint = "OS"),
  METABRIC  = list(time_col = "dss_time_months", event_col = "dss_event", label = "METABRIC",  endpoint = "DSS"),
  GSE20685  = list(time_col = "os_time_months",  event_col = "os_event",  label = "GSE20685",  endpoint = "OS")
)

scatter_data_all <- list()
comparison_rows  <- list()

for (cohort in names(cohort_config)) {
  cfg <- cohort_config[[cohort]]
  message(sprintf("[%s] Processing cohort: %s | endpoint: %s", SCRIPT_NAME, cohort, cfg$endpoint))

  # Load expression (gene x samples parquet)
  expr_path <- file.path(proc_cohort(cohort), "expression_genelevel_preZ.parquet")
  expr_mat  <- strict_parquet(expr_path)

  gene_col   <- "gene"
  sample_ids <- setdiff(names(expr_mat), gene_col)
  genes_expr <- expr_mat[[gene_col]]

  # Transpose to samples x genes matrix
  expr_vals          <- t(as.matrix(expr_mat[, sample_ids]))
  colnames(expr_vals) <- genes_expr

  # Intra-cohort Z-score per gene
  old_warn  <- getOption("warn"); options(warn = 0)
  gene_means <- colMeans(expr_vals, na.rm = TRUE)
  gene_sds   <- apply(expr_vals, 2, sd, na.rm = TRUE)
  options(warn = old_warn)

  zero_sd <- gene_sds == 0 | is.na(gene_sds)
  if (any(zero_sd)) {
    message(sprintf("[%s] %s: %d genes with SD=0 -> z=0 for those", SCRIPT_NAME, cohort, sum(zero_sd)))
  }
  gene_sds[zero_sd] <- 1

  z_mat          <- sweep(sweep(expr_vals, 2, gene_means, "-"), 2, gene_sds, "/")
  z_mat[, zero_sd] <- 0

  # --- CorePAM score ---
  corepam_res   <- compute_score(z_mat, corepam_weights)
  score_corepam <- corepam_res$score

  # --- PAM50full score ---
  pam50full_res  <- compute_score(z_mat, pam50full_weights)
  score_pam50full <- pam50full_res$score

  message(sprintf("[%s] %s CorePAM genes present: %d/%d (%.1f%%)",
                  SCRIPT_NAME, cohort, corepam_res$n_present,
                  length(corepam_weights), corepam_res$frac_present * 100))
  message(sprintf("[%s] %s PAM50full genes present: %d/%d (%.1f%%)",
                  SCRIPT_NAME, cohort, pam50full_res$n_present,
                  length(pam50full_weights), pam50full_res$frac_present * 100))

  # Standardize (z-score) for HR per 1 SD — determine direction via clinical data
  clin_path <- file.path(proc_cohort(cohort), "clinical_FINAL.parquet")
  clin_df   <- strict_parquet(clin_path)

  # Build score data frame with sample IDs
  score_tbl <- tibble(
    patient_id     = sample_ids,
    score_corepam  = score_corepam,
    score_pam50full = score_pam50full
  )

  # Join clinical — determine best key by actual overlap
  # GSE20685: expression uses GEO sample IDs (GSMxxxxxxx); clinical 'sample_id' matches, 'patient_id' does not
  # SCANB, TCGA_BRCA, METABRIC: 'patient_id' overlaps with expression sample IDs
  sample_ids_set <- sample_ids   # character vector of expression column names

  join_key <- "patient_id"  # default
  if ("patient_id" %in% names(clin_df)) {
    ov_pid <- length(intersect(clin_df$patient_id, sample_ids_set))
    if (ov_pid == 0 && "sample_id" %in% names(clin_df)) {
      ov_sid <- length(intersect(clin_df$sample_id, sample_ids_set))
      if (ov_sid > 0) {
        join_key <- "sample_id"
        message(sprintf("[%s] %s: joining on 'sample_id' (%d overlap)", SCRIPT_NAME, cohort, ov_sid))
      }
    } else {
      message(sprintf("[%s] %s: joining on 'patient_id' (%d overlap)", SCRIPT_NAME, cohort, ov_pid))
    }
  } else if ("sample_id" %in% names(clin_df)) {
    join_key <- "sample_id"
    message(sprintf("[%s] %s: no patient_id; joining on 'sample_id'", SCRIPT_NAME, cohort))
  }

  # score_tbl always has patient_id = sample_ids; rename if needed
  if (join_key == "sample_id") {
    score_tbl_join <- score_tbl %>% rename(sample_id = patient_id)
  } else {
    score_tbl_join <- score_tbl
  }

  clin_score <- inner_join(clin_df, score_tbl_join, by = join_key)
  # Ensure patient_id column exists in result for downstream use
  if (!"patient_id" %in% names(clin_score) && "sample_id" %in% names(clin_score)) {
    clin_score <- clin_score %>% mutate(patient_id = sample_id)
  }
  n_joined <- nrow(clin_score)
  message(sprintf("[%s] %s: joined %d samples", SCRIPT_NAME, cohort, n_joined))
  if (n_joined == 0) {
    message(sprintf("[%s] %s: ERROR - 0 samples after join; skipping cohort", SCRIPT_NAME, cohort))
    next
  }

  # Filter time <= 0
  time_col  <- cfg$time_col
  event_col <- cfg$event_col

  if (!(time_col %in% names(clin_score))) {
    message(sprintf("[%s] %s: time column '%s' not found; skipping", SCRIPT_NAME, cohort, time_col))
    next
  }
  if (!(event_col %in% names(clin_score))) {
    message(sprintf("[%s] %s: event column '%s' not found; skipping", SCRIPT_NAME, cohort, event_col))
    next
  }

  n_before <- nrow(clin_score)
  clin_score <- clin_score %>%
    filter(!is.na(.data[[time_col]]), .data[[time_col]] > 0,
           !is.na(.data[[event_col]]))
  n_after <- nrow(clin_score)
  if (n_before > n_after) {
    message(sprintf("[%s] %s: dropped %d rows with time<=0 or NA", SCRIPT_NAME, cohort, n_before - n_after))
  }
  if (n_after == 0) {
    message(sprintf("[%s] %s: 0 valid samples after time filter; skipping", SCRIPT_NAME, cohort))
    next
  }

  time_vec  <- clin_score[[time_col]]
  event_vec <- clin_score[[event_col]]

  # Standardize scores — check direction for each score independently
  old_warn <- getOption("warn"); options(warn = 0)

  # CorePAM direction (use analysis_ready score_direction)
  ready_path  <- file.path(proc_cohort(cohort), "analysis_ready.parquet")
  ready_df    <- strict_parquet(ready_path)
  # Use the stored score_direction from the previous pipeline step
  score_dir_corepam <- if ("score_direction" %in% names(ready_df)) {
    unique(ready_df$score_direction)[1]
  } else {
    "original"
  }
  options(warn = old_warn)

  # Apply same direction to our freshly computed CorePAM score
  corepam_sign <- if (score_dir_corepam == "inverted") -1 else 1

  old_warn2 <- getOption("warn"); options(warn = 0)
  score_corepam_z <- as.vector(scale(corepam_sign * clin_score$score_corepam))

  # PAM50full direction: check Cox HR
  score_pam50full_raw <- clin_score$score_pam50full
  old_warn3 <- getOption("warn"); options(warn = 0)
  score_pam50full_z_tmp <- as.vector(scale(score_pam50full_raw))
  cox_pam50full_dir <- tryCatch(
    coxph(Surv(time_vec, event_vec) ~ score_pam50full_z_tmp),
    error = function(e) NULL
  )
  options(warn = old_warn3)

  pam50full_sign <- 1
  if (!is.null(cox_pam50full_dir)) {
    hr_check <- exp(coef(cox_pam50full_dir)[1])
    if (hr_check < 1) {
      pam50full_sign <- -1
      message(sprintf("[%s] %s: PAM50full direction inverted (HR raw = %.4f)", SCRIPT_NAME, cohort, hr_check))
    } else {
      message(sprintf("[%s] %s: PAM50full direction original (HR raw = %.4f)", SCRIPT_NAME, cohort, hr_check))
    }
  }
  score_pam50full_z <- as.vector(scale(pam50full_sign * score_pam50full_raw))
  options(warn = old_warn2)

  # --- Cox HR for each score ---
  old_warn4 <- getOption("warn"); options(warn = 0)
  cox_cp <- tryCatch(
    coxph(Surv(time_vec, event_vec) ~ score_corepam_z),
    error = function(e) NULL
  )
  cox_pf <- tryCatch(
    coxph(Surv(time_vec, event_vec) ~ score_pam50full_z),
    error = function(e) NULL
  )
  options(warn = old_warn4)

  extract_cox <- function(fit) {
    if (is.null(fit)) return(list(hr = NA, lo = NA, hi = NA, p = NA, loghr = NA, se = NA))
    s   <- summary(fit)
    est <- s$coefficients[1, ]  # coef row
    ci  <- s$conf.int[1, ]
    list(
      hr    = ci["exp(coef)"],
      lo    = ci["lower .95"],
      hi    = ci["upper .95"],
      p     = est["Pr(>|z|)"],
      loghr = est["coef"],
      se    = est["se(coef)"]
    )
  }

  res_cp <- extract_cox(cox_cp)
  res_pf <- extract_cox(cox_pf)

  message(sprintf("[%s] %s CorePAM HR = %.4f (%.4f–%.4f), p = %.4g",
                  SCRIPT_NAME, cohort,
                  res_cp$hr, res_cp$lo, res_cp$hi, res_cp$p))
  message(sprintf("[%s] %s PAM50full HR = %.4f (%.4f–%.4f), p = %.4g",
                  SCRIPT_NAME, cohort,
                  res_pf$hr, res_pf$lo, res_pf$hi, res_pf$p))

  # --- C-index ---
  ci_cp <- bootstrap_cindex(time_vec, event_vec, score_corepam_z,  n_boot = 500)
  ci_pf <- bootstrap_cindex(time_vec, event_vec, score_pam50full_z, n_boot = 500)
  message(sprintf("[%s] %s CorePAM C-index = %.4f (%.4f–%.4f)",
                  SCRIPT_NAME, cohort, ci_cp$c_index, ci_cp$ci_low, ci_cp$ci_high))
  message(sprintf("[%s] %s PAM50full C-index = %.4f (%.4f–%.4f)",
                  SCRIPT_NAME, cohort, ci_pf$c_index, ci_pf$ci_low, ci_pf$ci_high))

  # --- Spearman correlation between scores ---
  rho <- cor(score_corepam_z, score_pam50full_z, method = "spearman", use = "complete.obs")
  message(sprintf("[%s] %s Spearman rho (CorePAM vs PAM50full) = %.4f", SCRIPT_NAME, cohort, rho))

  # Collect scatter data
  scatter_data_all[[cohort]] <- tibble(
    cohort           = cfg$label,
    corepam_z        = score_corepam_z,
    pam50full_z      = score_pam50full_z,
    endpoint         = cfg$endpoint
  )

  comparison_rows[[cohort]] <- tibble(
    cohort              = cfg$label,
    endpoint            = cfg$endpoint,
    n_samples           = nrow(clin_score),
    n_events            = sum(event_vec),
    spearman_rho        = rho,
    hr_corepam          = res_cp$hr,
    hr_corepam_lo95     = res_cp$lo,
    hr_corepam_hi95     = res_cp$hi,
    p_corepam           = res_cp$p,
    loghr_corepam       = res_cp$loghr,
    se_loghr_corepam    = res_cp$se,
    c_corepam           = ci_cp$c_index,
    c_corepam_lo95      = ci_cp$ci_low,
    c_corepam_hi95      = ci_cp$ci_high,
    hr_pam50full        = res_pf$hr,
    hr_pam50full_lo95   = res_pf$lo,
    hr_pam50full_hi95   = res_pf$hi,
    p_pam50full         = res_pf$p,
    loghr_pam50full     = res_pf$loghr,
    se_loghr_pam50full  = res_pf$se,
    c_pam50full         = ci_pf$c_index,
    c_pam50full_lo95    = ci_pf$ci_low,
    c_pam50full_hi95    = ci_pf$ci_high,
    genes_corepam       = corepam_res$n_present,
    genes_pam50full     = pam50full_res$n_present
  )
}

# Combine comparison results
comp_df <- bind_rows(comparison_rows)

# Save comparison CSV
write_csv(comp_df, out_comp_csv)
h_c  <- sha256_file(out_comp_csv)
sz_c <- file.info(out_comp_csv)$size / 1e6
registry_append("ALL", "pam50full_comparison", out_comp_csv, h_c, "ok",
                SCRIPT_NAME, sz_c)
message(sprintf("[%s] Comparison CSV saved: %s (SHA256: %s)", SCRIPT_NAME, out_comp_csv, h_c))

# =============================================================================
# STEP 4: Figures
# =============================================================================
message(sprintf("[%s] === Step 4: Generating figures ===", SCRIPT_NAME))

scatter_df <- bind_rows(scatter_data_all)

# Build annotation label per cohort
rho_labels <- comp_df %>%
  transmute(
    cohort = cohort,
    label  = sprintf("rho = %.3f", spearman_rho)
  )

scatter_df <- left_join(scatter_df, rho_labels, by = "cohort")

# Ordered factor so facets appear consistently
cohort_order <- c("SCAN-B", "TCGA-BRCA", "METABRIC", "GSE20685")
scatter_df$cohort <- factor(scatter_df$cohort, levels = cohort_order)

# --------------------------------------------------------------------------
# Figure A: Scatter (2x2 facets)
# --------------------------------------------------------------------------
message(sprintf("[%s] Figure A: Scatter CorePAM vs PAM50full", SCRIPT_NAME))

p_scatter <- ggplot(scatter_df, aes(x = pam50full_z, y = corepam_z)) +
  geom_point(size = 0.4, alpha = 0.35, color = "steelblue4") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = "firebrick3") +
  geom_text(
    data = scatter_df %>% distinct(cohort, label),
    aes(x = -Inf, y = Inf, label = label),
    hjust = -0.1, vjust = 1.3, size = 3.2, fontface = "italic", inherit.aes = FALSE
  ) +
  facet_wrap(~ cohort, nrow = 2, ncol = 2) +
  labs(
    title    = "CorePAM vs PAM50-full score",
    subtitle = sprintf("PAM50-full: %d genes (max df model) | CorePAM: %d genes",
                       nrow(pam50full_df), nrow(corepam_df)),
    x        = "PAM50-full score (z-standardized)",
    y        = "CorePAM score (z-standardized)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey90"),
    strip.text       = element_text(face = "bold"),
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

save_fig(p_scatter, "FigS_CorePAM_vs_PAM50full_EN", w = 9, h = 8)

# --------------------------------------------------------------------------
# Figure B: Forest plot — HR per 1 SD comparison
# --------------------------------------------------------------------------
message(sprintf("[%s] Figure B: Forest HR comparison CorePAM vs PAM50full", SCRIPT_NAME))

# Build pooled estimate via fixed-effects meta-analysis (inverse-variance weighting)
pool_row <- function(loghr_col, se_col, label_col = "cohort") {
  rows <- comp_df %>% select(!!sym(label_col), all_of(c(loghr_col, se_col))) %>%
    filter(!is.na(.data[[se_col]]), .data[[se_col]] > 0)
  w  <- 1 / rows[[se_col]]^2
  lw <- sum(w * rows[[loghr_col]]) / sum(w)
  se <- 1 / sqrt(sum(w))
  tibble(
    cohort  = "Pooled (FE)",
    loghr   = lw,
    se_loghr = se,
    hr      = exp(lw),
    hr_lo   = exp(lw - 1.96 * se),
    hr_hi   = exp(lw + 1.96 * se)
  )
}

forest_df_cp <- comp_df %>%
  transmute(
    cohort   = cohort,
    model    = "CorePAM",
    loghr    = loghr_corepam,
    se_loghr = se_loghr_corepam,
    hr       = hr_corepam,
    hr_lo    = hr_corepam_lo95,
    hr_hi    = hr_corepam_hi95,
    p_val    = p_corepam,
    n        = n_samples,
    endpoint = endpoint
  )

forest_df_pf <- comp_df %>%
  transmute(
    cohort   = cohort,
    model    = "PAM50-full",
    loghr    = loghr_pam50full,
    se_loghr = se_loghr_pam50full,
    hr       = hr_pam50full,
    hr_lo    = hr_pam50full_lo95,
    hr_hi    = hr_pam50full_hi95,
    p_val    = p_pam50full,
    n        = n_samples,
    endpoint = endpoint
  )

# Add pooled row for each model
pool_cp <- pool_row("loghr_corepam",  "se_loghr_corepam")  %>%
  mutate(model = "CorePAM",   p_val = NA, n = sum(comp_df$n_samples), endpoint = "pooled")
pool_pf <- pool_row("loghr_pam50full","se_loghr_pam50full") %>%
  mutate(model = "PAM50-full",p_val = NA, n = sum(comp_df$n_samples), endpoint = "pooled")

forest_df <- bind_rows(forest_df_cp, forest_df_pf, pool_cp, pool_pf)

# Ordered factor
cohort_levels <- c(cohort_order, "Pooled (FE)")
forest_df$cohort <- factor(forest_df$cohort, levels = rev(cohort_levels))
forest_df$model  <- factor(forest_df$model, levels = c("CorePAM", "PAM50-full"))

# Label: n and endpoint
forest_df <- forest_df %>%
  mutate(
    cohort_label = case_when(
      cohort == "Pooled (FE)" ~ "Pooled (FE)",
      TRUE ~ sprintf("%s\n(n=%d, %s)", cohort, n, endpoint)
    )
  )
forest_df$cohort_label <- factor(forest_df$cohort_label,
                                  levels = rev(c(
                                    sprintf("SCAN-B\n(n=%d, OS)",     comp_df$n_samples[comp_df$cohort == "SCAN-B"]),
                                    sprintf("TCGA-BRCA\n(n=%d, OS)",  comp_df$n_samples[comp_df$cohort == "TCGA-BRCA"]),
                                    sprintf("METABRIC\n(n=%d, DSS)",  comp_df$n_samples[comp_df$cohort == "METABRIC"]),
                                    sprintf("GSE20685\n(n=%d, OS)",   comp_df$n_samples[comp_df$cohort == "GSE20685"]),
                                    "Pooled (FE)"
                                  )))

# Dodge offset for side-by-side
dodge_width <- 0.5

p_forest <- ggplot(
  forest_df,
  aes(x = hr, y = cohort_label, color = model, shape = model)
) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.7) +
  geom_pointrange(
    aes(xmin = hr_lo, xmax = hr_hi),
    position = position_dodge(width = dodge_width),
    linewidth = 0.7,
    size = 0.5
  ) +
  scale_color_manual(
    name   = "Score",
    values = c("CorePAM" = "#2166AC", "PAM50-full" = "#888888")
  ) +
  scale_shape_manual(
    name   = "Score",
    values = c("CorePAM" = 16, "PAM50-full" = 17)
  ) +
  scale_x_log10(
    breaks = c(0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0),
    labels = c("0.5","0.8","1.0","1.2","1.5","2.0","2.5","3.0")
  ) +
  labs(
    title    = "HR per 1 SD: CorePAM vs PAM50-full",
    subtitle = sprintf("CorePAM: %d genes | PAM50-full: %d genes | Univariate Cox",
                       nrow(corepam_df), nrow(pam50full_df)),
    x        = "Hazard ratio per 1 SD (95% CI, log scale)",
    y        = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(linetype = "dotted", color = "grey80"),
    axis.text.y      = element_text(size = 9)
  )

save_fig(p_forest, "FigS_CorePAM_vs_PAM50full_HR_EN", w = 9, h = 6)

# =============================================================================
# DONE
# =============================================================================
message(sprintf("[%s] === All done ===", SCRIPT_NAME))
message(sprintf("[%s] PAM50full genes: %d", SCRIPT_NAME, nrow(pam50full_df)))
message(sprintf("[%s] Cohorts processed: %s", SCRIPT_NAME, paste(names(cohort_config), collapse = ", ")))
message(sprintf("[%s] Outputs:", SCRIPT_NAME))
message(sprintf("[%s]   %s", SCRIPT_NAME, out_weights))
message(sprintf("[%s]   %s", SCRIPT_NAME, out_comp_csv))
message(sprintf("[%s]   %s", SCRIPT_NAME, file.path(fig_supp, "FigS_CorePAM_vs_PAM50full_EN.pdf")))
message(sprintf("[%s]   %s", SCRIPT_NAME, file.path(fig_supp, "FigS_CorePAM_vs_PAM50full_HR_EN.pdf")))
