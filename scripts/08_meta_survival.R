# =============================================================================
# SCRIPT: 08_meta_survival.R
# PURPOSE: Random-effects meta-analysis of univariate and CORE-A log(HR) per
#          cohort; I2/tau2; leave-one-out (sensitivity); forest plot Fig4.
#          Cohorts: SCANB (training), TCGA_BRCA, METABRIC, GSE20685, GSE1456.
#          Meta-analysis restricted to validation cohorts only (K=4).
#          Follows Memorial v6.1 sec.8.2.
# PROJETO: Core-PAM (Memorial v6.1 / Freeze Core-PAM)
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "08_meta_survival.R"

suppressPackageStartupMessages({
  library(metafor)
  library(ggplot2)
})

message(sprintf("[%s] Starting CorePAM survival meta-analysis", SCRIPT_NAME))

COHORTS <- c("SCANB", "TCGA_BRCA", "METABRIC", "GSE20685", "GSE1456")

# --------------------------------------------------------------------------
# 1) Read survival results from all cohorts
# --------------------------------------------------------------------------
supp_dir <- PATHS$results$supp

res_list <- lapply(COHORTS, function(coh) {
  path <- file.path(supp_dir, sprintf("survival_results_%s.csv", coh))
  if (!file.exists(path)) {
    message(sprintf("[%s] WARNING: file not found for %s: %s", SCRIPT_NAME, coh, path))
    return(NULL)
  }
  df <- strict_csv(path)
  df$cohort_label <- coh
  df
})

# Filter cohorts with data
res_list <- Filter(Negate(is.null), res_list)
if (length(res_list) == 0) {
  stop(sprintf("[%s] No survival_results files found. Run scripts 07 first.", SCRIPT_NAME))
}

all_res <- bind_rows(res_list)
message(sprintf("[%s] Cohorts with results: %d", SCRIPT_NAME, nrow(all_res)))

# Check required columns
needed_cols <- c("cohort", "loghr_uni", "se_loghr_uni", "hr_uni", "hr_uni_lo95", "hr_uni_hi95")
missing_cols <- setdiff(needed_cols, names(all_res))
if (length(missing_cols) > 0) {
  stop(sprintf("[%s] Columns missing in results: %s",
               SCRIPT_NAME, paste(missing_cols, collapse = ", ")))
}

# Remove rows with missing logHR or SE
all_res <- all_res[!is.na(all_res$loghr_uni) & !is.na(all_res$se_loghr_uni), ]
message(sprintf("[%s] All cohorts with valid logHR+SE: %d rows", SCRIPT_NAME, nrow(all_res)))

# --------------------------------------------------------------------------
# Subset: EXTERNAL VALIDATION only (exclude SCAN-B training cohort)
# Primary endpoint: METABRIC DSS (primary); TCGA OS; GSE20685 OS
# --------------------------------------------------------------------------
VALIDATION_COHORTS <- c("TCGA_BRCA", "METABRIC", "GSE20685", "GSE1456")

# For METABRIC: keep only the DSS row (primary endpoint) for the primary analysis
val_res <- all_res %>%
  filter(cohort %in% VALIDATION_COHORTS) %>%
  filter(!(cohort == "METABRIC" & !is.null(endpoint) & "endpoint" %in% names(.) & endpoint != "DSS"))

# If the survival CSV has an 'endpoint' column (METABRIC has DSS + OS rows),
# keep only the primary endpoint row per cohort for the primary analysis
if ("endpoint" %in% names(val_res)) {
  val_res <- val_res %>%
    filter(
      (cohort == "METABRIC"   & endpoint == "DSS") |
      (cohort == "TCGA_BRCA"  & endpoint == "OS")  |
      (cohort == "GSE20685"   & endpoint == "OS")  |
      (cohort == "GSE1456"    & endpoint == "OS")  |
      !(cohort %in% c("METABRIC", "TCGA_BRCA", "GSE20685", "GSE1456"))
    )
}

# Deduplicate — keep one row per cohort (the first, which is the primary endpoint)
val_res <- val_res %>% group_by(cohort) %>% slice(1) %>% ungroup()

message(sprintf("[%s] Validation cohorts for primary meta-analysis: %d (K=%d)",
                SCRIPT_NAME, nrow(val_res), length(unique(val_res$cohort))))

if (nrow(val_res) < 2) {
  stop(sprintf("[%s] Fewer than 2 validation cohorts with valid data. Cannot perform meta-analysis.", SCRIPT_NAME))
}

# --------------------------------------------------------------------------
# 2) Random-effects meta-analysis — univariate (logHR), validation only
# --------------------------------------------------------------------------
message(sprintf("[%s] Random-effects meta-analysis (uni), K=%d validation cohorts...",
                SCRIPT_NAME, nrow(val_res)))

old_warn <- getOption("warn"); options(warn = 0)
meta_uni <- tryCatch(
  rma(yi = loghr_uni, sei = se_loghr_uni, data = val_res, method = "REML"),
  error = function(e) {
    message(sprintf("[%s] rma REML failed; trying DL: %s", SCRIPT_NAME, e$message))
    tryCatch(
      rma(yi = loghr_uni, sei = se_loghr_uni, data = val_res, method = "DL"),
      error = function(e2) { stop(sprintf("[%s] Meta-analysis failed: %s", SCRIPT_NAME, e2$message)) }
    )
  }
)
options(warn = old_warn)

message(sprintf("[%s] Meta HR uni = %.3f (%.3f-%.3f) | I2=%.1f%% | tau2=%.4f | p=%.4g",
                SCRIPT_NAME,
                exp(meta_uni$b[1]),
                exp(meta_uni$ci.lb),
                exp(meta_uni$ci.ub),
                meta_uni$I2,
                meta_uni$tau2,
                meta_uni$pval))

# Hartung-Knapp (HK) sensitivity — more conservative CIs when I² > 0
old_warn <- getOption("warn"); options(warn = 0)
meta_uni_hk <- tryCatch(
  rma(yi = loghr_uni, sei = se_loghr_uni, data = val_res, method = "REML", test = "knha"),
  error = function(e) {
    message(sprintf("[%s] HK REML failed; trying DL: %s", SCRIPT_NAME, e$message))
    tryCatch(
      rma(yi = loghr_uni, sei = se_loghr_uni, data = val_res, method = "DL", test = "knha"),
      error = function(e2) { message(sprintf("[%s] HK failed: %s", SCRIPT_NAME, e2$message)); NULL }
    )
  }
)
options(warn = old_warn)

if (!is.null(meta_uni_hk)) {
  message(sprintf("[%s] HK sensitivity HR = %.3f (%.3f-%.3f) | p=%.4g",
                  SCRIPT_NAME,
                  exp(meta_uni_hk$b[1]),
                  exp(meta_uni_hk$ci.lb),
                  exp(meta_uni_hk$ci.ub),
                  meta_uni_hk$pval))
}

# --------------------------------------------------------------------------
# 3) Random-effects meta-analysis — multivariate CORE-A (when available)
# --------------------------------------------------------------------------
has_multi <- "loghr_multi" %in% names(val_res) && "se_loghr_multi" %in% names(val_res)

# Compute loghr_multi / se_loghr_multi from hr_multi if not already present
if (!has_multi && all(c("hr_multi", "hr_multi_lo95", "hr_multi_hi95") %in% names(val_res))) {
  val_res <- val_res %>%
    mutate(
      loghr_multi    = ifelse(!is.na(hr_multi), log(hr_multi), NA_real_),
      se_loghr_multi = ifelse(!is.na(hr_multi_lo95) & !is.na(hr_multi_hi95),
                              (log(hr_multi_hi95) - log(hr_multi_lo95)) / (2 * 1.96),
                              NA_real_)
    )
  has_multi <- TRUE
}

meta_multi <- NULL
if (has_multi) {
  df_m <- val_res[!is.na(val_res$loghr_multi) & !is.na(val_res$se_loghr_multi), ]
  if (nrow(df_m) >= 2) {
    old_warn <- getOption("warn"); options(warn = 0)
    meta_multi <- tryCatch(
      rma(yi = loghr_multi, sei = se_loghr_multi, data = df_m, method = "REML"),
      error = function(e) {
        tryCatch(
          rma(yi = loghr_multi, sei = se_loghr_multi, data = df_m, method = "DL"),
          error = function(e2) NULL
        )
      }
    )
    options(warn = old_warn)
    if (!is.null(meta_multi)) {
      message(sprintf("[%s] Meta HR multi (CORE-A) = %.3f (%.3f-%.3f) | I2=%.1f%% | p=%.4g",
                      SCRIPT_NAME,
                      exp(meta_multi$b[1]),
                      exp(meta_multi$ci.lb),
                      exp(meta_multi$ci.ub),
                      meta_multi$I2,
                      meta_multi$pval))
    }
  }
}

# --------------------------------------------------------------------------
# 4) Leave-one-out (sensitivity)
# --------------------------------------------------------------------------
message(sprintf("[%s] Leave-one-out...", SCRIPT_NAME))

old_warn <- getOption("warn"); options(warn = 0)
loo_uni <- tryCatch(
  leave1out(meta_uni),
  error = function(e) { message("[", SCRIPT_NAME, "] LOO failed:", e$message); NULL }
)
options(warn = old_warn)

loo_df <- NULL
if (!is.null(loo_uni)) {
  loo_df <- as.data.frame(loo_uni)
  loo_df$cohort_left_out <- val_res$cohort
  message(sprintf("[%s] LOO: %d analyses", SCRIPT_NAME, nrow(loo_df)))
}

# --------------------------------------------------------------------------
# 4a) OS-harmonised sensitivity meta-analysis
# Use OS endpoint for all four validation cohorts (replaces METABRIC DSS
# with METABRIC OS sensitivity results). Removes endpoint heterogeneity.
# --------------------------------------------------------------------------
message(sprintf("[%s] OS-harmonised sensitivity meta-analysis...", SCRIPT_NAME))

meta_os_harm <- NULL

# Get METABRIC OS sensitivity row from all_res (has all endpoints)
if ("endpoint" %in% names(all_res)) {
  os_harm_df <- all_res %>%
    filter(cohort %in% VALIDATION_COHORTS) %>%
    filter(
      (cohort == "METABRIC"   & endpoint == "OS")  |
      (cohort == "TCGA_BRCA"  & endpoint == "OS")  |
      (cohort == "GSE20685"   & endpoint == "OS")  |
      (cohort == "GSE1456"    & endpoint == "OS")
    ) %>%
    group_by(cohort) %>% slice(1) %>% ungroup() %>%
    filter(!is.na(loghr_uni) & !is.na(se_loghr_uni))
} else {
  # If no endpoint column, use the same val_res rows (already primary endpoints)
  os_harm_df <- val_res
}

if (nrow(os_harm_df) >= 2) {
  old_warn <- getOption("warn"); options(warn = 0)
  meta_os_harm <- tryCatch(
    rma(yi = loghr_uni, sei = se_loghr_uni, data = os_harm_df, method = "REML"),
    error = function(e) {
      tryCatch(
        rma(yi = loghr_uni, sei = se_loghr_uni, data = os_harm_df, method = "DL"),
        error = function(e2) {
          message(sprintf("[%s] OS-harm meta failed: %s", SCRIPT_NAME, e2$message))
          NULL
        }
      )
    }
  )
  options(warn = old_warn)
  if (!is.null(meta_os_harm)) {
    message(sprintf("[%s] OS-harm HR = %.3f (%.3f-%.3f) K=%d I2=%.1f%% p=%.4g",
                    SCRIPT_NAME,
                    exp(meta_os_harm$b[1]),
                    exp(meta_os_harm$ci.lb),
                    exp(meta_os_harm$ci.ub),
                    nrow(os_harm_df),
                    meta_os_harm$I2,
                    meta_os_harm$pval))
  }
} else {
  message(sprintf("[%s] OS-harmonised analysis: insufficient data (rows=%d)",
                  SCRIPT_NAME, nrow(os_harm_df)))
}

# --------------------------------------------------------------------------
# 5) Save meta-analysis results
# --------------------------------------------------------------------------
main_dir <- PATHS$results$main
dir.create(main_dir, showWarnings = FALSE, recursive = TRUE)

# Helper: build a single row for meta_summary tibble (full precision)
make_row <- function(analysis_name, rma_obj, k) {
  tibble(
    analysis     = analysis_name,
    n_cohorts    = k,
    meta_logHR   = as.numeric(rma_obj$b[1]),
    meta_HR      = exp(as.numeric(rma_obj$b[1])),
    meta_HR_lo95 = exp(rma_obj$ci.lb),
    meta_HR_hi95 = exp(rma_obj$ci.ub),
    p_value      = rma_obj$pval,
    I2_pct       = rma_obj$I2,
    tau2         = rma_obj$tau2,
    Q_stat       = rma_obj$QE,
    Q_p          = rma_obj$QEp
  )
}

meta_summary <- bind_rows(
  make_row("primary_validation_re",    meta_uni,       nrow(val_res)),
  if (!is.null(meta_uni_hk))
    make_row("primary_validation_re_hk", meta_uni_hk,  nrow(val_res)) else NULL,
  if (!is.null(meta_os_harm))
    make_row("os_harmonized_re",         meta_os_harm, nrow(os_harm_df)) else NULL
)

meta_path <- file.path(main_dir, "meta_survival_summary.csv")
readr::write_csv(meta_summary, meta_path)
h_meta <- sha256_file(meta_path)
registry_append("ALL", "meta_survival_summary", meta_path, h_meta, "ok", SCRIPT_NAME,
                file.info(meta_path)$size / 1e6)

# Save full-precision .rds with metafor objects for exact reproducibility
meta_rds_path <- file.path(main_dir, "meta_survival_objects.rds")
saveRDS(list(
  meta_uni    = meta_uni,
  meta_uni_hk = meta_uni_hk,
  meta_multi  = meta_multi,
  meta_os_harm = meta_os_harm,
  val_res     = val_res
), meta_rds_path)
h_rds <- sha256_file(meta_rds_path)
registry_append("ALL", "meta_survival_objects_rds", meta_rds_path, h_rds, "ok", SCRIPT_NAME,
                file.info(meta_rds_path)$size / 1e6)

# LOO output (leave-one-out)
if (!is.null(loo_df)) {
  loo_path <- file.path(PATHS$results$supp, "meta_survival_loo.csv")
  readr::write_csv(loo_df, loo_path)
  h_loo <- sha256_file(loo_path)
  registry_append("ALL", "meta_survival_loo", loo_path, h_loo, "ok", SCRIPT_NAME,
                  file.info(loo_path)$size / 1e6)
}

# --------------------------------------------------------------------------
# 6) Forest plot — Fig4
# --------------------------------------------------------------------------
fp_pdf <- file.path(PATHS$figures$main_en_pdf, "Fig4_Meta_Forest_HR_per1SD_CorePAM.pdf")
fp_png <- file.path(PATHS$figures$main_en_png, "Fig4_Meta_Forest_HR_per1SD_CorePAM.png")

# Prepare data for manual forest plot with ggplot2 — validation cohorts only
plot_df <- val_res %>%
  select(cohort, hr_uni, hr_uni_lo95, hr_uni_hi95, n_samples, n_events) %>%
  mutate(
    label = sprintf("%s (n=%d, ev=%d)", cohort,
                    ifelse("n_samples" %in% names(val_res), n_samples, NA_integer_),
                    ifelse("n_events" %in% names(val_res), n_events, NA_integer_)),
    cohort_f = factor(cohort, levels = rev(VALIDATION_COHORTS))
  )

# Add meta-analysis summary row
meta_row <- tibble(
  cohort    = "Meta (RE)",
  hr_uni    = exp(meta_uni$b[1]),
  hr_uni_lo95 = exp(meta_uni$ci.lb),
  hr_uni_hi95 = exp(meta_uni$ci.ub),
  n_samples = NA_integer_,
  n_events  = NA_integer_,
  label     = sprintf("Meta RE (K=%d, I2=%.0f%%, tau2=%.3f)",
                      nrow(val_res), meta_uni$I2, meta_uni$tau2),
  cohort_f  = factor("Meta (RE)", levels = c(rev(VALIDATION_COHORTS), "Meta (RE)"))
)

levels_full <- c(rev(VALIDATION_COHORTS), "Meta (RE)")
plot_df$cohort_f <- factor(plot_df$cohort, levels = levels_full)
meta_row$cohort_f <- factor("Meta (RE)", levels = levels_full)

fp_data <- bind_rows(plot_df, meta_row)

fp <- ggplot(fp_data, aes(x = hr_uni, y = cohort_f)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(xmin = hr_uni_lo95, xmax = hr_uni_hi95),
                orientation = "y", width = 0.2, linewidth = 0.8) +
  geom_point(aes(shape = ifelse(cohort == "Meta (RE)", 18L, 15L),
                 size  = ifelse(cohort == "Meta (RE)", 5, 3),
                 color = ifelse(cohort == "Meta (RE)", "Meta (RE)", "Cohort"))) +
  scale_shape_identity() +
  scale_size_identity() +
  scale_color_manual(values = c("Cohort" = "#2980B9", "Meta (RE)" = "#E74C3C"),
                     name = NULL) +
  scale_x_log10(breaks = c(0.5, 1, 2, 3, 4)) +
  labs(
    title = "CorePAM Meta-analysis: HR per 1 SD (random-effects)",
    x     = "HR per 1 SD of score (log scale)",
    y     = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title  = element_text(face = "bold", size = 13),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom"
  )

old_warn <- getOption("warn"); options(warn = 0)
pdf(fp_pdf, width = 9, height = 5)
print(fp)
dev.off()
png(fp_png, width = 900, height = 500, res = 100)
print(fp)
dev.off()
options(warn = old_warn)

h_fp_pdf <- sha256_file(fp_pdf)
h_fp_png <- sha256_file(fp_png)
registry_append("ALL", "figure_forest_meta", fp_pdf, h_fp_pdf, "ok", SCRIPT_NAME,
                file.info(fp_pdf)$size / 1e6)
registry_append("ALL", "figure_forest_meta_png", fp_png, h_fp_png, "ok", SCRIPT_NAME,
                file.info(fp_png)$size / 1e6)

message(sprintf("[%s] COMPLETED | Meta HR=%.3f (%.3f-%.3f) I2=%.1f%%",
                SCRIPT_NAME,
                exp(meta_uni$b[1]),
                exp(meta_uni$ci.lb),
                exp(meta_uni$ci.ub),
                meta_uni$I2))
