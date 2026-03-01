# =============================================================================
# SCRIPT: 08_meta_survival.R
# PURPOSE: Random-effects meta-analysis of univariate and CORE-A log(HR) per
#          cohort; I2/tau2; leave-one-out (sensitivity); forest plot Fig4.
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

COHORTS <- c("SCANB", "TCGA_BRCA", "METABRIC", "GSE20685")

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
message(sprintf("[%s] Cohorts with valid logHR+SE: %d", SCRIPT_NAME, nrow(all_res)))

if (nrow(all_res) < 2) {
  stop(sprintf("[%s] Fewer than 2 cohorts with valid data. Cannot perform meta-analysis.", SCRIPT_NAME))
}

# --------------------------------------------------------------------------
# 2) Random-effects meta-analysis — univariate (logHR)
# --------------------------------------------------------------------------
message(sprintf("[%s] Random-effects meta-analysis (uni)...", SCRIPT_NAME))

old_warn <- getOption("warn"); options(warn = 0)
meta_uni <- tryCatch(
  rma(yi = loghr_uni, sei = se_loghr_uni, data = all_res, method = "REML"),
  error = function(e) {
    message(sprintf("[%s] rma REML failed; trying DL: %s", SCRIPT_NAME, e$message))
    tryCatch(
      rma(yi = loghr_uni, sei = se_loghr_uni, data = all_res, method = "DL"),
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

# --------------------------------------------------------------------------
# 3) Random-effects meta-analysis — multivariate CORE-A (when available)
# --------------------------------------------------------------------------
has_multi <- "loghr_multi" %in% names(all_res) && "se_loghr_multi" %in% names(all_res)

# Compute loghr_multi / se_loghr_multi from hr_multi if not already present
if (!has_multi && all(c("hr_multi", "hr_multi_lo95", "hr_multi_hi95") %in% names(all_res))) {
  all_res <- all_res %>%
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
  df_m <- all_res[!is.na(all_res$loghr_multi) & !is.na(all_res$se_loghr_multi), ]
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
  loo_df$cohort_left_out <- all_res$cohort
  message(sprintf("[%s] LOO: %d analyses", SCRIPT_NAME, nrow(loo_df)))
}

# --------------------------------------------------------------------------
# 5) Save meta-analysis results
# --------------------------------------------------------------------------
main_dir <- PATHS$results$main
dir.create(main_dir, showWarnings = FALSE, recursive = TRUE)

meta_summary <- tibble(
  analysis        = c("univariate", if (!is.null(meta_multi)) "CORE-A_multivariate" else character(0)),
  n_cohorts       = c(nrow(all_res), if (!is.null(meta_multi)) nrow(df_m) else integer(0)),
  meta_logHR      = c(round(meta_uni$b[1], 6),
                      if (!is.null(meta_multi)) round(meta_multi$b[1], 6) else numeric(0)),
  meta_HR         = c(round(exp(meta_uni$b[1]), 4),
                      if (!is.null(meta_multi)) round(exp(meta_multi$b[1]), 4) else numeric(0)),
  meta_HR_lo95    = c(round(exp(meta_uni$ci.lb), 4),
                      if (!is.null(meta_multi)) round(exp(meta_multi$ci.lb), 4) else numeric(0)),
  meta_HR_hi95    = c(round(exp(meta_uni$ci.ub), 4),
                      if (!is.null(meta_multi)) round(exp(meta_multi$ci.ub), 4) else numeric(0)),
  p_value         = c(signif(meta_uni$pval, 4),
                      if (!is.null(meta_multi)) signif(meta_multi$pval, 4) else numeric(0)),
  I2_pct          = c(round(meta_uni$I2, 2),
                      if (!is.null(meta_multi)) round(meta_multi$I2, 2) else numeric(0)),
  tau2            = c(round(meta_uni$tau2, 6),
                      if (!is.null(meta_multi)) round(meta_multi$tau2, 6) else numeric(0)),
  Q_stat          = c(round(meta_uni$QE, 4),
                      if (!is.null(meta_multi)) round(meta_multi$QE, 4) else numeric(0)),
  Q_p             = c(signif(meta_uni$QEp, 4),
                      if (!is.null(meta_multi)) signif(meta_multi$QEp, 4) else numeric(0))
)

meta_path <- file.path(main_dir, "meta_survival_summary.csv")
readr::write_csv(meta_summary, meta_path)
h_meta <- sha256_file(meta_path)
registry_append("ALL", "meta_survival_summary", meta_path, h_meta, "ok", SCRIPT_NAME,
                file.info(meta_path)$size / 1e6)

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

# Prepare data for manual forest plot with ggplot2
plot_df <- all_res %>%
  select(cohort, hr_uni, hr_uni_lo95, hr_uni_hi95, n_samples, n_events) %>%
  mutate(
    label = sprintf("%s (n=%d, ev=%d)", cohort,
                    ifelse("n_samples" %in% names(all_res), n_samples, NA_integer_),
                    ifelse("n_events" %in% names(all_res), n_events, NA_integer_)),
    cohort_f = factor(cohort, levels = rev(COHORTS))
  )

# Add meta-analysis summary row
meta_row <- tibble(
  cohort    = "Meta (RE)",
  hr_uni    = exp(meta_uni$b[1]),
  hr_uni_lo95 = exp(meta_uni$ci.lb),
  hr_uni_hi95 = exp(meta_uni$ci.ub),
  n_samples = NA_integer_,
  n_events  = NA_integer_,
  label     = sprintf("Meta RE (I2=%.0f%%, tau2=%.3f)", meta_uni$I2, meta_uni$tau2),
  cohort_f  = factor("Meta (RE)", levels = c(rev(COHORTS), "Meta (RE)"))
)

levels_full <- c(rev(COHORTS), "Meta (RE)")
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
