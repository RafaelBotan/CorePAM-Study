# =============================================================================
# SCRIPT: 34_subgroup_er_stratified.R
# PURPOSE: ER-stratified subgroup analysis in SCAN-B (N=3069, ER available ~94%).
#          Computes Cox HR + C-index per ER subgroup (ER+ vs ER-),
#          tests interaction (score * ER), and generates forest plot.
#          Addresses reviewer request for subgroup analysis by receptor status.
# OUTPUT:
#   results/supp/er_stratified_analysis.csv
#   figures/supp/en/pdf/FigS_ER_Stratified_Forest.pdf
#   figures/supp/en/png/FigS_ER_Stratified_Forest.png
# PROJECT: Core-PAM (Memorial CorePAM — Major Revision)
# =============================================================================

source("scripts/00_setup.R")
source("scripts/00_colors.R")
SCRIPT_NAME <- "34_subgroup_er_stratified.R"

suppressPackageStartupMessages({
  library(survival)
  library(arrow)
  library(ggplot2)
})

set.seed(FREEZE$seed_folds)
message(sprintf("[%s] ========== ER-STRATIFIED ANALYSIS ==========", SCRIPT_NAME))

# --------------------------------------------------------------------------
# 1) Load SCAN-B analysis_ready
# --------------------------------------------------------------------------
ar_path <- file.path(proc_cohort("SCANB"), "analysis_ready.parquet")
df <- as.data.frame(strict_parquet(ar_path))

time_col  <- "os_time_months"
event_col <- "os_event"

stopifnot(all(c(time_col, event_col, "score_z", "er_status") %in% names(df)))

df <- df[!is.na(df[[time_col]]) & df[[time_col]] > 0 &
           !is.na(df[[event_col]]) & !is.na(df$score_z) &
           !is.na(df$er_status), ]

message(sprintf("[%s] SCAN-B N=%d with ER status | ER+: %d | ER-: %d",
                SCRIPT_NAME, nrow(df),
                sum(df$er_status == 1 | df$er_status == "positive" | df$er_status == "ER+"),
                sum(df$er_status == 0 | df$er_status == "negative" | df$er_status == "ER-")))

# Standardize ER status to binary (handle different encodings)
if (is.numeric(df$er_status)) {
  df$er_binary <- ifelse(df$er_status == 1, "ER+", "ER-")
} else {
  df$er_binary <- ifelse(
    tolower(df$er_status) %in% c("positive", "er+", "1"),
    "ER+", "ER-"
  )
}

# --------------------------------------------------------------------------
# 2) Bootstrap C-index helper
# --------------------------------------------------------------------------
bootstrap_cindex <- function(time, event, score_z, n_boot = 1000, seed = 42) {
  set.seed(seed)
  n <- length(time)
  cvals <- numeric(n_boot)
  old_warn <- getOption("warn"); options(warn = 0)
  for (i in seq_len(n_boot)) {
    idx <- sample(n, n, replace = TRUE)
    cx <- tryCatch(
      concordance(Surv(time[idx], event[idx]) ~ score_z[idx])$concordance,
      error = function(e) NA_real_
    )
    cvals[i] <- max(cx, 1 - cx, na.rm = TRUE)
  }
  options(warn = old_warn)
  c_raw <- concordance(Surv(time, event) ~ score_z)$concordance
  list(
    c_index = max(c_raw, 1 - c_raw),
    ci_low  = quantile(cvals, 0.025, na.rm = TRUE),
    ci_high = quantile(cvals, 0.975, na.rm = TRUE)
  )
}

# --------------------------------------------------------------------------
# 3) Overall analysis (for reference)
# --------------------------------------------------------------------------
old_warn <- getOption("warn"); options(warn = 0)
cox_overall <- coxph(Surv(df[[time_col]], df[[event_col]]) ~ score_z, data = df)
options(warn = old_warn)
sm_all <- summary(cox_overall)

boot_all <- bootstrap_cindex(
  df[[time_col]], df[[event_col]], df$score_z,
  n_boot = as.integer(FREEZE$bootstrap_n),
  seed   = as.integer(FREEZE$seed_folds)
)

# --------------------------------------------------------------------------
# 4) Stratified analysis: ER+ and ER- separately
# --------------------------------------------------------------------------
results <- list()

for (er_group in c("ER+", "ER-")) {
  df_sub <- df[df$er_binary == er_group, ]
  n_sub  <- nrow(df_sub)
  n_ev   <- sum(df_sub[[event_col]])

  message(sprintf("[%s] %s: N=%d, Events=%d", SCRIPT_NAME, er_group, n_sub, n_ev))

  if (n_sub < 30 || n_ev < 10) {
    message(sprintf("[%s] %s: insufficient sample size, skipping", SCRIPT_NAME, er_group))
    results[[er_group]] <- data.frame(
      subgroup = er_group, n = n_sub, n_events = n_ev,
      hr = NA, hr_lo95 = NA, hr_hi95 = NA, p_value = NA,
      c_index = NA, c_lo95 = NA, c_hi95 = NA,
      stringsAsFactors = FALSE
    )
    next
  }

  # Cox HR
  old_warn <- getOption("warn"); options(warn = 0)
  cox_sub <- tryCatch(
    coxph(Surv(df_sub[[time_col]], df_sub[[event_col]]) ~ score_z, data = df_sub),
    error = function(e) NULL
  )
  options(warn = old_warn)

  hr <- lo <- hi <- p_val <- NA_real_
  if (!is.null(cox_sub)) {
    sm <- summary(cox_sub)
    hr    <- sm$conf.int[1, "exp(coef)"]
    lo    <- sm$conf.int[1, "lower .95"]
    hi    <- sm$conf.int[1, "upper .95"]
    p_val <- sm$coefficients[1, "Pr(>|z|)"]
  }

  # Bootstrap C-index
  boot_sub <- bootstrap_cindex(
    df_sub[[time_col]], df_sub[[event_col]], df_sub$score_z,
    n_boot = as.integer(FREEZE$bootstrap_n),
    seed   = as.integer(FREEZE$seed_folds)
  )

  results[[er_group]] <- data.frame(
    subgroup  = er_group,
    n         = n_sub,
    n_events  = n_ev,
    hr        = round(hr, 4),
    hr_lo95   = round(lo, 4),
    hr_hi95   = round(hi, 4),
    p_value   = signif(p_val, 4),
    c_index   = round(boot_sub$c_index, 4),
    c_lo95    = round(boot_sub$ci_low, 4),
    c_hi95    = round(boot_sub$ci_high, 4),
    stringsAsFactors = FALSE
  )

  message(sprintf("[%s] %s: HR=%.3f (%.3f-%.3f) p=%.4g | C=%.4f",
                  SCRIPT_NAME, er_group, hr, lo, hi, p_val, boot_sub$c_index))
}

# --------------------------------------------------------------------------
# 5) Interaction test: score_z * er_binary
# --------------------------------------------------------------------------
df$er_numeric <- ifelse(df$er_binary == "ER+", 1, 0)

old_warn <- getOption("warn"); options(warn = 0)
cox_interaction <- tryCatch(
  coxph(Surv(df[[time_col]], df[[event_col]]) ~ score_z * er_numeric, data = df),
  error = function(e) NULL
)
options(warn = old_warn)

interaction_hr <- interaction_p <- NA_real_
if (!is.null(cox_interaction)) {
  sm_int <- summary(cox_interaction)
  int_idx <- grep("score_z:er_numeric", rownames(sm_int$coefficients))
  if (length(int_idx) > 0) {
    interaction_hr <- round(exp(sm_int$coefficients[int_idx, "coef"]), 4)
    interaction_p  <- signif(sm_int$coefficients[int_idx, "Pr(>|z|)"], 4)
    message(sprintf("[%s] Interaction HR=%.4f, p=%.4g", SCRIPT_NAME, interaction_hr, interaction_p))
  }
}

# --------------------------------------------------------------------------
# 6) Compile results table
# --------------------------------------------------------------------------
res_df <- do.call(rbind, results)

# Add overall row
res_df <- rbind(
  data.frame(
    subgroup = "Overall",
    n        = nrow(df),
    n_events = sum(df[[event_col]]),
    hr       = round(sm_all$conf.int[1, "exp(coef)"], 4),
    hr_lo95  = round(sm_all$conf.int[1, "lower .95"], 4),
    hr_hi95  = round(sm_all$conf.int[1, "upper .95"], 4),
    p_value  = signif(sm_all$coefficients[1, "Pr(>|z|)"], 4),
    c_index  = round(boot_all$c_index, 4),
    c_lo95   = round(boot_all$ci_low, 4),
    c_hi95   = round(boot_all$ci_high, 4),
    stringsAsFactors = FALSE
  ),
  res_df
)

# Add interaction test as metadata
res_df$interaction_hr <- interaction_hr
res_df$interaction_p  <- interaction_p

out_path <- file.path(PATHS$results$supp, "er_stratified_analysis.csv")
write.csv(res_df, out_path, row.names = FALSE)
h <- sha256_file(out_path)
registry_append("SCANB", "er_stratified_analysis", out_path, h, "ok", SCRIPT_NAME,
                file.info(out_path)$size / 1e6)
message(sprintf("[%s] Results saved: %s", SCRIPT_NAME, out_path))

# --------------------------------------------------------------------------
# 7) Forest plot: ER-stratified
# --------------------------------------------------------------------------
plot_df <- res_df[!is.na(res_df$hr), ]
plot_df$label <- sprintf("%s (n=%d, ev=%d)", plot_df$subgroup, plot_df$n, plot_df$n_events)
plot_df$hr_label <- sprintf("%.2f (%.2f\u2013%.2f)", plot_df$hr, plot_df$hr_lo95, plot_df$hr_hi95)

# Order: Overall at bottom
plot_df$subgroup_f <- factor(plot_df$subgroup,
                              levels = rev(c("Overall", "ER+", "ER-")))

p_forest <- ggplot(plot_df, aes(x = hr, y = subgroup_f)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(xmin = hr_lo95, xmax = hr_hi95),
                width = 0.2, linewidth = 0.8, color = COL$primary) +
  geom_point(aes(shape = ifelse(subgroup == "Overall", 18L, 15L),
                 size  = ifelse(subgroup == "Overall", 5, 3.5)),
             color = COL$primary) +
  geom_text(aes(x = max(hr_hi95, na.rm = TRUE) * 1.15, label = hr_label),
            hjust = 0, size = 3.2, color = "grey30") +
  scale_shape_identity() +
  scale_size_identity() +
  scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3)) +
  coord_cartesian(xlim = c(0.3, max(plot_df$hr_hi95, na.rm = TRUE) * 2)) +
  labs(
    title    = "CorePAM HR per 1 SD by ER Status (SCAN-B)",
    subtitle = sprintf("Interaction test: HR=%.2f, p=%.3f",
                        ifelse(is.na(interaction_hr), NA, interaction_hr),
                        ifelse(is.na(interaction_p), NA, interaction_p)),
    x        = "Hazard Ratio (log scale)",
    y        = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, color = "grey40"),
    axis.text.y   = element_text(size = 11)
  )

fig_save(p_forest, "FigS_ER_Stratified_Forest", lang = "en", section = "supp",
         w = 9, h = 5)

message(sprintf("[%s] Forest plot saved", SCRIPT_NAME))
message(sprintf("[%s] ========== DONE ==========", SCRIPT_NAME))
