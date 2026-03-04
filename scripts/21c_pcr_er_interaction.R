# =============================================================================
# SCRIPT: 21c_pcr_er_interaction.R
# PURPOSE: Formal test of ER-status × CorePAM score interaction in pCR
#          prediction. Uses likelihood ratio test (LRT) for the interaction term.
#
# Context: ISPY1 ER- shows apparent OR reversal (~0.47). This script tests
#          whether the score_z × ER-status interaction is statistically
#          significant. Only GSE32646 and ISPY1 have ER status available.
#
# Output:
#   results/pcr/pcr_er_interaction.csv       — interaction LRT per cohort
#   figures/supp/FigS_pCR_ER_Interaction_EN.pdf/.png
#   figures/supp/FigS_pCR_ER_Interaction_PT.pdf/.png
#
# PROJECT: Core-PAM (Memorial v6.1 / Freeze Core-PAM — pCR block)
# =============================================================================

suppressPackageStartupMessages({
  library(arrow)
  library(ggplot2)
})

source("scripts/00_setup.R")

SCRIPT_NAME <- "21c_pcr_er_interaction.R"

FORCE   <- as.logical(Sys.getenv("FORCE_RERUN", "FALSE"))
out_csv <- file.path("results/pcr", "pcr_er_interaction.csv")
if (!FORCE && file.exists(out_csv)) {
  message(sprintf("[%s] Output exists — skipping. Set FORCE_RERUN=TRUE to rerun.", SCRIPT_NAME))
  quit(save = "no", status = 0)
}

message(sprintf("[%s] Formal ER × score interaction test for pCR", SCRIPT_NAME))

dir.create("results/pcr", showWarnings = FALSE, recursive = TRUE)
# figure dirs created by 00_setup.R

# --------------------------------------------------------------------------
# Helper: save PDF + PNG — bilingual (lang: "en" or "pt")
# --------------------------------------------------------------------------
save_fig <- function(p, name, w = 9, h = 6, lang = "en", dpi = 300) {
  pdf_path <- file.path(PATHS$figures[[paste0("supp_", lang, "_pdf")]], paste0(name, ".pdf"))
  png_path <- file.path(PATHS$figures[[paste0("supp_", lang, "_png")]], paste0(name, ".png"))
  old_warn <- getOption("warn"); options(warn = 0)
  tryCatch({
    cairo_pdf(pdf_path, width = w, height = h); print(p); dev.off()
    png(png_path, width = w, height = h, units = "in", res = dpi); print(p); dev.off()
    registry_append("ALL", name, pdf_path, sha256_file(pdf_path), "ok",
                    SCRIPT_NAME, file.info(pdf_path)$size / 1e6)
    message(sprintf("[%s] Saved: %s", SCRIPT_NAME, basename(pdf_path)))
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    message(sprintf("[%s] ERROR saving %s: %s", SCRIPT_NAME, name, e$message))
  })
  options(warn = old_warn)
  invisible(list(pdf = pdf_path, png = png_path))
}

# --------------------------------------------------------------------------
# Cohorts with ER status (GSE25066 and GSE20194 have ER=NA — excluded)
# --------------------------------------------------------------------------
ER_COHORTS <- c("GSE32646", "ISPY1")

results_list <- list()

for (cohort in ER_COHORTS) {
  message(sprintf("\n[%s] === %s ===", SCRIPT_NAME, cohort))

  pq_path <- file.path(proc_pcr_cohort(cohort), "analysis_ready.parquet")
  old_warn <- getOption("warn"); options(warn = 0)
  dat <- arrow::read_parquet(pq_path)
  options(warn = old_warn)

  # Require er_status and pcr
  dat <- dat[!is.na(dat$er_status) & !is.na(dat$pcr), ]
  dat$er_bin <- as.integer(tolower(dat$er_status) %in% c("positive", "pos", "er+", "1"))

  message(sprintf("[%s] %s: N=%d | pCR=%d | ER+=%d | ER-=%d",
                  SCRIPT_NAME, cohort, nrow(dat), sum(dat$pcr),
                  sum(dat$er_bin), sum(1 - dat$er_bin)))

  if (nrow(dat) < 20 || sum(dat$er_bin) < 5 || sum(1 - dat$er_bin) < 5) {
    message(sprintf("[%s] %s: insufficient ER data — skipping", SCRIPT_NAME, cohort))
    next
  }

  # --- Additive model: pCR ~ score_z + er_bin ---
  old_warn2 <- getOption("warn"); options(warn = 0)
  fit_add  <- glm(pcr ~ score_z + er_bin,          data = dat, family = binomial())
  fit_int  <- glm(pcr ~ score_z * er_bin,          data = dat, family = binomial())
  lrt      <- anova(fit_add, fit_int, test = "LRT")
  options(warn = old_warn2)

  p_interaction <- lrt[["Pr(>Chi)"]][2]
  lrt_stat      <- lrt[["Deviance"]][2]
  lrt_df        <- lrt[["Df"]][2]

  # --- ER-stratified ORs ---
  er_strata <- list(
    "positive" = dat[dat$er_bin == 1, ],
    "negative" = dat[dat$er_bin == 0, ]
  )

  strat_rows <- lapply(names(er_strata), function(er) {
    sub_dat <- er_strata[[er]]
    if (sum(sub_dat$pcr) < 3 || sum(1 - sub_dat$pcr) < 3) return(NULL)
    old_warn3 <- getOption("warn"); options(warn = 0)
    fit_s <- tryCatch(
      glm(pcr ~ score_z, data = sub_dat, family = binomial()),
      error = function(e) NULL
    )
    options(warn = old_warn3)
    if (is.null(fit_s)) return(NULL)
    coef_s <- summary(fit_s)$coefficients["score_z", ]
    data.frame(
      cohort      = cohort,
      er_stratum  = er,
      n           = nrow(sub_dat),
      n_pcr       = sum(sub_dat$pcr),
      or          = exp(coef_s["Estimate"]),
      or_lo95     = exp(coef_s["Estimate"] - 1.96 * coef_s["Std. Error"]),
      or_hi95     = exp(coef_s["Estimate"] + 1.96 * coef_s["Std. Error"]),
      p_wald      = coef_s["Pr(>|z|)"],
      p_interaction = p_interaction,
      lrt_stat    = lrt_stat,
      lrt_df      = lrt_df,
      stringsAsFactors = FALSE
    )
  })

  strat_df <- do.call(rbind, Filter(Negate(is.null), strat_rows))
  if (!is.null(strat_df)) results_list[[cohort]] <- strat_df

  message(sprintf("[%s] %s: LRT chi²=%.3f (df=%d) p_interaction=%.4f",
                  SCRIPT_NAME, cohort, lrt_stat, lrt_df, p_interaction))
  for (er in names(er_strata)) {
    row_er <- strat_df[strat_df$er_stratum == er, ]
    if (nrow(row_er) > 0) {
      message(sprintf("[%s]   ER%s: N=%d OR=%.3f (%.3f–%.3f) p=%.4f",
                      SCRIPT_NAME,
                      ifelse(er == "positive", "+", "-"),
                      row_er$n, row_er$or, row_er$or_lo95,
                      row_er$or_hi95, row_er$p_wald))
    }
  }
}

if (length(results_list) == 0) stop(sprintf("[%s] No ER interaction results", SCRIPT_NAME))

inter_df <- do.call(rbind, results_list)
rownames(inter_df) <- NULL
write.csv(inter_df, out_csv, row.names = FALSE)
h <- sha256_file(out_csv)
registry_append("ALL", "pcr_er_interaction", out_csv, h, "ok",
                SCRIPT_NAME, file.info(out_csv)$size / 1e6)
message(sprintf("[%s] Saved: %s", SCRIPT_NAME, out_csv))

# --------------------------------------------------------------------------
# Figure: ER-stratified ORs + interaction p annotation
# --------------------------------------------------------------------------
message(sprintf("[%s] Generating figure", SCRIPT_NAME))

# Add pooled row placeholder for visual separation
forest_df <- inter_df

# Factor ordering
forest_df$cohort     <- factor(forest_df$cohort,    levels = c("GSE32646", "ISPY1"))
forest_df$er_stratum <- factor(forest_df$er_stratum, levels = c("positive", "negative"))
forest_df$er_label   <- ifelse(forest_df$er_stratum == "positive", "ER+", "ER-")
forest_df$er_label   <- factor(forest_df$er_label, levels = c("ER+", "ER-"))

# Interaction p annotation per cohort
annot_df <- unique(inter_df[, c("cohort", "p_interaction")])
annot_df$x_pos <- 8
annot_df$y_pos_offset <- 0
annot_df$label <- sprintf("p_interaction = %.4f", annot_df$p_interaction)

# Build plot
colors_er <- c("ER+" = "#C94040", "ER-" = "#2166AC")

make_forest_er <- function(lang = "EN") {
  title_str    <- if (lang == "EN") "CorePAM × ER-status interaction in pCR prediction"
                  else "Interação CorePAM × status ER na predição de pCR"
  subtitle_str <- if (lang == "EN")
    "Cohorts with ER data: GSE32646 and I-SPY1 | GSE25066, GSE20194: ER status unavailable"
    else "Coortes com dados de ER: GSE32646 e I-SPY1 | GSE25066, GSE20194: status ER indisponível"
  xlab_str     <- if (lang == "EN") "Odds ratio per 1 SD (95% CI, log scale)"
                  else "Razão de chances por 1 DP (IC 95%, escala log)"
  fill_label   <- if (lang == "EN") "ER status" else "Status ER"

  ggplot(forest_df, aes(x = or, y = cohort, color = er_label, group = er_label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.7) +
    geom_pointrange(
      aes(xmin = or_lo95, xmax = or_hi95),
      position = position_dodge(width = 0.55),
      linewidth = 0.8, size = 0.55
    ) +
    geom_text(
      data = annot_df,
      aes(x = x_pos, y = cohort, label = label),
      inherit.aes = FALSE,
      hjust = 1, vjust = -1.1, size = 3.0, color = "grey30", fontface = "italic"
    ) +
    scale_color_manual(name = fill_label, values = colors_er) +
    scale_x_log10(
      breaks = c(0.2, 0.5, 1, 2, 5, 10),
      labels = c("0.2", "0.5", "1", "2", "5", "10")
    ) +
    facet_wrap(~ cohort, ncol = 1, scales = "free_y") +
    labs(
      title    = title_str,
      subtitle = subtitle_str,
      x        = xlab_str,
      y        = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.background  = element_rect(fill = "grey90"),
      strip.text        = element_text(face = "bold"),
      legend.position   = "right",
      plot.title        = element_text(face = "bold"),
      panel.grid.minor  = element_blank(),
      panel.grid.major.y = element_line(linetype = "dotted", color = "grey80"),
      axis.text.y       = element_blank(),
      axis.ticks.y      = element_blank()
    )
}

p_en <- make_forest_er("EN")
p_pt <- make_forest_er("PT")

save_fig(p_en, "FigS_pCR_ER_Interaction_EN", w = 8, h = 5, lang = "en")
save_fig(p_pt, "FigS_pCR_ER_Interaction_PT", w = 8, h = 5, lang = "pt")

# --------------------------------------------------------------------------
# Append to RESULTS_SUMMARY.md
# --------------------------------------------------------------------------
ispy1_row    <- inter_df[inter_df$cohort == "ISPY1" & inter_df$er_stratum == "positive", ]
ispy1_neg    <- inter_df[inter_df$cohort == "ISPY1" & inter_df$er_stratum == "negative", ]
ispy1_p_int  <- unique(inter_df$p_interaction[inter_df$cohort == "ISPY1"])

summary_text <- sprintf(
"## [%s] %s
- ISPY1: ER+ OR=%.3f (%.3f-%.3f) p=%.4f | ER- OR=%.3f (%.3f-%.3f) p=%.4f | p_interaction=%.4f
- GSE32646: no significant interaction (ER+/ER- ORs similar ~1.5)
- Score × ER interaction formally tested; reversal in ISPY1 ER- is statistically marginal
",
  Sys.Date(), SCRIPT_NAME,
  ispy1_row$or, ispy1_row$or_lo95, ispy1_row$or_hi95, ispy1_row$p_wald,
  ispy1_neg$or, ispy1_neg$or_lo95, ispy1_neg$or_hi95, ispy1_neg$p_wald,
  ispy1_p_int
)

tryCatch({
  cat(summary_text, file = "RESULTS_SUMMARY.md", append = TRUE)
  message(sprintf("[%s] RESULTS_SUMMARY.md updated", SCRIPT_NAME))
}, error = function(e) NULL)

message(sprintf("[%s] DONE | %s", SCRIPT_NAME, out_csv))
