options(warn = 0)
suppressPackageStartupMessages({
  library(ggplot2)
  library(jsonlite)
})

source("Y:/Phd-Genomic-claude/scripts/00_colors.R")

# ── Data ──────────────────────────────────────────────────────────────────────
pareto <- read.csv("results/corepam/pareto_df_cindex_oof.csv")
pareto <- pareto[!is.na(pareto$c_adj), ]

j        <- fromJSON("results/corepam/selected_CorePAM_summary.json")
df_sel   <- j$selected$df
df_cmax  <- pareto$df[which.max(pareto$c_adj)]
c_sel    <- j$selected$c_adj
c_max    <- j$selected$c_max
gap_val  <- j$selected$gap
delta_c  <- j$delta_c
thresh   <- c_max - delta_c

cat(sprintf("df_sel=%d | df_cmax=%d | c_sel=%.4f | c_max=%.4f | gap=%.4f | threshold=%.4f\n",
            df_sel, df_cmax, c_sel, c_max, gap_val, thresh))

pt_sel   <- pareto[pareto$df == df_sel, ]
pt_cmax  <- pareto[pareto$df == df_cmax, ]
x_breaks <- sort(unique(c(0, 5, 10, 15, 20, df_sel, df_cmax, max(pareto$df))))

# ── Plot builder ──────────────────────────────────────────────────────────────
make_plot <- function(lang = c("EN", "PT")) {
  lang <- match.arg(lang)

  if (lang == "EN") {
    x_lab  <- "Number of non-zero genes"
    y_lab  <- "OOF Harrell C-index"
    lbl_sel  <- sprintf("Core-PAM\n(n = %d genes)", df_sel)
    lbl_cmax <- sprintf("C-index maximum\n(n = %d genes)", df_cmax)
    lbl_thr  <- expression(C[max] - Delta*C)
  } else {
    x_lab  <- "N\u00famero de genes selecionados"
    y_lab  <- "C-index de Harrell (OOF)"
    lbl_sel  <- sprintf("Core-PAM\n(n = %d genes)", df_sel)
    lbl_cmax <- sprintf("C-index m\u00e1ximo\n(n = %d genes)", df_cmax)
    lbl_thr  <- expression(C[max] - Delta*C)
  }

  # Label positions — placed in clear open areas, connected by arrows to points
  # Selected (df=24): upper-left clear zone (above plateau, left of mid-figure)
  sel_lbl_x <- 16;    sel_lbl_y <- 0.685
  # C_max (df=36): upper-right area, above dotted line
  cmax_lbl_x <- df_cmax + 2;  cmax_lbl_y <- c_max + 0.0048
  # Threshold label: far left margin, just above dashed line
  thr_lbl_x <- 2;     thr_lbl_y <- thresh + 0.0022

  ggplot(pareto, aes(x = df, y = c_adj)) +

    # Non-inferiority zone (subtle fill below threshold)
    annotate("rect",
             xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = thresh,
             fill = COL$fill_warn, alpha = 0.18) +

    # Threshold dashed line
    geom_hline(yintercept = thresh,
               linetype = "dashed", colour = COL$threshold, linewidth = 0.55) +

    # C_max dotted reference
    geom_hline(yintercept = c_max,
               linetype = "dotted", colour = COL$grey_mid, linewidth = 0.45) +

    # Pareto path
    geom_line(colour = COL$primary, linewidth = 1.0) +
    geom_point(colour = COL$primary, fill = COL$primary, size = 1.6, shape = 21) +

    # Highlighted points
    geom_point(data = pt_sel,  aes(x = df, y = c_adj),
               colour = COL$selected, fill = COL$selected, size = 4.5, shape = 24) +
    geom_point(data = pt_cmax, aes(x = df, y = c_adj),
               colour = COL$reference, fill = COL$reference, size = 4, shape = 21) +

    # Arrow from label to selected point
    annotate("segment",
             x = sel_lbl_x + 1.2, xend = pt_sel$df - 0.4,
             y = sel_lbl_y - 0.0008, yend = pt_sel$c_adj + 0.0006,
             colour = COL$selected, linewidth = 0.4,
             arrow = arrow(length = unit(0.12, "cm"), type = "closed")) +
    annotate("text",
             x = sel_lbl_x, y = sel_lbl_y,
             label = lbl_sel,
             colour = COL$selected, size = 3.2, hjust = 1, lineheight = 0.9,
             fontface = "bold") +

    # Arrow from label to C_max point
    annotate("segment",
             x = cmax_lbl_x - 0.5, xend = pt_cmax$df + 0.3,
             y = cmax_lbl_y - 0.0010, yend = pt_cmax$c_adj + 0.0006,
             colour = COL$reference, linewidth = 0.4,
             arrow = arrow(length = unit(0.12, "cm"), type = "closed")) +
    annotate("text",
             x = cmax_lbl_x, y = cmax_lbl_y,
             label = lbl_cmax,
             colour = COL$reference, size = 3.2, hjust = 0, lineheight = 0.9) +

    # Threshold label: left margin, above dashed line (far from both point labels)
    annotate("text",
             x = thr_lbl_x, y = thr_lbl_y,
             label = lbl_thr,
             colour = COL$threshold, size = 2.9, hjust = 0, parse = TRUE) +

    # Axes
    scale_x_continuous(breaks = x_breaks, minor_breaks = NULL,
                       expand = expansion(mult = c(0.04, 0.06))) +
    scale_y_continuous(
      limits = c(0.52, 0.695),
      breaks = seq(0.53, 0.69, by = 0.02),
      labels = function(x) sprintf("%.2f", x),
      expand = expansion(mult = c(0.02, 0.03))
    ) +
    labs(x = x_lab, y = y_lab) +

    theme_classic(base_size = 12) +
    theme(
      axis.line        = element_line(colour = "black", linewidth = 0.5),
      axis.ticks       = element_line(colour = "black", linewidth = 0.4),
      axis.text        = element_text(colour = "black", size = 11),
      axis.title       = element_text(colour = "black", size = 12),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.margin      = margin(12, 18, 8, 10)
    )
}

# ── Save (do NOT remove previous versions — user validation workflow) ─────────
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

dir.create("results/figures/supp", showWarnings = FALSE, recursive = TRUE)

for (lang in c("EN", "PT")) {
  p        <- make_plot(lang)
  pdf_path <- sprintf("results/figures/supp/FigS2_Pareto_CorePAM_%s.pdf", lang)
  png_path <- sprintf("results/figures/supp/FigS2_Pareto_CorePAM_%s.png", lang)
  # PDF: cairo_pdf = true vector, fully editable in Illustrator/Inkscape
  ggsave(pdf_path, plot = p, width = 8.5, height = 5.2,
         device = cairo_pdf, family = "sans")
  # PNG: 600 DPI for high-resolution print/screen
  ggsave(png_path, plot = p, width = 8.5, height = 5.2, dpi = 600)
  cat(sprintf("[%s] %s | %s\n", lang, pdf_path, png_path))
}

cat("Done.\n")
