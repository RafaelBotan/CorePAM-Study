# =============================================================================
# SCRIPT: 39_fix_figures_F2_F3.R
# PURPOSE: Generate the two figures whose captions did not match the included
#          PNG in thesis v1.0 (post-defense audit, 2026-06):
#   F2) FigS_Cindex_Comparison  — CorePAM vs ROR-S vs OncotypeDX-RS C-index by
#       cohort, from results/supp/EXP_genefu_comparison.csv (script 26 produced
#       only the CSV, never a figure).
#   F3) FigS_Frozen_vs_Intra    — per-patient scatter of frozen vs intra-cohort
#       CorePAM score, reproducing the exact scoring of 30_frozen_zscore_
#       sensitivity.R (script 30 produced only the CSV, never a figure).
# Reads only validated intermediate outputs; recomputes nothing analytical
# except the per-sample scores (identical formula to script 30).
# Writes PT + EN, PNG + PDF, into figures/supp/{lang}/{ext}.
# =============================================================================

suppressPackageStartupMessages({
  library(arrow); library(ggplot2); library(dplyr); library(tidyr)
})

ROOT <- "Y:/Doutorado Botan/Estudos/CorePAM_Study_accepted"
SUPP <- function(lang, ext) file.path(ROOT, "figures", "supp", lang, ext)
for (l in c("pt","en")) for (e in c("png","pdf")) dir.create(SUPP(l,e), showWarnings=FALSE, recursive=TRUE)

save_both <- function(plot, name, w=9, h=5.5, dpi=300) {
  for (lang in c("pt","en")) {} # placeholder (lang handled by caller)
}

# ---------------------------------------------------------------------------
# F2 — 3-way C-index comparison
# ---------------------------------------------------------------------------
g <- read.csv(file.path(ROOT, "results/supp/EXP_genefu_comparison.csv"),
              stringsAsFactors=FALSE, check.names=FALSE)
coh_lab <- c(SCANB="SCAN-B", TCGA_BRCA="TCGA-BRCA", METABRIC="METABRIC", GSE20685="GSE20685")
f2 <- data.frame(
  cohort = factor(coh_lab[g$cohort], levels=unname(coh_lab)),
  CorePAM = g$c_corepam, `ROR-S` = g$c_ror_s, `OncotypeDX-RS` = g$c_odx_rs,
  check.names=FALSE
) |>
  pivot_longer(-cohort, names_to="Assinatura", values_to="C")
f2$Assinatura <- factor(f2$Assinatura, levels=c("CorePAM","ROR-S","OncotypeDX-RS"))

make_f2 <- function(lang) {
  ttl <- if (lang=="pt") "Discriminação (C-index de Harrell) por coorte" else "Discrimination (Harrell C-index) by cohort"
  sub <- if (lang=="pt") "CorePAM (24 genes) vs ROR-S (PAM50) vs OncotypeDX-RS (genefu)" else "CorePAM (24 genes) vs ROR-S (PAM50) vs OncotypeDX-RS (genefu)"
  ylb <- if (lang=="pt") "C-index de Harrell" else "Harrell C-index"
  ref <- if (lang=="pt") "Acaso (0,5)" else "Chance (0.5)"
  ggplot(f2, aes(cohort, C, fill=Assinatura)) +
    geom_col(position=position_dodge(0.8), width=0.75, colour="grey25", linewidth=0.2) +
    geom_text(aes(label=formatC(C, format="f", digits=3)),
              position=position_dodge(0.8), vjust=-0.4, size=2.7) +
    geom_hline(yintercept=0.5, linetype="dashed", colour="grey50") +
    annotate("text", x=0.6, y=0.505, label=ref, hjust=0, size=2.6, colour="grey45") +
    scale_fill_manual(values=c("CorePAM"="#2c7fb8","ROR-S"="#7fcdbb","OncotypeDX-RS"="#c7e9b4")) +
    coord_cartesian(ylim=c(0.5, 0.72)) +
    labs(title=ttl, subtitle=sub, x=NULL, y=ylb, fill=NULL) +
    theme_classic(base_size=11) +
    theme(legend.position="top", plot.title=element_text(face="bold", size=12))
}
for (lang in c("pt","en")) {
  p <- make_f2(lang)
  ggsave(file.path(SUPP(lang,"png"), sprintf("FigS_Cindex_Comparison_%s.png", toupper(lang))),
         p, width=9, height=5.5, dpi=300)
  ggsave(file.path(SUPP(lang,"pdf"), sprintf("FigS_Cindex_Comparison_%s.pdf", toupper(lang))),
         p, width=9, height=5.5, device=cairo_pdf)
}
cat("F2 done\n")

# ---------------------------------------------------------------------------
# F3 — per-patient frozen vs intra-cohort score scatter
#      (exact scoring formula from 30_frozen_zscore_sensitivity.R)
# ---------------------------------------------------------------------------
w_df <- read.csv(file.path(ROOT, "results/corepam/CorePAM_weights.csv"), stringsAsFactors=FALSE)
w_df <- w_df[w_df$weight != 0, ]
panel_genes   <- w_df$gene
panel_weights <- setNames(w_df$weight, w_df$gene)

ref_df <- read.csv(file.path(ROOT, "results/corepam/SCANB_reference_meanSD.csv"), stringsAsFactors=FALSE)
ref_m  <- setNames(ref_df$mean_scanb, ref_df$gene)
ref_s  <- setNames(ref_df$sd_scanb,   ref_df$gene); ref_s[ref_s==0 | is.na(ref_s)] <- 1

cohorts <- c("SCANB","TCGA_BRCA","METABRIC","GSE20685","GSE1456")
coh_lab2 <- c(SCANB="SCAN-B", TCGA_BRCA="TCGA-BRCA", METABRIC="METABRIC",
              GSE20685="GSE20685", GSE1456="GSE1456")
plat <- c(SCANB="RNA-seq", TCGA_BRCA="RNA-seq", METABRIC="Microarray",
          GSE20685="Microarray", GSE1456="Microarray")

rows <- list(); rlab <- list()
for (coh in cohorts) {
  ep <- file.path(ROOT, "01_Base_Pura_CorePAM/PROCESSED", coh, "expression_genelevel_preZ.parquet")
  if (!file.exists(ep)) next
  ex <- read_parquet(ep)
  ex <- ex[ex$gene %in% panel_genes, , drop=FALSE]       # keep only panel genes (memory-safe)
  gn <- ex$gene
  sids <- setdiff(names(ex), "gene")
  M <- t(as.matrix(ex[, sids])); colnames(M) <- gn       # samples x genes
  gp <- intersect(panel_genes, intersect(gn, ref_df$gene))
  # frozen z (SCAN-B reference)
  zf <- sweep(M[, gp, drop=FALSE], 2, ref_m[gp], "-"); zf <- sweep(zf, 2, ref_s[gp], "/")
  # intra z (cohort mean/sd)
  im <- colMeans(M[, gp, drop=FALSE], na.rm=TRUE)
  is <- apply(M[, gp, drop=FALSE], 2, sd, na.rm=TRUE); is[is==0 | is.na(is)] <- 1
  zi <- sweep(M[, gp, drop=FALSE], 2, im, "-"); zi <- sweep(zi, 2, is, "/")
  wp <- panel_weights[gp]; den <- sum(abs(wp))
  s_fro <- as.vector(zf %*% wp) / den
  s_int <- as.vector(zi %*% wp) / den
  r <- cor(s_fro, s_int, use="complete.obs")
  # standardize each for display so the y=x diagonal is interpretable (r is scale-invariant)
  rows[[coh]] <- data.frame(cohort=coh_lab2[coh], platform=plat[coh],
                            intra=as.numeric(scale(s_int)), frozen=as.numeric(scale(s_fro)))
  rlab[[coh]] <- data.frame(cohort=coh_lab2[coh], r=r)
  cat(sprintf("F3 %s: n=%d genes=%d r=%.3f\n", coh, length(s_fro), length(gp), r))
}
sc <- bind_rows(rows); sc$cohort <- factor(sc$cohort, levels=unname(coh_lab2))
labs_r <- bind_rows(rlab); labs_r$cohort <- factor(labs_r$cohort, levels=unname(coh_lab2))

make_f3 <- function(lang) {
  ttl <- if (lang=="pt") "Concordância: escore intra-coorte vs z-score congelado (parâmetros SCAN-B)" else "Agreement: intra-cohort vs frozen z-score (SCAN-B parameters)"
  sub <- if (lang=="pt") "Cada ponto = uma paciente (escores padronizados); linha = concordância perfeita (y=x)" else "Each point = one patient (standardized scores); line = perfect agreement (y=x)"
  xlb <- if (lang=="pt") "Escore intra-coorte (padronizado)" else "Intra-cohort score (standardized)"
  ylb <- if (lang=="pt") "Escore congelado (padronizado)" else "Frozen score (standardized)"
  labs_r$txt <- sprintf("r = %.3f", labs_r$r)
  ggplot(sc, aes(intra, frozen)) +
    geom_abline(slope=1, intercept=0, colour="grey55", linewidth=0.5) +
    geom_point(aes(colour=platform), alpha=0.35, size=0.7) +
    geom_text(data=labs_r, aes(x=-2.6, y=2.6, label=txt), hjust=0, size=3, inherit.aes=FALSE) +
    facet_wrap(~cohort, nrow=2) +
    scale_colour_manual(values=c("RNA-seq"="#2c7fb8","Microarray"="#d95f0e")) +
    coord_cartesian(xlim=c(-3,3), ylim=c(-3,3)) +
    labs(title=ttl, subtitle=sub, x=xlb, y=ylb, colour=NULL) +
    theme_classic(base_size=11) +
    theme(legend.position="top", plot.title=element_text(face="bold", size=11.5),
          strip.background=element_rect(fill="grey92", colour=NA))
}
for (lang in c("pt","en")) {
  p <- make_f3(lang)
  ggsave(file.path(SUPP(lang,"png"), sprintf("FigS_Frozen_vs_Intra_%s.png", toupper(lang))),
         p, width=9, height=6, dpi=300)
  ggsave(file.path(SUPP(lang,"pdf"), sprintf("FigS_Frozen_vs_Intra_%s.pdf", toupper(lang))),
         p, width=9, height=6, device=cairo_pdf)
}
cat("F3 done\n")
cat("ALL FIGURES WRITTEN\n")
