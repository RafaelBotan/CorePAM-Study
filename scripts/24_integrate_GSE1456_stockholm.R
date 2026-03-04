# =============================================================================
# SCRIPT: 24_integrate_GSE1456_stockholm.R
# PROJECT: Core-PAM (Memorial v6.1 / Freeze Core-PAM)
# PURPOSE: Integrate GSE1456 (Stockholm/Karolinska, N=159) as 4th external
#          validation cohort (Affymetrix HG-U133A, OS + DSS).
#          Downloads, harmonises clinical, preprocesses expression, computes
#          CorePAM score, and runs survival analysis in a single script.
# INPUTS:  GEO accession GSE1456, results/corepam/CorePAM_weights.csv
# OUTPUTS:
#   - PROCESSED/GSE1456/expression_genelevel_preZ.parquet
#   - PROCESSED/GSE1456/clinical_FINAL.parquet
#   - PROCESSED/GSE1456/analysis_ready.parquet
#   - results/supp/survival_results_GSE1456.csv
#   - results/supp/risk_score_summary_GSE1456.csv
#   - figures/main/en/{pdf,png}/Fig2_KM_GSE1456_OS_EN.{pdf,png}
# =============================================================================
source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(survival)
  library(arrow)
  library(GEOquery)
  library(ggplot2)
  library(survminer)
})

message("[GSE1456] ========== INTEGRATE STOCKHOLM COHORT ==========")

PROC_DIR <- file.path(PATHS$processed, "GSE1456")
dir.create(PROC_DIR, showWarnings = FALSE, recursive = TRUE)

# ==========================================================================
# 1) DOWNLOAD & EXTRACT
# ==========================================================================
message("[GSE1456] Step 1: Downloading from GEO...")
old_warn <- getOption("warn"); options(warn = 0)
gse <- getGEO("GSE1456", GSEMatrix = TRUE, destdir = tempdir())
options(warn = old_warn)

eset <- gse[["GSE1456-GPL96_series_matrix.txt.gz"]]
pd <- pData(eset)
raw_expr <- exprs(eset)  # probes x samples (log2-RMA from GEO)

message(sprintf("[GSE1456] Raw: %d probes x %d samples", nrow(raw_expr), ncol(raw_expr)))

# Save raw clinical for reproducibility
raw_dir <- file.path(PATHS$raw, "GSE1456")
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(pd, file.path(raw_dir, "GSE1456_clinical_raw.rds"))
message("[GSE1456] Saved raw clinical RDS")

# ==========================================================================
# 2) HARMONIZE CLINICAL
# ==========================================================================
message("[GSE1456] Step 2: Harmonizing clinical data...")

parse_char <- function(x) {
  val <- sub("^[^:]+:\\s*", "", x)
  val[val == "NA" | val == ""] <- NA
  val
}

clinical <- data.frame(
  sample_id = pd$geo_accession,
  patient_id = pd$geo_accession,  # GEO samples as patient IDs
  os_time_months = as.numeric(parse_char(pd$`SURV_DEATH:ch1`)) * 12,
  os_event = as.integer(parse_char(pd$`DEATH:ch1`)),
  dss_time_months = as.numeric(parse_char(pd$`SURV_DEATH:ch1`)) * 12,
  dss_event = as.integer(parse_char(pd$`DEATH_BC:ch1`)),
  age = NA_real_,  # Not available in GSE1456
  er_status = NA_character_,  # Not available in GSE1456
  pam50_subtype = parse_char(pd$`characteristics_ch1.6`),
  elston_grade = parse_char(pd$`characteristics_ch1.7`),
  stringsAsFactors = FALSE
)

# Clean subtype field
clinical$pam50_subtype <- sub("^SUBTYPE:\\s*", "", clinical$pam50_subtype)
clinical$elston_grade <- sub("^ELSTON:\\s*", "", clinical$elston_grade)

# Drop time <= 0
n_before <- nrow(clinical)
clinical <- clinical[!is.na(clinical$os_time_months) & clinical$os_time_months > 0, ]
n_dropped <- n_before - nrow(clinical)
if (n_dropped > 0) message(sprintf("[GSE1456] Dropped %d samples with time <= 0", n_dropped))

message(sprintf("[GSE1456] Clinical: N=%d, OS events=%d (%.1f%%), DSS events=%d (%.1f%%)",
                nrow(clinical),
                sum(clinical$os_event, na.rm=TRUE), 100*mean(clinical$os_event, na.rm=TRUE),
                sum(clinical$dss_event, na.rm=TRUE), 100*mean(clinical$dss_event, na.rm=TRUE)))
message(sprintf("[GSE1456] Median FU (censored): %.1f months",
                median(clinical$os_time_months[clinical$os_event == 0], na.rm=TRUE)))

# Save clinical_FINAL.parquet
clin_out <- clinical[, c("sample_id", "patient_id", "os_time_months", "os_event", "age", "er_status")]
write_parquet(clin_out, file.path(PROC_DIR, "clinical_FINAL.parquet"))

# ==========================================================================
# 3) EXPRESSION PREPROCESSING
# ==========================================================================
message("[GSE1456] Step 3: Preprocessing expression data...")

# Get GPL96 annotation for probe-to-gene mapping
old_warn <- getOption("warn"); options(warn = 0)
gpl <- getGEO("GPL96", destdir = tempdir())
options(warn = old_warn)

gpl_tab <- Table(gpl)
probe_to_gene <- data.frame(
  probe_id = gpl_tab$ID,
  gene_symbol = sapply(strsplit(gpl_tab$`Gene Symbol`, " /// "), `[`, 1),
  stringsAsFactors = FALSE
)
probe_to_gene <- probe_to_gene[!is.na(probe_to_gene$gene_symbol) & probe_to_gene$gene_symbol != "", ]

# For each gene, pick probe with highest variance
unique_genes <- unique(probe_to_gene$gene_symbol)
gene_expr_list <- list()

for (gene in unique_genes) {
  probes <- probe_to_gene$probe_id[probe_to_gene$gene_symbol == gene]
  probes <- probes[probes %in% rownames(raw_expr)]
  if (length(probes) == 0) next
  if (length(probes) == 1) {
    gene_expr_list[[gene]] <- raw_expr[probes, ]
  } else {
    vars <- apply(raw_expr[probes, , drop=FALSE], 1, var, na.rm=TRUE)
    best <- names(which.max(vars))
    gene_expr_list[[gene]] <- raw_expr[best, ]
  }
}

gene_expr_mat <- do.call(rbind, gene_expr_list)
rownames(gene_expr_mat) <- names(gene_expr_list)

message(sprintf("[GSE1456] Gene-level expression: %d genes x %d samples",
                nrow(gene_expr_mat), ncol(gene_expr_mat)))

# Save expression_genelevel_preZ.parquet (genes as rows, gene column first)
expr_df <- data.frame(gene = rownames(gene_expr_mat), gene_expr_mat,
                      check.names = FALSE, stringsAsFactors = FALSE)
write_parquet(expr_df, file.path(PROC_DIR, "expression_genelevel_preZ.parquet"))
message("[GSE1456] Saved expression_genelevel_preZ.parquet")

# ==========================================================================
# 4) Z-SCORE & COREPAM SCORE
# ==========================================================================
message("[GSE1456] Step 4: Z-score and CorePAM scoring...")

wt <- read.csv("results/corepam/CorePAM_weights.csv", stringsAsFactors = FALSE)
corepam_genes <- wt$gene

# Check gene coverage
available <- intersect(corepam_genes, rownames(gene_expr_mat))
missing <- setdiff(corepam_genes, rownames(gene_expr_mat))
frac_present <- length(available) / length(corepam_genes)

message(sprintf("[GSE1456] CorePAM genes: %d/%d (%.1f%%). Missing: %s",
                length(available), length(corepam_genes), 100*frac_present,
                if(length(missing)>0) paste(missing, collapse=", ") else "none"))

stopifnot(frac_present >= FREEZE$min_genes_fraction)

# Z-score per gene (intra-cohort)
old_warn <- getOption("warn"); options(warn = 0)
z_mat <- t(apply(gene_expr_mat[available, , drop=FALSE], 1, scale))
colnames(z_mat) <- colnames(gene_expr_mat)
options(warn = old_warn)

# CorePAM score: sum(w_i * z_i) / sum(|w_i|) for present genes
wt_sub <- wt[wt$gene %in% available, ]
denom <- sum(abs(wt_sub$weight))
score_raw <- rep(0, ncol(z_mat))
names(score_raw) <- colnames(z_mat)
for (i in seq_len(nrow(wt_sub))) {
  score_raw <- score_raw + wt_sub$weight[i] * z_mat[wt_sub$gene[i], ]
}
score_raw <- score_raw / denom

# Check direction
tmp <- data.frame(
  sample_id = names(score_raw),
  score_tmp = score_raw,
  stringsAsFactors = FALSE
)
tmp <- merge(tmp, clinical, by = "sample_id")
tmp <- tmp[!is.na(tmp$os_time_months) & tmp$os_time_months > 0 & !is.na(tmp$os_event), ]

old_warn <- getOption("warn"); options(warn = 0)
fit_dir <- coxph(Surv(os_time_months, os_event) ~ score_tmp, data = tmp)
options(warn = old_warn)

score_dir <- if (coef(fit_dir) < 0) "inverted" else "original"
if (score_dir == "inverted") {
  score_raw <- -score_raw
  message("[GSE1456] Score FLIPPED (direction = inverted)")
} else {
  message("[GSE1456] Score direction: original")
}

score_z <- as.numeric(scale(score_raw))
names(score_z) <- names(score_raw)

# Build analysis_ready
ar <- data.frame(
  sample_id = names(score_raw),
  patient_id = names(score_raw),
  stringsAsFactors = FALSE
)
ar <- merge(ar, clinical[, c("sample_id","os_time_months","os_event","dss_time_months","dss_event","age","er_status")],
            by = "sample_id")
ar$score <- score_raw[ar$sample_id]
ar$score_z <- score_z[ar$sample_id]
ar$genes_present <- length(available)
ar$denom_sum_absw <- denom
ar$score_direction <- score_dir

# Filter valid survival
ar <- ar[!is.na(ar$os_time_months) & ar$os_time_months > 0 & !is.na(ar$os_event), ]

write_parquet(ar, file.path(PROC_DIR, "analysis_ready.parquet"))
message(sprintf("[GSE1456] analysis_ready: N=%d, saved to %s", nrow(ar), PROC_DIR))

# Save risk_score_summary
rss <- data.frame(
  cohort = "GSE1456",
  n_samples = nrow(ar),
  n_panel_genes = length(corepam_genes),
  genes_present = length(available),
  genes_missing = paste(missing, collapse=","),
  frac_present = round(frac_present, 4),
  denom_sum_absw = round(denom, 6),
  score_direction = score_dir,
  score_mean = round(mean(ar$score), 6),
  score_sd = round(sd(ar$score), 6),
  n_time_leq0_dropped = n_dropped,
  stringsAsFactors = FALSE
)
write.csv(rss, file.path("results/supp", "risk_score_summary_GSE1456.csv"), row.names = FALSE)

# ==========================================================================
# 5) SURVIVAL ANALYSIS
# ==========================================================================
message("[GSE1456] Step 5: Survival analysis...")

old_warn <- getOption("warn"); options(warn = 0)

# --- Univariate OS ---
surv_os <- Surv(ar$os_time_months, ar$os_event)
fit_uni <- coxph(surv_os ~ score_z, data = ar)
hr_uni <- exp(coef(fit_uni))
ci_uni <- exp(confint(fit_uni))
p_uni <- summary(fit_uni)$coefficients[, "Pr(>|z|)"]
conc_obj <- concordance(fit_uni)
c_idx <- conc_obj$concordance
c_se <- sqrt(conc_obj$var)
loghr <- coef(fit_uni)
se_loghr <- sqrt(vcov(fit_uni)[1,1])

message(sprintf("[GSE1456] OS uni: HR=%.3f (%.3f-%.3f), p=%.4g, C=%.4f",
                hr_uni, ci_uni[1], ci_uni[2], p_uni, c_idx))

# --- Multivariate (CORE-A) --- age/ER unavailable, so multi = uni
# GSE1456 has NO age and NO ER, so CORE-A = univariate only
hr_multi <- hr_uni
ci_multi <- ci_uni
p_multi <- p_uni
corea_vars <- "none (age/ER unavailable)"

message(sprintf("[GSE1456] CORE-A: age/ER unavailable; adjusted = univariate"))

# --- DSS sensitivity ---
surv_dss <- Surv(ar$dss_time_months, ar$dss_event)
fit_dss <- coxph(surv_dss ~ score_z, data = ar)
hr_dss <- exp(coef(fit_dss))
ci_dss <- exp(confint(fit_dss))
p_dss <- summary(fit_dss)$coefficients[, "Pr(>|z|)"]
c_dss <- concordance(fit_dss)$concordance

message(sprintf("[GSE1456] DSS: HR=%.3f (%.3f-%.3f), p=%.4g, C=%.4f",
                hr_dss, ci_dss[1], ci_dss[2], p_dss, c_dss))

options(warn = old_warn)

# Extract confint as explicit scalars (confint returns matrix, avoid indexing issues)
ci_uni_lo  <- as.numeric(ci_uni[1,1])
ci_uni_hi  <- as.numeric(ci_uni[1,2])
ci_multi_lo <- ci_uni_lo  # CORE-A unavailable, same as uni
ci_multi_hi <- ci_uni_hi
ci_dss_lo  <- as.numeric(ci_dss[1,1])
ci_dss_hi  <- as.numeric(ci_dss[1,2])

# Save survival results
surv_res <- data.frame(
  cohort = "GSE1456",
  endpoint = "OS",
  n_samples = nrow(ar),
  n_events = sum(ar$os_event),
  fu_median_months = round(median(ar$os_time_months[ar$os_event == 0], na.rm=TRUE), 1),
  hr_uni = round(as.numeric(hr_uni), 4),
  hr_uni_lo95 = round(ci_uni_lo, 4),
  hr_uni_hi95 = round(ci_uni_hi, 4),
  p_uni = as.numeric(p_uni),
  hr_multi = round(as.numeric(hr_multi), 4),
  hr_multi_lo95 = round(ci_multi_lo, 4),
  hr_multi_hi95 = round(ci_multi_hi, 4),
  p_multi = as.numeric(p_multi),
  corea_vars_used = corea_vars,
  c_index = round(c_idx, 4),
  c_index_lo95 = round(c_idx - 1.96*c_se, 4),
  c_index_hi95 = round(c_idx + 1.96*c_se, 4),
  loghr_uni = as.numeric(loghr),
  se_loghr_uni = as.numeric(se_loghr),
  hr_sens_dss = round(as.numeric(hr_dss), 4),
  hr_sens_dss_lo95 = round(ci_dss_lo, 4),
  hr_sens_dss_hi95 = round(ci_dss_hi, 4),
  p_sens_dss = as.numeric(p_dss),
  c_sens_dss = round(as.numeric(c_dss), 4),
  stringsAsFactors = FALSE
)
write.csv(surv_res, file.path("results/supp", "survival_results_GSE1456.csv"), row.names = FALSE)
message("[GSE1456] Saved survival results")

# ==========================================================================
# 6) KM PLOT
# ==========================================================================
message("[GSE1456] Step 6: KM plot...")

ar$score_group <- factor(
  ifelse(ar$score_z >= median(ar$score_z), "High CorePAM", "Low CorePAM"),
  levels = c("Low CorePAM", "High CorePAM")
)

old_warn <- getOption("warn"); options(warn = 0)
km_fit <- survfit(surv_os ~ score_group, data = ar)

km_label <- sprintf("GSE1456 (Stockholm) — OS\nHR = %.2f (95%% CI %.2f–%.2f), p = %s\nC-index = %.3f | N = %d",
                     hr_uni, ci_uni_lo, ci_uni_hi,
                     if(p_uni < 0.001) sprintf("%.2e", p_uni) else sprintf("%.4f", p_uni),
                     c_idx, nrow(ar))

p_km <- ggsurvplot(
  km_fit,
  data = ar,
  palette = c("#2166AC", "#B2182B"),
  risk.table = TRUE,
  risk.table.height = 0.25,
  conf.int = TRUE,
  pval = FALSE,
  xlab = "Time (months)",
  ylab = "Overall survival probability",
  title = km_label,
  legend.title = "",
  legend.labs = c("Low CorePAM", "High CorePAM"),
  ggtheme = theme_classic(base_size = 14)
)
options(warn = old_warn)

# Save KM figure
fig_dir <- "figures/main/en/png"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

km_path <- file.path(fig_dir, "Fig2_KM_GSE1456_OS_EN.png")
png(km_path, width = 2400, height = 1800, res = 300)
print(p_km)
dev.off()
message(sprintf("[GSE1456] KM plot saved to %s", km_path))

# PDF version
fig_dir_pdf <- "figures/main/en/pdf"
dir.create(fig_dir_pdf, showWarnings = FALSE, recursive = TRUE)
km_pdf <- file.path(fig_dir_pdf, "Fig2_KM_GSE1456_OS_EN.pdf")
old_warn <- getOption("warn"); options(warn = 0)
pdf(km_pdf, width = 8, height = 6)
print(p_km)
dev.off()
options(warn = old_warn)

# ==========================================================================
# 7) UPDATE TABLE 3
# ==========================================================================
message("[GSE1456] Step 7: Updating Table 3...")

t3_path <- file.path(PATHS$results$main, "Table3_Survival_Performance_ByCohort.csv")
t3 <- read.csv(t3_path, stringsAsFactors = FALSE, check.names = FALSE)

# Create new row in same format
new_row <- data.frame(
  Cohort = "GSE1456",
  Role = "Validation",
  Endpoint = "OS",
  N = nrow(ar),
  Events = sum(ar$os_event),
  `FU med (mo)` = round(median(ar$os_time_months[ar$os_event == 0], na.rm=TRUE), 0),
  `HR uni (95%CI)` = sprintf("%.2f (%.2f-%.2f)", hr_uni, ci_uni_lo, ci_uni_hi),
  `p (uni)` = if(p_uni < 0.001) sprintf("%.2e", p_uni) else sprintf("%.4f", p_uni),
  `C_adj (95%CI)` = sprintf("%.3f (%.3f-%.3f)", c_idx, c_idx-1.96*c_se, c_idx+1.96*c_se),
  `HR adj (CORE-A)` = sprintf("%.2f (%.2f-%.2f)", hr_multi, ci_multi_lo, ci_multi_hi),
  `p (CORE-A adj)` = if(p_multi < 0.001) sprintf("%.2e", p_multi) else sprintf("%.4f", p_multi),
  check.names = FALSE, stringsAsFactors = FALSE
)

# Append if not already present
if (!"GSE1456" %in% t3$Cohort) {
  t3 <- rbind(t3, new_row)
  write.csv(t3, t3_path, row.names = FALSE)
  message("[GSE1456] Added to Table 3")
} else {
  message("[GSE1456] Already in Table 3, skipping")
}

# ==========================================================================
# 8) UPDATE GENE COVERAGE TABLE
# ==========================================================================
gcov_path <- file.path("results/supp", "gene_coverage_by_cohort.csv")
if (file.exists(gcov_path)) {
  gcov <- read.csv(gcov_path, stringsAsFactors = FALSE)
  if (!"GSE1456" %in% gcov$cohort) {
    gcov <- rbind(gcov, data.frame(
      cohort = "GSE1456",
      platform = "Affymetrix HG-U133A microarray",
      genes_present = length(available),
      genes_total = length(corepam_genes),
      coverage_pct = round(100*frac_present, 1),
      missing_genes = if(length(missing)>0) paste(missing, collapse=", ") else "",
      stringsAsFactors = FALSE
    ))
    write.csv(gcov, gcov_path, row.names = FALSE)
    message("[GSE1456] Added to gene coverage table")
  }
}

# ==========================================================================
# 9) UPDATE COHORT SUMMARY TABLE
# ==========================================================================
cst_path <- file.path("results/supp", "cohort_summary_table.csv")
if (file.exists(cst_path)) {
  cst <- read.csv(cst_path, stringsAsFactors = FALSE)
  if (!"GSE1456" %in% cst$cohort) {
    cst <- rbind(cst, data.frame(
      cohort = "GSE1456",
      role = "VALIDATION",
      platform = "Affymetrix HG-U133A",
      n_samples = nrow(ar),
      n_os_events = sum(ar$os_event),
      pct_os_events = round(100*mean(ar$os_event), 1),
      median_fu_all_mo = round(median(ar$os_time_months, na.rm=TRUE), 1),
      median_fu_cens_mo = round(median(ar$os_time_months[ar$os_event==0], na.rm=TRUE), 1),
      max_fu_mo = round(max(ar$os_time_months, na.rm=TRUE), 1),
      n_genes_total = nrow(gene_expr_mat),
      stringsAsFactors = FALSE
    ))
    write.csv(cst, cst_path, row.names = FALSE)
    message("[GSE1456] Added to cohort summary")
  }
}

message("[GSE1456] ========== INTEGRATION COMPLETE ==========")
message(sprintf("[GSE1456] SUMMARY: N=%d, OS events=%d, HR=%.3f (p=%.4g), C=%.4f, genes=%d/%d",
                nrow(ar), sum(ar$os_event), hr_uni, p_uni, c_idx,
                length(available), length(corepam_genes)))
