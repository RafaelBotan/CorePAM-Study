# Audit: Pipeline Changes and New Findings for R2 Revision

**Date:** 2026-04-12
**Context:** CorePAM manuscript revision (Breast Cancer Research, Round 2)

---

## 1. Bug Fix: GSE20685 Platform Misannotation

### What was wrong
Script `03_expression_preprocess_GSE20685.R` used case-sensitive regex patterns against `tolower()` input, causing the Affymetrix HG-U133 Plus 2.0 platform (GPL570, ~54k probes) to fall through to the HGU133A annotation database (~22k probes).

### Root cause
```r
# BEFORE (broken): uppercase patterns never match tolower() input
platform_pkg <- dplyr::case_when(
  grepl("GPL96|hgu133a|133a$",  tolower(platform)) ~ "hgu133a",    # <-- matched first!
  grepl("GPL570|hgu133plus2",   tolower(platform)) ~ "hgu133plus2", # <-- never reached
  ...
)
```

### Fix applied
```r
# AFTER (correct): lowercased patterns, GPL570 checked first
platform_pkg <- dplyr::case_when(
  grepl("gpl570|hgu133plus2",   tolower(platform)) ~ "hgu133plus2",
  grepl("gpl96|hgu133a|133a$",  tolower(platform)) ~ "hgu133a",
  ...
)
```

### Impact on results

| Metric | Before (HGU133A) | After (HG-U133 Plus 2.0) | Change |
|--------|-------------------|---------------------------|--------|
| Probes mapped | 12,970 | 20,529 | +7,559 |
| Unique HGNC genes | 9,397 | 13,889 | +4,492 |
| CorePAM coverage | 22/24 (91.7%) | **24/24 (100%)** | +2 genes |
| Missing genes | GPR160, CXXC5 | none | Fixed |
| HR univariate | 1.40 (1.12-1.76) | 1.40 (1.12-1.75) | Negligible |
| C-index (Harrell) | 0.623 (0.562-0.683) | 0.626 (0.565-0.684) | +0.003 |
| HR adjusted (age) | 1.40 (1.12-1.76) | 1.40 (1.12-1.75) | Negligible |
| C frozen z-score | 0.581 | 0.591 | +0.010 |
| Delta-C frozen-intra | -0.042 | -0.036 | Improved |
| r (frozen vs intra) | 0.921 | 0.923 | Stable |
| Meta-analysis HR (RE) | 1.373 (I^2=38.3%) | 1.373 (I^2=38.2%) | Negligible |

**Conclusion:** Factually important correction (platform name, gene coverage), but quantitative impact is minimal. All conclusions hold. The correction actually strengthens the frozen z-score performance (smaller gap: -0.036 vs -0.042).

### GSE1456 NOT affected
GSE1456 correctly uses GPL96 / HGU133A. Verified by regex: `grepl("gpl96", tolower("GPL96"))` = TRUE. No changes needed.

---

## 2. New Analysis: Head-to-Head External Comparison (24 vs 49 genes)

### What was done
Script `38_R2_headtohead_24vs50.R` compares CorePAM (24 genes) against the PAM50-full elastic-net model (49 genes with non-zero weights at max-df lambda) in all 4 external validation cohorts.

Two scoring modes:
- **Intra-cohort z-score** (framework-matched): each cohort z-scored independently
- **Frozen z-score** (deployment-matched): uses SCAN-B reference mean/SD

### Key finding: NDC80 has zero coefficient
The elastic-net Cox model retains 49 of 50 PAM50 genes at the maximum-df lambda. NDC80 is shrunk to exactly zero. This is expected elastic-net behavior (L1 penalty).

### Results

| Cohort | Mode | C_24 | C_50 | Delta-C (95% CI) | Verdict |
|--------|------|------|------|-------------------|---------|
| TCGA-BRCA | intra | 0.624 | 0.631 | -0.008 (-0.038, +0.022) | Non-inferior (CI includes 0) |
| TCGA-BRCA | frozen | 0.624 | 0.634 | -0.010 (-0.039, +0.020) | Non-inferior (CI includes 0) |
| METABRIC | intra | 0.638 | 0.616 | **+0.022** (+0.011, +0.032) | **24-gene SUPERIOR** |
| METABRIC | frozen | 0.638 | 0.622 | **+0.017** (+0.007, +0.027) | **24-gene SUPERIOR** |
| GSE20685 | intra | 0.626 | 0.571 | **+0.055** (+0.010, +0.101) | **24-gene SUPERIOR** |
| GSE20685 | frozen | 0.591 | 0.544 | **+0.047** (+0.007, +0.085) | **24-gene SUPERIOR** |
| GSE1456 | intra | 0.661 | 0.564 | **+0.097** (+0.032, +0.166) | **24-gene SUPERIOR** |
| GSE1456 | frozen | 0.622 | 0.510 | **+0.111** (+0.045, +0.181) | **24-gene SUPERIOR** |

### Interpretation
- In 3/4 external cohorts, the 24-gene CorePAM is **statistically superior** to the 49-gene model (95% CI entirely above zero).
- Only TCGA-BRCA shows a slight edge for the 50-gene model (Delta-C = -0.008 to -0.010), but the CI comfortably includes zero.
- This is consistent with regularization theory: sparser models generalize better cross-platform. The Pareto frontier (OOF C-index vs df) already showed diminishing returns beyond df=24.
- Pearson correlation between 24-gene and 49-gene scores ranges from 0.62 (GSE1456 frozen) to 0.91 (METABRIC frozen), confirming they capture overlapping but not identical signal.

### Output files
- `results/supp/headtohead_24vs50_external.csv` (data table)
- `figures/supp/en/pdf/FigS_HeadToHead_24vs50.pdf` (forest plot)
- `figures/supp/en/png/FigS_HeadToHead_24vs50.png` (forest plot)

---

## 3. Multi-Panel KM Figure (Script 37)

### What was done
Replaced 5 individual KM figures (old Figures 2-6) with a single consolidated multi-panel figure (new Figure 2, panels A-E).

### Changes
- Each panel: one cohort, dichotomised at intra-cohort median
- Single log-rank p-value + HR annotation per panel
- Consistent colour scheme and formatting
- 3x2 grid layout (5 panels + methodology note)

### Output
- `figures/main/en/pdf/Fig2_KM_MultiPanel_EN.pdf`
- `figures/main/en/png/Fig2_KM_MultiPanel_EN.png`

---

## 4. Supplementary Renumbering

All supplementary materials renumbered to BCR's unified "Additional file" system (AF1-AF17), with physical reordering in the LaTeX source. Cross-references updated throughout.

---

## 5. Summary of Manuscript Changes Needed

### Factual corrections (must change)
1. Line 75: "Affymetrix HGU133A" --> "Affymetrix HG-U133 Plus 2.0 (GPL570)"
2. Line 115: Update platform description for GSE20685 (hgu133plus2.db, not hgu133a.db)
3. Line 117: Remove "GSE20685: 22/24 genes, 91.7%" example (now 24/24 = 100%)
4. Line 128 (Table): GSE20685 row: platform, coverage, missing genes
5. Lines 229, 285, 575, 591, 606, 645: Update GSE20685 numbers (minimal changes)
6. Line 306: Update Uno's C and tdAUC for GSE20685

### New content (add)
7. Head-to-head comparison results (new supplementary table + figure)
8. GSE1456 row in frozen z-score table (currently omitted)
9. Brief mention of head-to-head results in Discussion

### NOT changing (confirmed by researcher)
- DCA remains supplementary (not promoted to main)
- Title keeps "cross-platform"
- No changes to methodology, scoring formula, or statistical framework
- No tone/framing changes beyond what's factually necessary
- Existing text structure and flow preserved

---

## 6. Figure Integrity Audit

### CURRENT (updated 2026-04-12, consistent with corrected pipeline)
- Fig2_KM_GSE20685_OS_EN/PT (pdf/png) — GSE20685 KM curves
- Fig4_Meta_Forest_HR_per1SD_CorePAM (pdf/png) — meta-analysis forest plot
- Fig5_DeltaCindex_COREA_vs_COREAplus_CorePAM (pdf/png) — incremental C
- Fig5_Calibration_60m_Panels_CorePAM (pdf/png) — calibration curves
- Fig2_KM_MultiPanel_EN (pdf/png) — new consolidated KM figure
- FigS_COREA_Sensitivity_EN (pdf/png) — frozen sensitivity
- FigS_DropGene_Tornado_EN (pdf/png) — leave-one-gene-out
- FigS_KM_GSE20685_OS_Quartis (EN/PT, pdf/png) — quartile KM
- FigS_HeadToHead_24vs50 (pdf/png) — NEW head-to-head forest plot

### REGENERATED (2026-04-12, previously stale)
- FigS_Heatmap_GSE20685 (pdf/png) — regenerated via 07x_extra_figures.R
- FigS_Forest_HR_ValidationCohorts (pdf/png) — regenerated via 07x_extra_figures.R
- FigS_ER_Stratified_Forest (pdf/png) — regenerated via 34_subgroup_er_stratified.R
- Fig4_Meta_Forest_HR_per1SD_CorePAM (pdf/png) — regenerated via 08_meta_survival.R
- Fig5_Calibration_60m_Panels_CorePAM (pdf/png) — regenerated via 11_incremental_value_and_dca.R
- Fig5_DeltaCindex_COREA_vs_COREAplus_CorePAM (pdf/png) — regenerated via 11_incremental_value_and_dca.R

**All figures containing GSE20685 data are now consistent with the corrected pipeline. No stale outputs remain.**

### UNAFFECTED (~74 files)
All pCR figures, study design, other cohort KM curves, DCA, and distribution figures are unaffected by the GSE20685 fix.

---

## 7. Registry and Traceability Audit

### Registry integrity: CONSISTENT
- All GSE20685 entries updated with new SHA256 hashes (2026-04-12)
- Gene coverage correctly shows 24/24 (was 22/24)
- Platform annotation: HG-U133 Plus 2.0 in all audit files
- Downstream files (survival, meta-analysis, figures) all have updated hashes
- Timestamps form consistent sequence: 15:00:10 -> 15:00:30 -> 15:01:08+

### Key registry entries verified
| Output | Old hash | New hash | Status |
|--------|----------|----------|--------|
| Gene_Mapping_Audit_GSE20685 | 1039b559... | eaf6cd3a... | Updated |
| Expression_preZ_GSE20685 | 9b69cc06... | 1bb926c4... | Updated |
| Analysis_ready_GSE20685 | 34121fe6... | d4b75d0b... | Updated (22->24 genes) |

### PAM50 gene coverage audit
- pam50_coverage_summary.csv: GSE20685 shows 50/50 PAM50 genes present, GO status
- pam50_gene_list_audit.csv: all 50 genes TRUE for GSE20685
- gene_mapping_audit_GSE20685.csv: complete audit trail with corrected platform

### Frozen z-score table
- GSE1456 successfully added (was previously omitted)
- All 5 cohorts now represented: SCANB, TCGA, METABRIC, GSE20685, GSE1456

---

## 8. Manuscript v3 Edits (2026-04-12)

### Factual corrections (completed)
- GSE20685 platform: HGU133A → HG-U133 Plus 2.0 (GPL570) — all references
- Gene coverage: 22/24 → 24/24 — tables and text
- All GSE20685 numerical results updated (HR, C-index, calibration, Uno's C, tdAUC)
- Meta-analysis: I²=38.3→38.2%, p=1.8e-9→1.6e-9
- Frozen ΔC: -0.042 → -0.036 in Discussion
- Calibration slope: 1.58 → 1.57 in Discussion
- Additional files: 1-17 → 1-18 (new Table S8 head-to-head)

### Editorial changes (R1/R2 response)
- Background: "demonstrated" → "is well established" (PAM50 original)
- Results: Head-to-head frase reescrita (comparador explicado, "non-inferior" per-coorte removido)
- Discussion opening: "We demonstrate" → anchor no framework + head-to-head integrado + caveat
- Limitations: reestruturado em 3 blocos (study design, analytical choices, clinical translation)
- Conclusions: non-inferiority ancorada ao framework, head-to-head mencionado, caveat de deployment
- Abstract Conclusions: "predicts pCR" → "associated with pCR", "resource-limited" → "simplify assay development + validation required"

### New content
- Additional file 18: Table S8 (head-to-head 24 vs 49 genes)
- 1 sentence in Results about head-to-head
- GSE1456 row in frozen z-score table

---

## 9. Submission Timeline

- **Deadline:** 2026-04-26 (14 days from audit date)
- **Requirements:** Point-by-point PDF, clean manuscript (no tracked changes), re-upload changed files
- **Status:** Manuscript v3 complete, all figures regenerated, point-by-point letter drafted
