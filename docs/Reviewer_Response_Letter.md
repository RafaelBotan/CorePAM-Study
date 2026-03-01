# Response to Reviewers
**Journal:** Breast Cancer Research
**Manuscript:** "Derivation and Multicohort Validation of Core-PAM: a Minimal Prognostic
Gene Signature from the PAM50 Panel in Breast Cancer"
**Date:** March 2026

---

## Cover Letter

To the Editors and Reviewers of *Breast Cancer Research*,

We would like to open this letter not with the usual formalities, but with something we
genuinely mean: **thank you**.

When we received the reviewers' comments on our original manuscript, we will not pretend
we were not initially discouraged — the critiques were pointed, and some of them struck at
the core of our methodology. But after carefully re-reading each point, what struck us most
was how insightful, rigorous, and — ultimately — generous the feedback truly was. The
reviewers did not simply criticize the work; they illuminated exactly where it fell short
and, in doing so, gave us a precise map for how to do it right.

We took that map seriously. In fact, we took it so seriously that **we rebuilt the entire
study from the ground up**.

---

### What changed, and why

Our original manuscript used a gene-intersection approach across microarray platforms to
identify a 40-gene subset of PAM50. While internally consistent, the reviewers rightly
noted that this criterion was largely subjective and that the resulting model showed poor
transferability to independent RNA-seq cohorts — a C-index of 0.42 in TCGA-BRCA, which,
frankly, is barely above chance and should have given us pause earlier. Endpoints were
mixed (OS and DRFS used interchangeably), there were no clinical covariates included in
the prognostic model, and the validation cohorts were all microarray-based, leaving the
RNA-seq generalizability question entirely unanswered.

Each of these points represented a genuine methodological weakness, and we are grateful
the reviewers had the rigor to name them clearly.

In the revised study, **every one of those concerns has been directly addressed**:

- **Gene selection is now fully data-driven.** We replaced the intersection criterion with
  an out-of-fold (OOF) non-inferiority test: we selected the *smallest* gene subset whose
  OOF C-index on the training cohort (SCAN-B) falls within ΔC = 0.010 of the full PAM50
  model. The number of genes was *derived*, not pre-specified or targeted. This approach
  is pre-specified, deterministic (based on SHA-256 patient-ID stratification, not random
  seeds), and fully reproducible.

- **The training cohort is now RNA-seq.** SCAN-B (GSE96058, N = 3,069 patients,
  322 OS events) was used as the discovery cohort — the largest publicly available
  breast cancer RNA-seq cohort with long-term follow-up. This fundamentally changed the
  transferability profile of the model.

- **Endpoints are now homogeneous.** All analyses are restricted to OS and DSS
  (disease-specific survival for METABRIC), eliminating the endpoint mixing that
  complicated interpretation in the original work.

- **Clinical covariates are explicitly modelled.** A CORE-A specification
  (CorePAM score + age + ER status) is reported alongside the univariate score, consistent
  with current prognostic modelling standards and with what a clinician would realistically
  use in practice.

- **Three independent external validation cohorts** were used: TCGA-BRCA (RNA-seq,
  N = 1,072), METABRIC (Illumina microarray, N = 1,978), and GSE20685 (Affymetrix
  microarray, N = 327). A secondary analysis block adds pCR prediction across five
  additional cohorts (N = 1,683 total), including I-SPY2 (N = 986) as an exploratory
  external dataset.

Perhaps most importantly — and we find this scientifically satisfying — the revised
methodology led us to a **24-gene panel**, smaller than the original 40-gene version, and
yet one that performs substantially better across all independent validation cohorts.
The reviewers' criticism did not just improve the paper; it led to a genuinely better
scientific result.

---

### Reproducibility and audit

We want to be transparent: the complete pipeline — from raw data download (GEO + GDC +
cBioPortal) through all statistical analyses, figure generation, and manuscript
compilation — is version-controlled and publicly available on GitHub. Every output file
carries a SHA-256 hash registered in a study registry, and the pipeline is designed for
full computational reproducibility by any independent reviewer. We believe this level of
transparency is what the reviewers' rigorous standards deserve, and we hope it also
serves as a useful resource for the community.

---

We address each reviewer comment individually in the point-by-point section below,
noting the specific changes made and the section/table/figure where each change can be
verified in the revised manuscript.

We are deeply grateful to the reviewers for the time and intellectual effort they
invested in our work. Their feedback did not merely improve a manuscript — it reshaped
an entire research program. We remain entirely open to further scrutiny and welcome any
additional questions or concerns.

Respectfully,

**Rafael Botan**
*PhD Candidate, Genomic Medicine*

---
---

## Point-by-Point Responses

> **Notation:** Reviewer comments are in italic. Our responses are in plain text.
> All numbers cited are drawn directly from results files; no values are hard-coded.

---

### Reviewer 1

---

**Reviewer 1, Comment 1:**
*"The model shows a C-index of 0.42 in TCGA-BRCA, which is below what would be expected
even from a random predictor with favourable event distribution. The authors should
address the poor cross-platform transferability of the signature."*

**Response:**

We thank the reviewer for this critical observation — it was the most consequential single
point in the entire review, and it was absolutely correct. A C-index of 0.42 in an
independent RNA-seq cohort signals not a modestly underperforming model, but a
fundamentally misaligned one.

The root cause was that our original model was trained exclusively on Illumina microarray
data (METABRIC) and then applied to RNA-seq without any accommodation for the
distributional differences between platforms. The gene-intersection criterion we used was
designed to ensure probe coverage, not platform-invariant signal.

In the revised manuscript, we addressed this at the source. **SCAN-B (GSE96058,
N = 3,069, RNA-seq, 322 OS events) is now the training cohort.** The Core-PAM score is
derived entirely on RNA-seq data using Cox elastic-net regularisation (α = 0.5, K = 10
folds), and validated on both an independent RNA-seq cohort (TCGA-BRCA) and two
independent microarray cohorts (METABRIC, GSE20685). The C-index values in the revised
study are:

| Cohort | Platform | N | Endpoint | C-index (95% CI) |
|---|---|---|---|---|
| SCAN-B *(training)* | RNA-seq | 3,069 | OS | 0.698 (0.668–0.729) |
| TCGA-BRCA *(validation)* | RNA-seq | 1,072 | OS | 0.624 (0.572–0.678) |
| METABRIC *(validation)* | Microarray | 1,978 | DSS | 0.638 (0.615–0.659) |
| GSE20685 *(validation)* | Microarray | 327 | OS | 0.623 (0.562–0.683) |

The validation C-index of **0.624 in TCGA-BRCA** represents a complete reversal of the
previous finding (0.42), and is now consistent with the discrimination observed across all
independent cohorts. We believe this directly resolves the transferability concern.

---

**Reviewer 1, Comment 2:**
*"The gene selection criterion — intersection of genes with probes across platforms — is
subjective and not statistically grounded. There is no principled reason to stop at 40
genes."*

**Response:**

This criticism is entirely valid and we thank the reviewer for articulating it so
precisely. The intersection criterion guaranteed probe coverage, but it did not answer
the scientifically meaningful question: *how many genes do you actually need?*

In the revised study, we reframe the gene selection problem as a non-inferiority
decision. Using OOF C-index curves across the regularisation path of a Cox elastic-net
model (α = 0.5, K = 10 stratified folds on SCAN-B), we select the **smallest** gene
subset whose OOF C-index is within ΔC = 0.010 of the full PAM50 model. The threshold
ΔC = 0.010 was pre-specified in the study protocol (*analysis_freeze.csv*) and aligns
with published benchmarks for clinically non-inferior discrimination.

This approach yielded a **24-gene Core-PAM panel** — smaller than the original 40-gene
version, derived without subjective judgment, and with a clearly interpretable stopping
criterion. The gene count is a result, not an input.

The panel, weights, and all regularisation parameters are frozen, version-controlled,
and verified by SHA-256 hash.

---

**Reviewer 1, Comment 3:**
*"No clinical covariates are included in the model. Age, ER status, and HER2 would be
expected confounders in any breast cancer prognostic analysis."*

**Response:**

We agree completely. The original manuscript reported only the univariate score, which is
insufficient for clinical interpretation and does not address the meaningful question of
whether the gene signature provides prognostic value *beyond* standard clinical variables.

In the revised manuscript, we report a **CORE-A specification** (CorePAM score + age +
ER status) alongside the univariate model in every survival cohort. Age-adjusted hazard
ratios are shown in Table 3 of the revised manuscript and remain consistent across all
cohorts. An incremental value analysis (ΔC-index relative to the CORE-A clinical-only
baseline) is reported in Section 3.3 and Supplementary Figure S5/S6, demonstrating that
the CorePAM score provides significant additive discrimination in all four validation
cohorts (ΔC range: +0.030 to +0.142).

We also note that Fine–Gray competing-risk models (treating non-breast-cancer deaths as
competing events) are reported for METABRIC DSS as a sensitivity analysis.

---

### Reviewer 2

---

**Reviewer 2, Comment 1:**
*"The authors mix OS and DRFS as if they were interchangeable endpoints. These have
fundamentally different biological and clinical meanings, and pooling them in a
meta-analysis is inappropriate."*

**Response:**

This is a very fair criticism and one we should have caught ourselves. Overall survival
and distant relapse-free survival are not interchangeable, and treating them as such
obscures the biological signal one is trying to measure.

In the revised manuscript, endpoints are fully harmonised:
- **OS** for SCAN-B, TCGA-BRCA, and GSE20685
- **DSS (disease-specific survival)** for METABRIC, where OS is reported as a
  sensitivity analysis

The meta-analysis of survival effects (Section 3.2, Figure 4) is restricted to OS
across the three validation cohorts with OS as primary endpoint; METABRIC contributes
its DSS result. The heterogeneity (I² = 90.6%) is discussed explicitly in the manuscript
and is expected given the multi-platform, multi-institution nature of the cohorts.

DRFS is no longer used in any primary analysis.

---

**Reviewer 2, Comment 2:**
*"Validation is performed only on microarray cohorts. The authors trained on microarray
(METABRIC) and validated on microarray. An RNA-seq validation cohort is essential to
demonstrate platform generalisability."*

**Response:**

We thank the reviewer for this point — it was the impetus for the most fundamental
redesign in the revised manuscript.

The training cohort is now **SCAN-B (RNA-seq, N = 3,069)**. Validation is performed on:
- **TCGA-BRCA (RNA-seq, N = 1,072)** — independent RNA-seq validation
- **METABRIC (Illumina microarray, N = 1,978)** — microarray validation
- **GSE20685 (Affymetrix microarray, N = 327)** — independent microarray validation

The signal is consistent across both platforms, with C-index values ranging from 0.623
to 0.638 in validation cohorts. This directly demonstrates the platform generalisability
the reviewer requested.

---

**Reviewer 2, Comment 3:**
*"The sample sizes in some reported cohorts appear inconsistent between the text and
supplementary tables. This undermines confidence in the reported statistics."*

**Response:**

We apologise for this inconsistency in the original submission. In the revised
manuscript, **all reported numbers are read directly from CSV/JSON results files at
manuscript compilation time**; no values are hard-coded in the text. Every table and
figure references a specific results file with a SHA-256 hash, making post-hoc
discrepancies structurally impossible.

The complete inclusion/exclusion flow for each cohort is documented in the supplementary
methods and in machine-readable inclusion logs (e.g., *ISPY2_inclusion_log.csv* for the
I-SPY2 exploratory cohort, which documents the transition from 988 enrolled patients to
986 with matched expression and clinical data).

---

### General Comments

---

**General Comment 1:**
*"The manuscript would benefit from a more principled assessment of discriminative
performance, including confidence intervals for C-index and formal tests of incremental
value."*

**Response:**

All C-index estimates in the revised manuscript are accompanied by **bootstrap 95%
confidence intervals** (B = 1,000 resamples). Incremental value (ΔC-index relative to
CORE-A clinical baseline) is tested using paired bootstrap, with results summarised in
Table S3. Calibration (intercept and slope from logistic calibration at the median
follow-up per cohort) and Brier scores are reported in Supplementary Figure S5/S6 for
all validation cohorts.

---

**General Comment 2:**
*"The secondary pCR analysis would be strengthened by inclusion of a larger, prospective
cohort if available."*

**Response:**

We share this aspiration. In response, we added **I-SPY2 (GSE194040, N = 986,
pCR = 32.4%)** as an exploratory external validation cohort for the pCR prediction
block. I-SPY2 is a multi-arm prospective neoadjuvant trial, representing the closest
available public dataset to a prospective external validation. The result (univariate
OR = 1.685, 95% CI: 1.454–1.953, AUC = 0.648; adjusted OR = 1.485, 95% CI:
1.244–1.771, AUC = 0.738) is consistent with the primary GEO cohort meta-analysis
(pooled OR = 1.686, I² = 0%, k = 4 cohorts). We classify I-SPY2 as *exploratory* to
maintain the pre-specified primary/secondary distinction of our protocol.

---

*We once again thank the reviewers for their exceptional engagement with this work.
The revised manuscript is, without question, a better scientific product because of it.*

---

*End of Response Letter*
