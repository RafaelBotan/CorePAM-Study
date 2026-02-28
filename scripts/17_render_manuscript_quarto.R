# =============================================================================
# SCRIPT: 17_render_manuscript_quarto.R
# PURPOSE: Render the Core-PAM manuscript Quarto document to PDF and HTML.
#          Requires all upstream scripts (01-16) to have completed.
#          Updates the RESULTS_SUMMARY.md after each successful render.
# PROJECT: Core-PAM (Memorial v6.1 §Manuscript)
#
# INPUTS:
#   manuscript/CorePAM_manuscript.qmd  (or any .qmd at ROOT_REPO level)
#   results/ and figures/ directories (all upstream outputs)
#
# OUTPUTS:
#   manuscript/CorePAM_manuscript.pdf
#   manuscript/CorePAM_manuscript.html
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "17_render_manuscript_quarto.R"

# Skip if outputs exist and QMD hasn't changed (set FORCE_RERUN=TRUE to rerender)
FORCE <- as.logical(Sys.getenv("FORCE_RERUN", "FALSE"))

# ---------------------------------------------------------------------------
# 1) Locate QMD file
# ---------------------------------------------------------------------------
qmd_candidates <- c(
  file.path(ROOT_REPO, "manuscript", "CorePAM_manuscript.qmd"),
  file.path(ROOT_REPO, "CorePAM_manuscript.qmd"),
  list.files(ROOT_REPO, pattern = "manuscript.*\\.qmd$",
             recursive = TRUE, full.names = TRUE)[1]
)
qmd_candidates <- qmd_candidates[!is.na(qmd_candidates)]

qmd_path <- qmd_candidates[file.exists(qmd_candidates)][1]

if (is.na(qmd_path)) {
  message("[17] No manuscript QMD found. Creating template scaffold...")

  # Create manuscript directory and template QMD
  ms_dir <- file.path(ROOT_REPO, "manuscript")
  dir.create(ms_dir, showWarnings = FALSE, recursive = TRUE)
  qmd_path <- file.path(ms_dir, "CorePAM_manuscript.qmd")

  template <- '---
title: "Core-PAM: A Minimum Breast Cancer Gene Panel Non-Inferior to PAM50 for Prognostic Risk Stratification"
author: "Rafael Botan"
date: today
format:
  pdf:
    documentclass: article
    toc: true
    number-sections: true
  html:
    toc: true
    toc-float: true
    embed-resources: true
bibliography: references.bib
---

```{r setup, include=FALSE}
source(here::here("scripts/00_setup.R"))
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
```

## Introduction

<!-- TODO: Introduction text -->

## Methods

<!-- TODO: Methods text — see METHODS_DRAFT.md for draft paragraphs -->

## Results

```{r load-card}
library(jsonlite)
card <- fromJSON(file.path(PATHS$results$corepam, "CorePAM_training_card.json"))
```

The Core-PAM panel consists of **`r card$n_genes`** genes selected from the PAM50 panel
via Cox elastic-net regression (α = `r card$alpha`, K = `r card$k_folds` stratified folds,
seed = `r card$seed_folds`).

```{r meta-results}
meta_df <- read_csv(file.path(PATHS$results$main, "meta_survival_summary.csv"),
                    show_col_types = FALSE)
uni_row <- meta_df[meta_df$analysis == "univariado", ]
```

The pooled hazard ratio across `r uni_row$n_cohorts` cohorts was
**HR = `r round(uni_row$meta_HR, 2)`
(95% CI: `r round(uni_row$meta_HR_lo95, 2)`–`r round(uni_row$meta_HR_hi95, 2)`)**.
Heterogeneity was I² = `r round(uni_row$I2_pct, 1)`%.

## Discussion

<!-- TODO: Discussion -->

## References

::: {#refs}
:::
'
  writeLines(template, qmd_path)
  message(sprintf("[17] Template QMD created: %s", qmd_path))
  message("[17] Fill in narrative sections before running again.")
}

message(sprintf("[17] QMD file: %s", qmd_path))

# Check if Quarto is available
quarto_check <- tryCatch(
  system2("quarto", "--version", stdout = TRUE, stderr = FALSE),
  error = function(e) NULL
)

if (is.null(quarto_check)) {
  message("[17] WARNING: Quarto not found in PATH.")
  message("[17] Install from https://quarto.org or add to PATH.")
  message("[17] Attempting via quarto R package...")
  if (!requireNamespace("quarto", quietly = TRUE)) {
    stop("[17] Quarto not available. Install quarto CLI or the R quarto package.")
  }
}

# Skip if PDF already exists and QMD not modified (unless FORCE)
pdf_out <- sub("\\.qmd$", ".pdf", qmd_path)
html_out <- sub("\\.qmd$", ".html", qmd_path)

if (!FORCE && file.exists(pdf_out)) {
  qmd_mtime <- file.info(qmd_path)$mtime
  pdf_mtime <- file.info(pdf_out)$mtime
  if (!is.na(pdf_mtime) && pdf_mtime > qmd_mtime) {
    message("[17] PDF is up to date. Set FORCE_RERUN=TRUE to re-render.")
    quit(save = "no", status = 0)
  }
}

# ---------------------------------------------------------------------------
# 2) Render to PDF
# ---------------------------------------------------------------------------
message("[17] Rendering manuscript to PDF...")

render_ok <- tryCatch({
  if (requireNamespace("quarto", quietly = TRUE)) {
    quarto::quarto_render(qmd_path, output_format = "pdf",
                          quiet = FALSE)
  } else {
    ret <- system2("quarto", c("render", shQuote(qmd_path),
                                "--to", "pdf"),
                   stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(ret, "status")) && attr(ret, "status") != 0) {
      stop(paste(ret, collapse = "\n"))
    }
  }
  TRUE
}, error = function(e) {
  message(sprintf("[17] PDF render failed: %s", e$message))
  FALSE
})

if (render_ok && file.exists(pdf_out)) {
  h_pdf <- sha256_file(pdf_out)
  sz_pdf <- file.info(pdf_out)$size / 1e6
  registry_append("ALL", "manuscript_pdf", pdf_out, h_pdf, "ok",
                  SCRIPT_NAME, sz_pdf)
  message(sprintf("[17] PDF rendered: %s (%.2f MB)", pdf_out, sz_pdf))
} else {
  message("[17] WARNING: PDF output not found after render attempt.")
}

# ---------------------------------------------------------------------------
# 3) Render to HTML
# ---------------------------------------------------------------------------
message("[17] Rendering manuscript to HTML...")

render_ok_html <- tryCatch({
  if (requireNamespace("quarto", quietly = TRUE)) {
    quarto::quarto_render(qmd_path, output_format = "html",
                          quiet = FALSE)
  } else {
    ret <- system2("quarto", c("render", shQuote(qmd_path),
                                "--to", "html"),
                   stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(ret, "status")) && attr(ret, "status") != 0) {
      stop(paste(ret, collapse = "\n"))
    }
  }
  TRUE
}, error = function(e) {
  message(sprintf("[17] HTML render failed: %s", e$message))
  FALSE
})

if (render_ok_html && file.exists(html_out)) {
  h_html <- sha256_file(html_out)
  sz_html <- file.info(html_out)$size / 1e6
  registry_append("ALL", "manuscript_html", html_out, h_html, "ok",
                  SCRIPT_NAME, sz_html)
  message(sprintf("[17] HTML rendered: %s (%.2f MB)", html_out, sz_html))
}

# ---------------------------------------------------------------------------
# 4) Append to RESULTS_SUMMARY.md
# ---------------------------------------------------------------------------
append_results_summary <- function(entry) {
  rs_path <- file.path(ROOT_REPO, "RESULTS_SUMMARY.md")
  entry_text <- sprintf(
    "\n---\n**Date:** %s | **Script:** %s\n%s\n",
    format(Sys.time(), "%Y-%m-%d %H:%M"),
    SCRIPT_NAME,
    entry
  )
  cat(entry_text, file = rs_path, append = TRUE)
}

append_results_summary(sprintf(
  "Manuscript rendered. PDF: %s | HTML: %s | Status: PDF=%s HTML=%s",
  if (file.exists(pdf_out)) basename(pdf_out) else "MISSING",
  if (file.exists(html_out)) basename(html_out) else "MISSING",
  if (render_ok) "OK" else "FAIL",
  if (render_ok_html) "OK" else "FAIL"
))

if (!render_ok && !render_ok_html) {
  stop("[17] Both PDF and HTML render failed. Check Quarto installation and QMD syntax.")
}

message("[17] Manuscript render complete.")
message("[17] Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
