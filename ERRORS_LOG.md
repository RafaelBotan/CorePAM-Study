# ERRORS_LOG.md — Core-PAM Pipeline Error Log

**Convention:** Before fixing any script error, append an entry here with:
- Date/time
- Script name
- Error message (verbatim or condensed)
- Root cause diagnosis
- Solution applied

---

## 2025-XX-XX | 00_setup.R | sprintf format error

**Script:** `00_setup.R`
**Error:**
```
Error in sprintf("%.3f", FREEZE$delta_c) :
  invalid format '%.3f'; use format %s for character objects
```
**Root cause:** `type.convert(as.character(...))` on a heterogeneous vector
returned a character class for the column, so numeric parameters were stored
as `"character"` type instead of `numeric`. `sprintf("%.3f", ...)` requires a
numeric argument.

**Solution:** Replaced `type.convert` with `lapply` + `suppressWarnings(as.numeric(x))`,
falling back to character if conversion fails. Now each element is individually
coerced before being stored in the `FREEZE` named list.

```r
# FIX: individually coerce each value to numeric where possible
FREEZE <- setNames(
  lapply(.freeze_raw$value, function(x) {
    num <- suppressWarnings(as.numeric(x))
    if (!is.na(num)) num else x
  }),
  .freeze_raw$parameter
)
```

---

## 2025-XX-XX | 03_expression_preprocess_TCGA_BRCA.R | unexpected 'else'

**Script:** `03_expression_preprocess_TCGA_BRCA.R`
**Error:**
```
Error: unexpected 'else' in "                   else"
```
**Root cause:** R's parser requires `if-else` to be on a single line when used
as an expression assigned to a variable. Splitting across two lines causes the
parser to treat the first line as a complete statement, making `else` orphaned.

**Solution:** Collapsed the two lines into one:
```r
# BROKEN (two lines):
count_assay <- if ("unstranded" %in% assayNames(se_obj)) "unstranded"
               else assayNames(se_obj)[1]

# FIXED (one line):
count_assay <- if ("unstranded" %in% assayNames(se_obj)) "unstranded" else assayNames(se_obj)[1]
```

---

## 2026-02-28 | 01_download_raw_data.R | GEOquery load fails under warn=2

**Script:** `01_download_raw_data.R`
**Error:**
```
Error: package or namespace load failed for 'GEOquery' in (function (n) :
 (converted from warning) internal error 1 in R_decompress1 with libdeflate
```
**Root cause:** Two overlapping issues:
1. `source("scripts/00_setup.R")` sets `options(warn=2)` (strict mode).
   `suppressPackageStartupMessages()` suppresses `message()` calls but NOT `warning()` calls.
   So internal warnings during package loading are converted to hard errors by `warn=2`.
2. GEOquery has an internal `libdeflate` decompression warning on this R installation,
   which is normally harmless but fatal under `warn=2`.

**Solution:**
- Wrap all `library()` calls in `01_download_raw_data.R` with `old_warn / options(warn=0) / restore`
  pattern — same pattern used in all other scripts for external package operations.
- Also reinstall GEOquery to clear the libdeflate issue:
  `BiocManager::install("GEOquery", force = TRUE)`

```r
# FIXED:
old_warn <- getOption("warn"); options(warn = 0)
suppressPackageStartupMessages({
  library(GEOquery)
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(httr)
  library(curl)
})
options(warn = old_warn)
```

---

## 2026-02-28 | 01_download_raw_data.R | HTTP 403 downloading METABRIC from cBioPortal

**Script:** `01_download_raw_data.R`
**Error:**
```
HTTP 403 downloading: https://cbioportal-datahub.s3.amazonaws.com/brca_metabric.tar.gz
```
**Root cause:** cBioPortal migrated their data hosting from `cbioportal-datahub.s3.amazonaws.com`
to a new CDN at `datahub.assets.cbioportal.org`. The old S3 bucket is no longer publicly
accessible (returns HTTP 403 Forbidden).

**Solution:** Update `METABRIC_BASE_URL` to use the new domain. Add old URL as commented
fallback for reference. The file path after the domain remains the same (`brca_metabric.tar.gz`).

```r
# FIXED:
METABRIC_BASE_URL <- "https://datahub.assets.cbioportal.org/brca_metabric.tar.gz"
# Old URL (403 as of 2026-02): "https://cbioportal-datahub.s3.amazonaws.com/brca_metabric.tar.gz"
```

**Reference:** cBioPortal datasets page (`https://www.cbioportal.org/datasets`) confirms
`study_download_url: "https://datahub.assets.cbioportal.org/"` as current base URL.

---

<!-- Add new errors below this line, most recent first -->
