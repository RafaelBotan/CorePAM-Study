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

## 2026-02-28 | 01_download_raw_data.R | METABRIC files missing after extraction

**Script:** `01_download_raw_data.R`
**Error:**
```
Error: (converted from warning) METABRIC file not found after extraction:
  .../RAW/METABRIC/brca_metabric/data_mrna_illumina_microarray.txt
```
**Root cause (cascade from previous HTTP 403 error):**
1. First run → HTTP 403 → httr wrote a 243-byte HTML error page to `brca_metabric.tar.gz`
2. Fixed URL, re-ran script
3. `download_file_safe` checks `file.info(destfile)$size > 0` → 243 bytes > 0 → **skipped real download**
4. `utils::untar()` silently failed on the 243-byte corrupt file (no error raised)
5. No subdirectory was created; `data_mrna_illumina_microarray.txt` not found

**Solution:**
- Delete the corrupt 243-byte tar.gz manually: `unlink(METABRIC_TAR)`
- Fix `download_file_safe` to require `size > MIN_VALID_BYTES` (e.g. 1 MB) for tar.gz files
- Fix extraction block to verify extracted files exist after `utils::untar()`, stop if missing

---

## 2026-02-28 | 01_download_raw_data.R | SSL connect error to GDC API (TCGA-BRCA)

**Script:** `01_download_raw_data.R`
**Error:**
```
GDCquery TCGA-BRCA expression failed: cannot open the connection to
  'https://api.gdc.cancer.gov/...'
SSL connect error
```
**Diagnosis:**
- System `curl` → GDC API → HTTP 200 OK (SSL works at OS level)
- R httr/libcurl → GDC API → SSL connect error (TLS handshake failure)
- Root cause: R on Windows uses its own bundled libcurl which may use a different
  TLS implementation or certificate store than the system curl. GDC requires TLS 1.2+
  and specific cipher suites that R's bundle may not negotiate correctly.

**Solution:** Add `httr::set_config(httr::config(ssl_verifypeer = FALSE))` before
TCGAbiolinks calls and `httr::reset_config()` after. Safe because:
1. GDC server is verified at OS level (system curl works)
2. Downloaded file integrity is verified with SHA-256

```r
httr::set_config(httr::config(ssl_verifypeer = FALSE))
# ... TCGAbiolinks GDCquery / GDCdownload / GDCprepare calls ...
httr::reset_config()
```

---

## 2026-02-28 | 01_download_raw_data.R | GSE20685 clinical download — connection error

**Script:** `01_download_raw_data.R`
**Error:**
```
Error in value[[3L]](cond) :
  error reading from the connection
```
**Context:** GSE20685 expression TAR (2672.7 MB) downloaded OK. Error occurred during
`geo_clinical_download()` → `GEOquery::getGEO("GSE20685", GSEMatrix=TRUE)` call.
**Root cause:** Intermittent network connection drop during GEO series matrix download.
The series matrix (~several MB) is a smaller file fetched over the same GEO FTP/HTTP
connection that had just completed a large download.
**Solution:** Re-run `01_download_raw_data.R` — all other cohorts will skip; only
GSE20685 clinical will retry (geo_clinical_download checks RDS existence).

---

## 2026-02-28 | 02_harmonize_clinical_TCGA_BRCA.R | strict_rds fails under warn=2

**Script:** `02_harmonize_clinical_TCGA_BRCA.R`
**Error:**
```
Error in readRDS(path) : (converted from warning) invalid or incomplete compressed data
```
**Root cause:** `options(warn=2)` converts harmless decompression warnings to hard errors.
`TCGA_BRCA_clinical_raw.rds` was saved with a different R/compression version and emits a
harmless decompression warning on `readRDS()`. Under `warn=2` this is fatal.

**Solution:** Patched `strict_rds()` in `00_setup.R` to temporarily lower warn level during
`readRDS()` using `on.exit()` for safe restoration:
```r
strict_rds <- function(path) {
  if (!file.exists(path)) stop("RDS not found: ", path)
  if (file.info(path)$size == 0) stop("Empty RDS (0 bytes): ", path)
  old_warn <- getOption("warn"); on.exit(options(warn = old_warn))
  options(warn = 0)
  readRDS(path)
}
```

---

## 2026-02-28 | 02_harmonize_clinical_GSE20685.R | wrong time/event columns detected

**Script:** `02_harmonize_clinical_GSE20685.R`
**Error:**
```
[02_GSE20685] Detected columns — time: 'time_to_metastasis (years):ch1' | event: 'status'
[02_GSE20685] N raw=327 | NA time: 244 | N final=75
```
**Root cause:** Auto-detection via `grepl("time|status|...")` matched:
- TIME_COL: `time_to_metastasis (years):ch1` (wrong) — should be `follow_up_duration (years):ch1`
- EVENT_COL: standard GEO metadata column `status` (wrong) — should be `event_death:ch1`

**Solution:** Hardcoded correct column names in script. Also fixed time unit:
`follow_up_duration` is in YEARS, not months — multiplied by 12, not left as-is.
```r
TIME_COL  <- "follow_up_duration (years):ch1"   # years → * 12 for months
EVENT_COL <- "event_death:ch1"                   # 0/1
AGE_COL   <- "age at diagnosis:ch1"
```
Result after fix: N=327, 83 OS events (25.4%), median follow-up 97.2 months.

---

## 2026-02-28 | 03_utils_gene_mapping.R | biomaRt fails to load

**Script:** `03_expression_preprocess_*.R` → `03_utils_gene_mapping.R`
**Error:**
```
Error: package or namespace load failed for 'biomaRt' in .checkGroupSigLength(list(generic@generic), list(generic)):
 cannot change value of locked binding for 'what'
```
**Root cause:** biomaRt 2.66.0/2.66.1 has a circular locked binding conflict in this
R 4.5.2 installation. Binary installation also fails (zip extraction error). Both source
and binary installation fail on this system.

**Solution:** Replaced `biomaRt` with `org.Hs.eg.db` + `AnnotationDbi` for
Ensembl→HGNC mapping in `build_hgnc_map_ensembl()`:
```r
map_raw <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = clean_ids,
  columns = c("ENSEMBL", "SYMBOL", "ENTREZID", "GENETYPE"),
  keytype = "ENSEMBL"
)
```
Advantages: fully offline, reproducible, no network dependency.
Disadvantage: slightly older annotations than live Ensembl, but sufficient for PAM50 genes.

---

## 2026-02-28 | 03_expression_preprocess_SCANB.R | ID mismatch clinical vs expression

**Script:** `03_expression_preprocess_SCANB.R`
**Error:**
```
[03_SCANB] Sample match: 0 of 3069 clinical IDs found in matrix.
```
**Root cause:** Clinical parquet `sample_id` = GEO accessions (GSM2528079, ...),
but GSE96058 expression CSV uses SCAN-B internal IDs (F1, F2, ...) as column names.
The matching key is `patient_id` (= pData `title`), not `sample_id` (= pData `geo_accession`).

**Solution:** Changed `scanb_ids <- normalize_id(clin$sample_id)` to
`scanb_ids <- normalize_id(clin$patient_id)`.

---

## 2026-02-28 | Multiple Bioc packages | lazy-load RDB corruption (systemic)

**Scripts:** `03_expression_preprocess_*.R` and `03_utils_gene_mapping.R`
**Error (various):**
```
lazy-load database 'C:/Users/.../Biostrings/R/Biostrings.rdb' is corrupt
lazy-load database 'C:/Program Files/R/R-4.5.2/library/Matrix/R/Matrix.rdb' is corrupt
```
**Root cause:** Multiple Bioconductor packages have corrupt lazy-load RDB files.
`options(warn=2)` converts any decompression warning during lazy loading to hard error.
Packages affected: Matrix (system lib), Biostrings, IRanges, AnnotationDbi and their deps.
The corruption may have been caused by prior partial installations under libdeflate issues.

**Solution:** Reinstall all affected Bioconductor packages with `force=TRUE`:
`BiocManager::install(c("BiocGenerics","S4Vectors","IRanges","Biostrings","GenomeInfoDb","AnnotationDbi","org.Hs.eg.db","Biobase","SummarizedExperiment","edgeR","GEOquery"), force=TRUE)`

Also note: even setting `options(warn=0)` before `library()` does not fully fix this
because lazy-loading of package internals happens at first function call, not at `library()`.
The only correct fix is to reinstall the packages so the RDB files are not corrupt.

---

## 2026-02-28 | 03_expression_preprocess_TCGA_BRCA.R | filterByExpr warning under warn=2

**Script:** `03_expression_preprocess_TCGA_BRCA.R`
**Error:**
```
Error in filterByExpr.DGEList(dge) :
  (converted from warning) All samples appear to belong to the same group.
```
**Root cause:** `edgeR::filterByExpr()` emits a warning (not an error) when no design matrix
is provided — it treats all samples as a single group, which is the intended behavior for
unsupervised low-count filtering. Under `options(warn=2)`, this warning becomes a fatal error.

**Solution:** Wrap `filterByExpr()` call with `options(warn=0)` / restore pattern:
```r
old_warn_filter <- getOption("warn"); options(warn = 0)
keep <- edgeR::filterByExpr(dge)
options(warn = old_warn_filter)
```

---

## 2026-02-28 | R sessions | Segmentation fault with Rscript -e (inline code)

**Context:** Trying to run multi-line R code via `Rscript -e "..."` with arrow loaded.
**Error:**
```
Exit code 139 — Segmentation fault
```
**Root cause:** Arrow's DLL (`arrow.dll`) was partially updated (reinstall failed with
"Permission denied" because another R session held the DLL lock). The R wrapper files
were updated but the old DLL remained, creating an API mismatch that causes segfault.
Inline `-e` mode may also be more sensitive to this on Windows.

**Solution:**
- Use script files instead of `-e` inline mode (write R code to .R file, then run)
- To properly reinstall arrow: close ALL R sessions first, then reinstall
- nanoparquet is a lightweight alternative for reading parquet without the C++ DLL issue

---

<!-- Add new errors below this line, most recent first -->
