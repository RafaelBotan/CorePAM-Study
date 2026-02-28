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

<!-- Add new errors below this line, most recent first -->
