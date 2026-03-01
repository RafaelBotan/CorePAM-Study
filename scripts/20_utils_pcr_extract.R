# =============================================================================
# UTILITY: 20_utils_pcr_extract.R
# PURPOSE: Shared helpers for extracting pCR (pathologic complete response),
#          age, and ER/HER2 status from GEO pData objects across all pCR cohorts.
#
# Detection strategy (extract_pcr_column):
#   Phase 1 — Column-name priority: look for columns whose NAME contains
#             "pcr", "patholog.*response", or similar. Parse 0/1 or yes/no values.
#   Phase 2 — Value search with coverage filter: scan all columns for values
#             matching pCR/RD patterns; require >= MIN_COVERAGE_FRAC valid calls.
#
# PROJETO: Core-PAM (Memorial v6.1 / Freeze Core-PAM — pCR block)
# =============================================================================

MIN_COVERAGE_FRAC <- 0.10   # minimum fraction of samples with valid pCR to accept column

# --------------------------------------------------------------------------
# Internal: extract value after colon from "key: value" GEO cells
# --------------------------------------------------------------------------
.geo_extract_value <- function(x) {
  trimws(sub("^[^:]+:", "", as.character(x)))
}

# --------------------------------------------------------------------------
# Internal: try to parse 0/1 or yes/no numeric pCR column
# --------------------------------------------------------------------------
.try_parse_01 <- function(vals_raw, n_samples, cohort, col) {
  val_part <- .geo_extract_value(vals_raw)
  lv       <- trimws(tolower(val_part))

  # Numeric 0 / 1 encoding (most common in I-SPY1-style datasets)
  num <- suppressWarnings(as.integer(lv))
  if (all(num %in% c(0L, 1L, NA_integer_))) {
    n_valid <- sum(!is.na(num))
    if (n_valid >= max(10, MIN_COVERAGE_FRAC * n_samples)) {
      message(sprintf(
        "  [pcr_extract:%s] Found pCR (0/1 encoding) in column '%s' — pCR+=%d, pCR-=%d, NA=%d",
        cohort, col, sum(num == 1L, na.rm = TRUE),
        sum(num == 0L, na.rm = TRUE), sum(is.na(num))
      ))
      return(list(pcr = num, col_used = col))
    }
  }

  # yes / no encoding
  has_yes <- any(grepl("^yes$|^pcr$|^pCR$", lv), na.rm = TRUE)
  has_no  <- any(grepl("^no$|^rd$|^RD$|^non-pcr$", lv, ignore.case = TRUE), na.rm = TRUE)
  if (has_yes && has_no) {
    pcr_vec <- dplyr::case_when(
      grepl("^yes$|^pcr$",  lv, ignore.case = TRUE) ~ 1L,
      grepl("^no$|^rd$|^non-pcr$", lv, ignore.case = TRUE) ~ 0L,
      TRUE ~ NA_integer_
    )
    n_valid <- sum(!is.na(pcr_vec))
    if (n_valid >= max(10, MIN_COVERAGE_FRAC * n_samples)) {
      message(sprintf(
        "  [pcr_extract:%s] Found pCR (yes/no) in column '%s' — pCR+=%d, pCR-=%d, NA=%d",
        cohort, col, sum(pcr_vec == 1L, na.rm = TRUE),
        sum(pcr_vec == 0L, na.rm = TRUE), sum(is.na(pcr_vec))
      ))
      return(list(pcr = pcr_vec, col_used = col))
    }
  }

  NULL
}

# --------------------------------------------------------------------------
# extract_pcr_column
#   Searches pData for a column containing pCR information.
#   Returns: list(pcr = integer 0/1/NA vector, col_used = character)
# --------------------------------------------------------------------------
extract_pcr_column <- function(pdata, cohort = "?") {

  n_samples <- nrow(pdata)

  # Patterns for column names that directly encode pCR
  name_re_pcr <- "pcr|patholog.*response|pathologic.*response|pCR|pcr_vs_rd|pcr.ncr|pcr.rd|response.*pcr"

  # pCR-positive value patterns (case-insensitive on extracted value)
  pos_re <- paste0(
    "\\bpcr\\b|patholog.*complete|complete.*response|",
    "pCR$|^pCR|\\byes\\b"
  )
  # pCR-negative value patterns — includes nCR (near complete response = NOT pCR)
  neg_re <- paste0(
    "\\brd\\b|residual.disease|\\bno\\b|not.pcr|non.pcr|non-pcr|",
    "\\bncr\\b|near.complete"
  )

  # ------------------------------------------------------------------
  # Phase 1: Column NAME priority
  # ------------------------------------------------------------------
  name_match_cols <- grep(name_re_pcr, names(pdata), ignore.case = TRUE, value = TRUE)

  for (col in name_match_cols) {
    vals_raw <- as.character(pdata[[col]])

    # Try 0/1 or yes/no parsing first (common in I-SPY1 style data)
    result <- .try_parse_01(vals_raw, n_samples, cohort, col)
    if (!is.null(result)) return(result)

    # Try pCR / RD / nCR pattern matching
    val_part <- .geo_extract_value(vals_raw)
    lv <- tolower(val_part)
    has_pos <- any(grepl(pos_re, lv, ignore.case = TRUE), na.rm = TRUE)
    has_neg <- any(grepl(neg_re, lv, ignore.case = TRUE), na.rm = TRUE)

    if (has_pos && has_neg) {
      pcr_vec <- dplyr::case_when(
        grepl(pos_re, lv, ignore.case = TRUE) ~ 1L,
        grepl(neg_re, lv, ignore.case = TRUE) ~ 0L,
        TRUE ~ NA_integer_
      )
      n_valid <- sum(!is.na(pcr_vec))
      if (n_valid >= max(10, MIN_COVERAGE_FRAC * n_samples)) {
        message(sprintf(
          "  [pcr_extract:%s] Found pCR (name-match) in column '%s' — pCR+=%d, pCR-=%d, NA=%d",
          cohort, col,
          sum(pcr_vec == 1L, na.rm = TRUE), sum(pcr_vec == 0L, na.rm = TRUE),
          sum(is.na(pcr_vec))
        ))
        return(list(pcr = pcr_vec, col_used = col))
      }
    }
  }

  # ------------------------------------------------------------------
  # Phase 2: Value search across all columns (with coverage filter)
  # ------------------------------------------------------------------
  for (col in names(pdata)) {
    vals_raw <- as.character(pdata[[col]])
    val_part <- .geo_extract_value(vals_raw)
    lv <- tolower(val_part)

    has_pos <- any(grepl(pos_re, lv, ignore.case = TRUE), na.rm = TRUE)
    has_neg <- any(grepl(neg_re, lv, ignore.case = TRUE), na.rm = TRUE)

    if (has_pos && has_neg) {
      pcr_vec <- dplyr::case_when(
        grepl(pos_re, lv, ignore.case = TRUE) ~ 1L,
        grepl(neg_re, lv, ignore.case = TRUE) ~ 0L,
        TRUE ~ NA_integer_
      )
      n_valid <- sum(!is.na(pcr_vec))
      # Coverage filter: skip if too few valid calls
      if (n_valid < max(10, MIN_COVERAGE_FRAC * n_samples)) {
        message(sprintf(
          "  [pcr_extract:%s] Column '%s' matched patterns but coverage too low (%d/%d) — skipping",
          cohort, col, n_valid, n_samples
        ))
        next
      }
      message(sprintf(
        "  [pcr_extract:%s] Found pCR in column '%s' — pCR+=%d, pCR-=%d, NA=%d",
        cohort, col,
        sum(pcr_vec == 1L, na.rm = TRUE), sum(pcr_vec == 0L, na.rm = TRUE),
        sum(is.na(pcr_vec))
      ))
      return(list(pcr = pcr_vec, col_used = col))
    }
  }

  # ------------------------------------------------------------------
  # Fallback: print all columns + unique values for debugging
  # ------------------------------------------------------------------
  message(sprintf("\n  [pcr_extract:%s] *** pCR column NOT found. Listing all pData columns: ***",
                  cohort))
  for (col in names(pdata)) {
    uvals <- unique(as.character(pdata[[col]]))
    if (length(uvals) <= 20) {
      message(sprintf("  %s: %s", col, paste(uvals, collapse = " | ")))
    } else {
      message(sprintf("  %s: [%d unique values] %s ...",
                      col, length(uvals), paste(head(uvals, 5), collapse = " | ")))
    }
  }
  stop(sprintf("[pcr_extract:%s] Could not auto-detect pCR column. Check pData printout above.",
               cohort))
}

# --------------------------------------------------------------------------
# extract_numeric_covariate
#   Extracts a numeric covariate (e.g., age) from pData characteristics.
# --------------------------------------------------------------------------
extract_numeric_covariate <- function(pdata, pattern = "age") {
  for (col in grep(pattern, names(pdata), ignore.case = TRUE, value = TRUE)) {
    vals <- as.character(pdata[[col]])
    val_part <- .geo_extract_value(vals)
    num <- suppressWarnings(as.numeric(val_part))
    if (mean(!is.na(num)) >= 0.5) return(num)
  }
  chars_cols <- grep("^characteristics_ch", names(pdata), value = TRUE)
  for (col in chars_cols) {
    vals     <- as.character(pdata[[col]])
    key_part <- trimws(sub(":.*$", "", vals))
    if (any(grepl(pattern, key_part, ignore.case = TRUE))) {
      val_part <- .geo_extract_value(vals)
      num <- suppressWarnings(as.numeric(val_part))
      if (mean(!is.na(num)) >= 0.5) return(num)
    }
  }
  return(rep(NA_real_, nrow(pdata)))
}

# --------------------------------------------------------------------------
# extract_binary_covariate
#   Extracts a binary covariate (e.g., ER status) from pData.
# --------------------------------------------------------------------------
extract_binary_covariate <- function(pdata, pattern = "er.status",
                                     pos_vals = c("positive", "pos", "er+", "1", "yes"),
                                     neg_vals = c("negative", "neg", "er-", "0", "no")) {
  pos_re <- paste(paste0("\\b", pos_vals, "\\b"), collapse = "|")
  neg_re <- paste(paste0("\\b", neg_vals, "\\b"), collapse = "|")

  for (col in grep(pattern, names(pdata), ignore.case = TRUE, value = TRUE)) {
    vals    <- tolower(.geo_extract_value(as.character(pdata[[col]])))
    has_pos <- any(grepl(pos_re, vals, ignore.case = TRUE), na.rm = TRUE)
    has_neg <- any(grepl(neg_re, vals, ignore.case = TRUE), na.rm = TRUE)
    if (has_pos || has_neg) {
      return(dplyr::case_when(
        grepl(pos_re, vals, ignore.case = TRUE) ~ "positive",
        grepl(neg_re, vals, ignore.case = TRUE) ~ "negative",
        TRUE ~ NA_character_
      ))
    }
  }
  chars_cols <- grep("^characteristics_ch", names(pdata), value = TRUE)
  for (col in chars_cols) {
    vals     <- as.character(pdata[[col]])
    key_part <- trimws(sub(":.*$", "", vals))
    if (any(grepl(pattern, key_part, ignore.case = TRUE))) {
      val_part <- tolower(.geo_extract_value(vals))
      return(dplyr::case_when(
        grepl(pos_re, val_part, ignore.case = TRUE) ~ "positive",
        grepl(neg_re, val_part, ignore.case = TRUE) ~ "negative",
        TRUE ~ NA_character_
      ))
    }
  }
  return(rep(NA_character_, nrow(pdata)))
}

message("[20_utils_pcr_extract] pCR extraction utilities loaded.")
