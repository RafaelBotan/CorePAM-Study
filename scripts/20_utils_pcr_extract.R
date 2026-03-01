# =============================================================================
# UTILITY: 20_utils_pcr_extract.R
# PURPOSE: Shared helpers for extracting pCR (pathologic complete response),
#          age, and ER status from GEO pData objects across all pCR-block cohorts.
#
# GEO pData characteristics columns follow the format "key: value" (e.g.,
# "pathological response: pCR") in columns named characteristics_ch1,
# characteristics_ch1.1, etc. These helpers search all columns dynamically.
#
# PROJETO: Core-PAM (Memorial v6.1 / Freeze Core-PAM — pCR block)
# =============================================================================

# --------------------------------------------------------------------------
# Internal: extract value after colon from "key: value" GEO cells
# --------------------------------------------------------------------------
.geo_extract_value <- function(x) {
  # x is a character vector like "pathological response: pCR"
  # Returns the part after the last colon (trimmed)
  trimws(sub("^[^:]+:", "", as.character(x)))
}

# --------------------------------------------------------------------------
# extract_pcr_column
#   Searches pData for a column containing pCR information.
#   Returns: list(pcr = integer 0/1/NA vector, col_used = character)
#
#   Detection order:
#   1. Column name matching "pcr|patholog|response|pcr_status"
#   2. Column VALUE matching pCR positive patterns in any characteristics col
#   3. Stop with informative error listing all columns + unique values
# --------------------------------------------------------------------------
extract_pcr_column <- function(pdata, cohort = "?") {

  chars_cols <- grep("^characteristics_ch1|^pcr|^response|^patholog",
                     names(pdata), ignore.case = TRUE, value = TRUE)

  # pCR-positive patterns (case-insensitive)
  pos_re <- "\\bpcr\\b|patholog.*complete|complete.*response|pCR$|^pCR|\\byes\\b"
  # pCR-negative patterns (case-insensitive)
  neg_re <- "\\brd\\b|residual.disease|\\bno\\b|not.pcr|non.pcr"

  for (col in names(pdata)) {
    vals <- as.character(pdata[[col]])
    # Extract value part if "key: value" format
    val_part <- .geo_extract_value(vals)
    lv <- tolower(val_part)

    has_pos <- any(grepl(pos_re, lv, ignore.case = TRUE), na.rm = TRUE)
    has_neg <- any(grepl(neg_re, lv, ignore.case = TRUE), na.rm = TRUE)

    if (has_pos && has_neg) {
      pcr_vec <- dplyr::case_when(
        grepl(pos_re, lv, ignore.case = TRUE) ~ 1L,
        grepl(neg_re, lv, ignore.case = TRUE) ~ 0L,
        TRUE ~ NA_integer_
      )
      message(sprintf(
        "  [pcr_extract:%s] Found pCR in column '%s' — pCR+=%d, pCR-=%d, NA=%d",
        cohort, col,
        sum(pcr_vec == 1L, na.rm = TRUE),
        sum(pcr_vec == 0L, na.rm = TRUE),
        sum(is.na(pcr_vec))
      ))
      return(list(pcr = pcr_vec, col_used = col))
    }
  }

  # Fallback: print all columns and their unique values for debugging
  message(sprintf("\n  [pcr_extract:%s] *** pCR column NOT found. Listing all pData columns: ***", cohort))
  for (col in names(pdata)) {
    uvals <- unique(as.character(pdata[[col]]))
    if (length(uvals) <= 20) {
      message(sprintf("  %s: %s", col, paste(uvals, collapse = " | ")))
    } else {
      message(sprintf("  %s: [%d unique values] %s ...",
                      col, length(uvals), paste(head(uvals, 5), collapse = " | ")))
    }
  }
  stop(sprintf(
    "[pcr_extract:%s] Could not auto-detect pCR column. Check pData printout above.\n",
    "Set the correct column in the prepare script with: pdata$pcr <- ...",
    cohort
  ))
}

# --------------------------------------------------------------------------
# extract_numeric_covariate
#   Extracts a numeric covariate (e.g., age) from pData characteristics.
#   pattern: regex applied to column names OR to the "key" part of values.
# --------------------------------------------------------------------------
extract_numeric_covariate <- function(pdata, pattern = "age") {
  for (col in grep(pattern, names(pdata), ignore.case = TRUE, value = TRUE)) {
    vals <- as.character(pdata[[col]])
    val_part <- .geo_extract_value(vals)
    num <- suppressWarnings(as.numeric(val_part))
    if (mean(!is.na(num)) >= 0.5) return(num)
  }
  # Try extracting value part from any characteristics column whose "key" matches
  chars_cols <- grep("^characteristics_ch1", names(pdata), value = TRUE)
  for (col in chars_cols) {
    vals <- as.character(pdata[[col]])
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
#   Returns character "positive"/"negative"/NA.
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
  # Try characteristics columns
  chars_cols <- grep("^characteristics_ch1", names(pdata), value = TRUE)
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
