# =============================================================================
# SCRIPT: 07z_table_figures_survival.R
# PURPOSE: Generate Table3_Survival_Performance_ByCohort as:
#          - CSV (results/main/)
#          - XLSX with formatting (results/main/)
#          - PNG figure for validation (figures/supp/)
# PROJECT: Core-PAM (Memorial v6.1)
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "07z_table_figures_survival.R"

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(openxlsx)
})

# --------------------------------------------------------------------------
# 1) Combine individual survival_results CSVs
# --------------------------------------------------------------------------
# Validation cohorts (Table 3 primary rows)
surv_files <- c(
  file.path(PATHS$results$supp, "survival_results_TCGA_BRCA.csv"),
  file.path(PATHS$results$supp, "survival_results_METABRIC.csv"),
  file.path(PATHS$results$supp, "survival_results_GSE20685.csv")
)
stopifnot(all(file.exists(surv_files)))

# Training cohort — C-index without optimism correction (no cross-validation)
scanb_path <- file.path(PATHS$results$supp, "survival_results_SCANB.csv")

tab_all <- bind_rows(lapply(surv_files, read_csv, show_col_types = FALSE))

# --------------------------------------------------------------------------
# 2) Table 3: primary survival performance
# --------------------------------------------------------------------------
format_row <- function(df) {
  df |>
    mutate(
      hr_uni_fmt    = sprintf("%.2f (%.2f–%.2f)", hr_uni, hr_uni_lo95, hr_uni_hi95),
      p_uni_fmt     = ifelse(p_uni < 0.001,
                             formatC(p_uni, format = "e", digits = 2),
                             sprintf("%.4f", p_uni)),
      c_index_fmt   = sprintf("%.3f (%.3f–%.3f)", c_index, c_index_lo95, c_index_hi95),
      hr_age_fmt    = ifelse(!is.na(hr_multi),
                             sprintf("%.2f (%.2f–%.2f)", hr_multi, hr_multi_lo95, hr_multi_hi95),
                             "—"),
      p_age_fmt     = ifelse(!is.na(p_multi),
                             ifelse(p_multi < 0.001,
                                    formatC(p_multi, format = "e", digits = 2),
                                    sprintf("%.4f", p_multi)),
                             "—")
    ) |>
    select(
      Cohort           = cohort,
      Role             = role,
      Endpoint         = endpoint,
      N                = n_samples,
      Events           = n_events,
      `FU med (mo)`    = fu_median_months,
      `HR uni (95%CI)` = hr_uni_fmt,
      `p (uni)`        = p_uni_fmt,
      `C_adj (95%CI)`  = c_index_fmt,
      `HR adj (CORE-A)`     = hr_age_fmt,
      `p (CORE-A adj)`    = p_age_fmt
    )
}

# Add Role column and build table: SCAN-B first (training), then validations
tab_scanb <- if (file.exists(scanb_path)) {
  read_csv(scanb_path, show_col_types = FALSE) |> mutate(role = "Training *")
} else NULL

tab_val <- tab_all |> mutate(role = "Validation")

tab3 <- bind_rows(tab_scanb, tab_val) |> format_row()

message("[", SCRIPT_NAME, "] Table 3 — ", nrow(tab3), " rows")
print(tab3)

# --------------------------------------------------------------------------
# 3) Save CSV (main results)
# --------------------------------------------------------------------------
main_dir <- PATHS$results$main
dir.create(main_dir, showWarnings = FALSE, recursive = TRUE)

csv_path <- file.path(main_dir, "Table3_Survival_Performance_ByCohort.csv")
write_csv(tab3, csv_path)
registry_append("ALL", "Table3_SurvivalPerformance_csv", csv_path,
                sha256_file(csv_path), "ok", SCRIPT_NAME, file.info(csv_path)$size / 1e6)
message("[", SCRIPT_NAME, "] CSV: ", csv_path)

# --------------------------------------------------------------------------
# 4) Save XLSX
# --------------------------------------------------------------------------
xlsx_path <- file.path(main_dir, "Table3_Survival_Performance_ByCohort.xlsx")
wb <- createWorkbook()
addWorksheet(wb, "Table3")
writeData(wb, "Table3", tab3)

# Header style
hs <- createStyle(fontColour = "#FFFFFF", fgFill = "#2C3E50",
                  halign = "CENTER", valign = "CENTER",
                  textDecoration = "Bold", border = "TopBottomLeftRight",
                  fontSize = 11)
addStyle(wb, "Table3", hs, rows = 1, cols = 1:ncol(tab3), gridExpand = TRUE)

# Row style: alternate shading
rs_even <- createStyle(fgFill = "#EBF5FB", halign = "CENTER",
                       border = "TopBottomLeftRight", fontSize = 10)
rs_odd  <- createStyle(fgFill = "#FDFEFE", halign = "CENTER",
                       border = "TopBottomLeftRight", fontSize = 10)
for (r in seq_len(nrow(tab3))) {
  sty <- if (r %% 2 == 0) rs_even else rs_odd
  addStyle(wb, "Table3", sty, rows = r + 1, cols = 1:ncol(tab3), gridExpand = TRUE)
}
setColWidths(wb, "Table3", cols = 1:ncol(tab3), widths = "auto")

# Footnote row: clarify CORE-A covariates per cohort
note_row   <- nrow(tab3) + 3
note_style <- createStyle(fontSize = 9, fontColour = "#5D6D7E",
                          textDecoration = "italic")
note_text  <- paste0(
  "* Training cohort (SCAN-B): C_adj = 0.698 (in-sample; no optimism correction). ",
  "OOF C-index used for model selection = 0.670 (PAM50-full OOF = 0.679; gap = 0.009 < \u03b4 = 0.010, non-inferiority confirmed). ",
  "** C_adj = max(C_raw, 1 − C_raw); direction-invariant; 95% CI by bootstrap (B=1,000). ",
  "*** HR adj (CORE-A): Cox model adjusted for available covariates. ",
  "CORE-A covariates by cohort: SCAN-B (age + ER status), METABRIC (age + ER status), ",
  "TCGA-BRCA (age only), GSE20685 (age only). Covariates included only if >80% non-missing."
)
writeData(wb, "Table3", note_text, startRow = note_row, startCol = 1)
addStyle(wb, "Table3", note_style, rows = note_row, cols = 1, gridExpand = FALSE)

saveWorkbook(wb, xlsx_path, overwrite = TRUE)
registry_append("ALL", "Table3_SurvivalPerformance_xlsx", xlsx_path,
                sha256_file(xlsx_path), "ok", SCRIPT_NAME, file.info(xlsx_path)$size / 1e6)
message("[", SCRIPT_NAME, "] XLSX: ", xlsx_path)

# --------------------------------------------------------------------------
# 5) PNG figure of the table (for validation)
# --------------------------------------------------------------------------
# Use base R table rendering via png + grid
# Table treated as figure → figures/supp/en/png/ (EN only; no bilingual table)
tbl_png <- file.path(PATHS$figures$supp_en_png, "FigT_Table3_Survival_Performance.png")

tryCatch({
  library(gridExtra)
  library(grid)
  grob <- tableGrob(
    tab3,
    rows = NULL,
    theme = ttheme_default(
      core    = list(fg_params = list(fontsize = 8),
                     bg_params = list(fill = c("#FDFEFE", "#EBF5FB"))),
      colhead = list(fg_params = list(fontsize = 9, fontface = "bold",
                                     col = "white"),
                     bg_params = list(fill = "#2C3E50"))
    )
  )
  png(tbl_png, width = 14, height = 2.5 + nrow(tab3) * 0.5,
      units = "in", res = 200)
  grid.newpage()
  grid.draw(grob)
  dev.off()
  registry_append("ALL", "FigT_Table3_png", tbl_png,
                  sha256_file(tbl_png), "ok", SCRIPT_NAME, file.info(tbl_png)$size / 1e6)
  message("[", SCRIPT_NAME, "] PNG figure: ", tbl_png)
}, error = function(e) {
  message("[", SCRIPT_NAME, "] PNG figure skipped (gridExtra): ", e$message)
})

message("[", SCRIPT_NAME, "] COMPLETED")
