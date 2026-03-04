#!/usr/bin/env Rscript
# Render thesis QMD to DOCX
# Run from repo root: Rscript tese/render_tese.R
# Or from tese/ dir:  Rscript render_tese.R

if (!file.exists("tese_corepam_v3.qmd")) {
  if (file.exists("tese/tese_corepam_v3.qmd")) {
    setwd("tese")
  } else {
    stop("Cannot find tese_corepam_v3.qmd. Run from repo root or tese/ directory.")
  }
}
cat("Working dir:", getwd(), "\n")
cat("Rendering tese_corepam_v3.qmd to DOCX...\n")
system2(
  "quarto",
  args = c("render", "tese_corepam_v3.qmd", "--to", "docx"),
  stdout = "", stderr = ""
)
if (file.exists("tese_corepam_v3.docx")) {
  cat("SUCCESS: tese_corepam_v3.docx created\n")
  cat("Size:", round(file.size("tese_corepam_v3.docx") / 1e6, 2), "MB\n")
} else {
  cat("FAILED: tese_corepam_v3.docx not created\n")
}
