#!/usr/bin/env Rscript
# Render thesis QMD to DOCX
setwd("Y:/Phd-Genomic-claude/tese")
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
