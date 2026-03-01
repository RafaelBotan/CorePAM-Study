library(arrow)
clin <- read_parquet("01_Base_Pura_CorePAM/PROCESSED/METABRIC/clinical_FINAL.parquet")
cat("Columns:", paste(names(clin), collapse=", "), "\n")
if ("er_status" %in% names(clin)) {
  cat("er_status values:\n")
  print(table(clin$er_status, useNA="always"))
} else {
  cat("er_status NOT found in clinical_FINAL.parquet\n")
}
# Also check ER_IHC raw file unique values
raw <- readr::read_tsv(
  "01_Base_Pura_CorePAM/RAW/METABRIC/brca_metabric/data_clinical_patient.txt",
  skip = 4, show_col_types = FALSE
)
cat("\nER_IHC unique values (raw cBioPortal):\n")
print(table(raw$ER_IHC, useNA="always"))
