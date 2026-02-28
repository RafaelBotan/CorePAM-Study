# =============================================================================
# SCRIPT: 03_expression_preprocess_METABRIC.R
# PURPOSE: Pré-processamento de expressão METABRIC.
#          Microarray Illumina HT-12 v3: log2 intensity as-is → HGNC map.
# PROJETO: Core-PAM (Memorial v6.1 §4.3)
#
# INPUT:
#   01_Base_Pura_CorePAM/RAW/METABRIC/brca_metabric/data_mrna_illumina_microarray.txt
#   01_Base_Pura_CorePAM/PROCESSED/METABRIC/clinical_FINAL.parquet
#
# OUTPUTS:
#   01_Base_Pura_CorePAM/PROCESSED/METABRIC/expression_genelevel_preZ.parquet
#   results/supp/gene_mapping_audit_METABRIC.csv
#
# REGRAS (Memorial v6.1 §4.3):
#   - Microarray: log2-intensity "as-is" (sem TMM).
#   - Padronização Z-score intra-coorte em 06_zscore_and_score_METABRIC.R.
#   - Múltiplas probes por gene → maior variância intra-coorte.
#   - HGNC mapping: arquivo cBioPortal já tem Hugo_Symbol → verificar via HGNChelper.
# =============================================================================

source("scripts/00_setup.R")
source("scripts/03_utils_gene_mapping.R")

SCRIPT_NAME <- "03_expression_preprocess_METABRIC.R"
COHORT      <- "METABRIC"

# =============================================================================
# 1) LEITURA DA MATRIZ (cBioPortal format)
#    Formato: Hugo_Symbol | Entrez_Gene_Id | SAMPLE1 | SAMPLE2 | ...
#    Linhas com "#" são metadados → pular
# =============================================================================
expr_path <- file.path(raw_cohort(COHORT), "brca_metabric",
                       "data_mrna_illumina_microarray.txt")
if (!file.exists(expr_path)) {
  stop("[03_METABRIC] Arquivo não encontrado: ", expr_path,
       "\nExecute 01_download_raw_data.R primeiro.")
}

message("[03_METABRIC] Lendo matriz de expressão cBioPortal...")
lines  <- readLines(expr_path, warn = FALSE)
skip_n <- sum(startsWith(lines, "#"))

old_warn <- getOption("warn"); options(warn = 0)
expr_raw <- read_tsv(expr_path, skip = skip_n, show_col_types = FALSE,
                     name_repair = "unique")
options(warn = old_warn)

message(sprintf("[03_METABRIC] Lido: %d linhas x %d colunas (incluindo Hugo_Symbol, Entrez)",
                nrow(expr_raw), ncol(expr_raw)))

# Separar colunas de genes e amostras
GENE_COL   <- "Hugo_Symbol"
ENTREZ_COL <- "Entrez_Gene_Id"
sample_cols <- setdiff(names(expr_raw), c(GENE_COL, ENTREZ_COL))

gene_symbols_raw <- expr_raw[[GENE_COL]]
entrez_ids_raw   <- expr_raw[[ENTREZ_COL]]

expr_mat <- as.matrix(expr_raw[, sample_cols])
rownames(expr_mat) <- gene_symbols_raw
storage.mode(expr_mat) <- "numeric"

message(sprintf("[03_METABRIC] Matriz: %d genes x %d amostras",
                nrow(expr_mat), ncol(expr_mat)))

# =============================================================================
# 2) FILTRAR AMOSTRAS COM CLÍNICA
# =============================================================================
clin_path <- file.path(proc_cohort(COHORT), "clinical_FINAL.parquet")
clin      <- strict_parquet(clin_path)
clin_ids  <- normalize_id(clin$sample_id)

col_ids_norm <- normalize_id(colnames(expr_mat))
keep_cols    <- col_ids_norm %in% clin_ids
n_matched    <- sum(keep_cols)

message(sprintf("[03_METABRIC] Match amostras: %d de %d IDs clínicos encontrados.",
                n_matched, length(clin_ids)))

if (n_matched == 0) {
  stop("[03_METABRIC] Nenhuma amostra clínica encontrada na matriz.\n",
       "Verifique IDs em clinical_FINAL.parquet vs sample columns.")
}

expr_mat <- expr_mat[, keep_cols, drop = FALSE]
colnames(expr_mat) <- col_ids_norm[keep_cols]

message(sprintf("[03_METABRIC] Matriz filtrada: %d genes x %d amostras",
                nrow(expr_mat), ncol(expr_mat)))

# =============================================================================
# 3) VERIFICAR ESCALA LOG2
#    Illumina HT-12 cBioPortal: tipicamente já em log2.
#    Heurística: se 95° percentil < 20 → assume log2; senão → aplicar log2.
# =============================================================================
p95 <- quantile(expr_mat, 0.95, na.rm = TRUE)
if (p95 > 20) {
  message(sprintf("[03_METABRIC] p95=%.1f > 20 → aplicando log2(x + 1).", p95))
  expr_mat <- log2(expr_mat + 1)
} else {
  message(sprintf("[03_METABRIC] p95=%.1f ≤ 20 → dados já em log2 (as-is).", p95))
}

# =============================================================================
# 4) MAPEAMENTO HGNC (Hugo_Symbol já disponível — validar via HGNChelper)
# =============================================================================
message("[03_METABRIC] Validando Hugo Symbols via HGNChelper...")
hgnc_res     <- build_hgnc_map_symbol(rownames(expr_mat))
gene_symbols <- hgnc_res$hgnc_symbol

# Remover genes sem mapeamento HGNC válido
keep_genes   <- !is.na(gene_symbols) & gene_symbols != ""
n_unmapped   <- sum(!keep_genes)
expr_mat     <- expr_mat[keep_genes, , drop = FALSE]
gene_symbols <- gene_symbols[keep_genes]

message(sprintf("[03_METABRIC] HGNC: %d válidos | %d descartados",
                sum(keep_genes), n_unmapped))

# =============================================================================
# 5) COLAPSO DE PROBES (maior variância intra-coorte)
# =============================================================================
n_before <- nrow(expr_mat)
expr_mat  <- collapse_probes_by_variance(expr_mat, gene_symbols)
n_after   <- nrow(expr_mat)
message(sprintf("[03_METABRIC] Colapso probes: %d → %d genes únicos HGNC",
                n_before, n_after))

# =============================================================================
# 6) QC: remover genes com missingness > 20%
# =============================================================================
miss_frac <- rowMeans(is.na(expr_mat))
n_high_miss <- sum(miss_frac > 0.20)
if (n_high_miss > 0) {
  message(sprintf("[03_METABRIC] Removendo %d genes com >20%% missingness.",
                  n_high_miss))
  expr_mat <- expr_mat[miss_frac <= 0.20, , drop = FALSE]
}

# =============================================================================
# 7) AUDIT PAM50
# =============================================================================
pam50_genes <- c(
  "ACTR3B","ANLN","BAG1","BCL2","BIRC5","BLVRA","CCNB1","CCNE1","CDC20",
  "CDC6","CDH3","CENPF","CEP55","CXXC5","EGFR","ERBB2","ESR1","EXO1",
  "FGFR4","FOXA1","FOXC1","GPR160","GRB7","KIF2C","KRT14","KRT17","KRT5",
  "MAPT","MDM2","MELK","MIA","MKI67","MLPH","MMP11","MYBL2","MYC","NAT1",
  "NDC80","NUF2","ORC6","PGR","PHGDH","PTTG1","RRM2","SFRP1","SLC39A6",
  "TMEM45B","TYMS","UBE2C","UBE2T"
)
pam50_present <- intersect(pam50_genes, rownames(expr_mat))
pam50_missing <- setdiff(pam50_genes, rownames(expr_mat))
message(sprintf("[03_METABRIC] PAM50: %d/50 presentes | ausentes: %s",
                length(pam50_present),
                if (length(pam50_missing) == 0) "nenhum"
                else paste(pam50_missing, collapse = ",")))

audit_df <- tibble(
  cohort   = COHORT,
  gene     = rownames(expr_mat),
  present  = TRUE,
  in_pam50 = rownames(expr_mat) %in% pam50_genes
) |>
  bind_rows(tibble(cohort = COHORT, gene = pam50_missing,
                   present = FALSE, in_pam50 = TRUE))
save_gene_mapping_audit(audit_df, COHORT, SCRIPT_NAME)

# =============================================================================
# 8) SALVAR
# =============================================================================
dest_dir <- proc_cohort(COHORT)
dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(dest_dir, "expression_genelevel_preZ.parquet")

expr_df  <- tibble::rownames_to_column(as.data.frame(expr_mat), var = "gene")
arrow::write_parquet(expr_df, out_path)
h    <- sha256_file(out_path)
size <- file.info(out_path)$size / 1024^2
registry_append(COHORT, "Expression_preZ", out_path, h,
                "INTEGRO", SCRIPT_NAME, size)
message(sprintf("[03_METABRIC] Salvo: %d genes x %d amostras | %.2f MB | SHA256: %s",
                nrow(expr_mat), ncol(expr_mat), size, h))

message("\n[03_METABRIC] Concluido: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
message("Proximo passo: scripts/03_expression_preprocess_GSE20685.R")
