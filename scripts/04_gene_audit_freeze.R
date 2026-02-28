# =============================================================================
# SCRIPT: 04_gene_audit_freeze.R
# PURPOSE: Auditoria transversal de genes PAM50 em todas as coortes.
#          Verifica cobertura do painel em treino e validações.
# PROJETO: Core-PAM (Memorial v6.1 §4.5 / Checklist §7)
#
# INPUTS:
#   01_Base_Pura_CorePAM/PROCESSED/<COHORT>/expression_genelevel_preZ.parquet
#
# OUTPUTS:
#   results/supp/gene_audit_by_cohort.csv   (auditoria principal)
#   results/supp/pam50_coverage_summary.csv (cobertura por coorte)
#
# CRITÉRIO GO/NO-GO (Memorial §0 Checklist):
#   SCANB (treino): todos os 50 genes PAM50 devem estar presentes.
#   Validações:     genes_present >= 80% do Core-PAM final.
#                   (neste script: auditamos PAM50; cobertura do Core-PAM
#                    é re-auditada em 06_zscore_and_score_<COHORT>.R)
# =============================================================================

source("scripts/00_setup.R")
SCRIPT_NAME <- "04_gene_audit_freeze.R"

# =============================================================================
# 1) LISTA CANÔNICA PAM50 (50 genes)
# =============================================================================
PAM50_GENES <- c(
  "ACTR3B","ANLN","BAG1","BCL2","BIRC5","BLVRA","CCNB1","CCNE1","CDC20",
  "CDC6","CDH3","CENPF","CEP55","CXXC5","EGFR","ERBB2","ESR1","EXO1",
  "FGFR4","FOXA1","FOXC1","GPR160","GRB7","KIF2C","KRT14","KRT17","KRT5",
  "MAPT","MDM2","MELK","MIA","MKI67","MLPH","MMP11","MYBL2","MYC","NAT1",
  "NDC80","NUF2","ORC6","PGR","PHGDH","PTTG1","RRM2","SFRP1","SLC39A6",
  "TMEM45B","TYMS","UBE2C","UBE2T"
)
stopifnot(length(PAM50_GENES) == 50)
message(sprintf("[04_audit] PAM50 canônico: %d genes", length(PAM50_GENES)))

# =============================================================================
# 2) COORTES A AUDITAR
# =============================================================================
COHORTS <- list(
  list(id = "SCANB",     role = "TRAIN",      platform = "RNAseq"),
  list(id = "GSE96058",  role = "VALIDATION", platform = "RNAseq"),
  list(id = "TCGA_BRCA", role = "VALIDATION", platform = "RNAseq"),
  list(id = "METABRIC",  role = "VALIDATION", platform = "Microarray"),
  list(id = "GSE20685",  role = "VALIDATION", platform = "Microarray")
)

# =============================================================================
# 3) AUDITORIA POR COORTE
# =============================================================================
audit_rows <- list()
coverage_rows <- list()

for (coh in COHORTS) {
  cohort   <- coh$id
  role     <- coh$role
  platform <- coh$platform

  expr_path <- file.path(proc_cohort(cohort), "expression_genelevel_preZ.parquet")

  if (!file.exists(expr_path)) {
    message(sprintf("[04_audit] AVISO: %s — arquivo nao encontrado: %s",
                    cohort, expr_path))
    message(sprintf("  Execute 03_expression_preprocess_%s.R primeiro.", cohort))

    # Registrar como ausente mesmo assim
    coverage_rows[[cohort]] <- tibble(
      cohort   = cohort,
      role     = role,
      platform = platform,
      n_pam50_present  = NA_integer_,
      n_pam50_missing  = NA_integer_,
      pct_pam50        = NA_real_,
      genes_missing    = NA_character_,
      go_nogo          = "NO-GO (arquivo ausente)",
      file_exists      = FALSE
    )
    next
  }

  # Leitura estrita
  expr <- strict_parquet(expr_path)
  genes_in_cohort <- expr$gene   # coluna "gene" adicionada em 03_*.R

  # Cobertura PAM50
  present <- intersect(PAM50_GENES, genes_in_cohort)
  missing <- setdiff(PAM50_GENES, genes_in_cohort)
  pct     <- 100 * length(present) / length(PAM50_GENES)

  message(sprintf("[04_audit] %-12s | %d genes totais | PAM50: %d/50 (%.1f%%) | ausentes: %s",
                  cohort, length(genes_in_cohort), length(present), pct,
                  if (length(missing) == 0) "nenhum" else paste(missing, collapse = ",")))

  # GO/NO-GO para treino: exigir 50/50
  if (role == "TRAIN") {
    go_nogo <- if (length(missing) == 0) "GO" else
      sprintf("NO-GO (faltam %d genes PAM50 no treino: %s)",
              length(missing), paste(missing, collapse = ","))
  } else {
    # Validação: GO se >=80% do PAM50 (proxy; será re-verificado para Core-PAM)
    go_nogo <- if (pct >= 80) "GO (>=80% PAM50)" else
      sprintf("AVISO (<80%% PAM50: %.1f%%)", pct)
  }

  # Linha de cobertura
  coverage_rows[[cohort]] <- tibble(
    cohort          = cohort,
    role            = role,
    platform        = platform,
    n_genes_total   = length(genes_in_cohort),
    n_pam50_present = length(present),
    n_pam50_missing = length(missing),
    pct_pam50       = round(pct, 2),
    genes_missing   = if (length(missing) == 0) "" else paste(missing, collapse = ";"),
    go_nogo         = go_nogo,
    file_exists     = TRUE
  )

  # Linha por gene (para tabela gene × coorte)
  for (g in PAM50_GENES) {
    audit_rows[[paste0(cohort, "_", g)]] <- tibble(
      cohort   = cohort,
      role     = role,
      platform = platform,
      gene     = g,
      present  = g %in% genes_in_cohort
    )
  }
}

# =============================================================================
# 4) TABELA GENE × COORTE (wide)
# =============================================================================
df_audit <- bind_rows(audit_rows)

df_wide <- df_audit |>
  select(gene, cohort, present) |>
  tidyr::pivot_wider(names_from = cohort, values_from = present,
                     values_fill = FALSE) |>
  mutate(n_cohorts_present = rowSums(across(where(is.logical)))) |>
  arrange(desc(n_cohorts_present), gene)

message("\n[04_audit] Tabela gene × coorte (PAM50):")
print(df_wide)

# =============================================================================
# 5) TABELA DE COBERTURA POR COORTE
# =============================================================================
df_coverage <- bind_rows(coverage_rows)

message("\n[04_audit] Resumo de cobertura:")
print(df_coverage |> select(cohort, role, n_pam50_present, pct_pam50, go_nogo))

# =============================================================================
# 6) SALVAR E REGISTRAR
# =============================================================================
# Tabela gene × coorte
out_audit <- file.path(PATHS$results$supp, "gene_audit_by_cohort.csv")
write_csv(df_wide, out_audit)
h    <- sha256_file(out_audit)
size <- file.info(out_audit)$size / 1024^2
registry_append("ALL", "Gene_Audit_PAM50", out_audit, h,
                "INTEGRO", SCRIPT_NAME, size)
message("\n[04_audit] Salvo: ", out_audit)

# Resumo de cobertura
out_cov <- file.path(PATHS$results$supp, "pam50_coverage_summary.csv")
write_csv(df_coverage, out_cov)
h    <- sha256_file(out_cov)
size <- file.info(out_cov)$size / 1024^2
registry_append("ALL", "PAM50_Coverage_Summary", out_cov, h,
                "INTEGRO", SCRIPT_NAME, size)
message("[04_audit] Salvo: ", out_cov)

# =============================================================================
# 7) GO / NO-GO FINAL
# =============================================================================
message("\n", strrep("=", 60))
message("[04_audit] VEREDICTO FINAL")
message(strrep("=", 60))

train_row <- df_coverage |> filter(role == "TRAIN")
if (nrow(train_row) > 0 && !is.na(train_row$n_pam50_missing)) {
  if (train_row$n_pam50_missing == 0) {
    message("GO — todos os 50 genes PAM50 presentes no SCANB (treino).")
    message("Proximo passo: scripts/05_reduce_pam50_to_corepam_FINAL.R")
  } else {
    message(sprintf("NO-GO — %d genes PAM50 ausentes no SCANB:",
                    train_row$n_pam50_missing))
    message("  ", train_row$genes_missing)
    message("Revisar 03_expression_preprocess_SCANB.R antes de prosseguir.")
  }
} else {
  message("INDETERMINADO — execute 03_expression_preprocess_SCANB.R primeiro.")
}

message("\n[04_audit] Concluido: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
