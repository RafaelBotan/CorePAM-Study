"""Fix remaining heading accent issues found in grep scan."""
QMD = 'TESE_CorePAM_v0.2.qmd'
t = open(QMD, encoding='utf-8').read()

fixes = [
    ('### Principios de validação externa', '### Princípios de validação externa'),
    ('### Atribuição deterministica de *folds* (SHA-256)', '### Atribuição determinística de *folds* (SHA-256)'),
    ('## Derivação reprodutível: inovacao metodologica', '## Derivação reprodutível: inovação metodológica'),
    ('## Coerencia biológica: por que esses 24 genes?', '## Coerência biológica: por que esses 24 genes?'),
    ('## Reflexoes metodológicas para estudos prognósticos futuros', '## Reflexões metodológicas para estudos prognósticos futuros'),
    ('## Relevancia para o contexto brasileiro é o SUS', '## Relevância para o contexto brasileiro e o SUS'),
    ('## Limitacoes', '## Limitações'),
    ('# REFERÊNCIAS BIBLIOGRAFICAS {.unnumbered}', '# REFERÊNCIAS BIBLIOGRÁFICAS {.unnumbered}'),
    ('# APENDICE A — Artigo submetido ao Breast Câncer Research {.unnumbered}',
     '# APÊNDICE A — Artigo submetido ao *Breast Cancer Research* {.unnumbered}'),
    ('# APENDICE B — Pipeline computacional {.unnumbered}', '# APÊNDICE B — Pipeline computacional {.unnumbered}'),
    ('## Pseudocodigo do procedimento de derivação', '## Pseudocódigo do procedimento de derivação'),
]
count = 0
for a, b in fixes:
    if a in t:
        t = t.replace(a, b)
        count += 1
        print('  FIXED:', a[:60])
    else:
        print('  SKIP (not found):', a[:60])

open(QMD, 'w', encoding='utf-8').write(t)
print(f'\n{count} heading fixes applied. size: {len(t)}')
