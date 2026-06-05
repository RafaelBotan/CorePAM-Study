"""Pass 5: comprehensive scan + fix for common unaccented words."""
import re
QMD = 'TESE_CorePAM_v0.2.qmd'
t = open(QMD, encoding='utf-8').read()

# Common PT words without accents that MUST be accented
# Focused on those found in body text (not already caught)
more = {
    # Parametros / generico
    'Parametro': 'Parâmetro',
    'parametro': 'parâmetro',
    'Parametros': 'Parâmetros',
    'parametros': 'parâmetros',
    # Deterministica
    'Deterministica': 'Determinística',
    'deterministica': 'determinística',
    'Deterministico': 'Determinístico',
    'deterministico': 'determinístico',
    # Covariaveis / indisponiveis
    'Covariaveis': 'Covariáveis',
    'covariaveis': 'covariáveis',
    'Covariavel': 'Covariável',
    'covariavel': 'covariável',
    'indisponiveis': 'indisponíveis',
    # Other -ível/-íveis
    'sensiveis': 'sensíveis',
    'sensivel': 'sensível',
    'reprodutiveis': 'reprodutíveis',
    'reprodutivel': 'reprodutível',
    'estavel': 'estável',
    'estaveis': 'estáveis',
    'instavel': 'instável',
    'instaveis': 'instáveis',
    # Estatistica (already handled but variations)
    'Estatistica': 'Estatística',
    'estatistica': 'estatística',
    'estatisticas': 'estatísticas',
    'estatistico': 'estatístico',
    'estatisticos': 'estatísticos',
    'estatisticamente': 'estatisticamente',
    # Numero (already)
    'Numero': 'Número',
    # Mudanca / Balanco / etc (cê-cedilha)
    'mudanca': 'mudança',
    'mudancas': 'mudanças',
    'Mudanca': 'Mudança',
    'Mudancas': 'Mudanças',
    'Balanco': 'Balanço',
    'balanco': 'balanço',
    'comecar': 'começar',
    'comecou': 'começou',
    'comeca': 'começa',
    'comecam': 'começam',
    # Redundancia / dominio
    'redundancia': 'redundância',
    'Redundancia': 'Redundância',
    'dominio': 'domínio',
    'Dominio': 'Domínio',
    'dominios': 'domínios',
    # degradacao / composicao / preparacao / integracao
    'degradacao': 'degradação',
    'Degradacao': 'Degradação',
    'composicao': 'composição',
    'Composicao': 'Composição',
    'composicoes': 'composições',
    'preparacao': 'preparação',
    'Preparacao': 'Preparação',
    'integracao': 'integração',
    'Integracao': 'Integração',
    'integracoes': 'integrações',
    'correlacoes': 'correlações',
    'Correlacoes': 'Correlações',
    'correlacao': 'correlação',
    'Correlacao': 'Correlação',
    # outras palavras
    'intersecao': 'intersecção',
    'Intersecao': 'Intersecção',
    'plato': 'platô',  # if noun, but may be verb "plato" — careful
    'continuo': 'contínuo',
    'continua': 'contínua',  # careful, verb form too
    'continuos': 'contínuos',
    'continuas': 'contínuas',
    'performatico': 'performático',
    'performatica': 'performática',
    'obtem-se': 'obtém-se',
    'mantem-se': 'mantém-se',
    'tem-se': 'tem-se',  # OK
    'atraves': 'através',
    'excluido': 'excluído',
    'excluida': 'excluída',
    'excluidos': 'excluídos',
    'excluidas': 'excluídas',
    'incluida': 'incluída',
    'incluido': 'incluído',
    'incluidos': 'incluídos',
    'incluidas': 'incluídas',
    'identica': 'idêntica',
    'identico': 'idêntico',
    'identicas': 'idênticas',
    'identicos': 'idênticos',
    'ocorrera': 'ocorrerá',
    'ocorrerao': 'ocorrerão',
    'aplicara': 'aplicará',
    'aplicarao': 'aplicarão',
    'sera': 'será',
    # CAREFUL: "sera" may be proper name or conjugation — in scientific text mostly future verb
    'serao': 'serão',
    'estara': 'estará',
    'estarao': 'estarão',
    'podera': 'poderá',
    'poderao': 'poderão',
    'devera': 'deverá',
    'deverao': 'deverão',
    'fara': 'fará',
    'farao': 'farão',
    'tera': 'terá',
    'terao': 'terão',
    # carga, silabas com ^
    'enfase': 'ênfase',
    'fenomeno': 'fenômeno',
    'fenomenos': 'fenômenos',
    'economia': 'economia',  # OK
    'economico': 'econômico',
    'economica': 'econômica',
    'economicos': 'econômicos',
    'economicas': 'econômicas',
    'propoem': 'propõem',
    'dispoem': 'dispõem',
    'contem': 'contêm',  # plural; singular is "contém"
    # Careful with singular/plural. Let me not touch "contem" generically.
}

# Filter to remove problematic ambiguous ones
skip = {'contem', 'continua', 'plato', 'sera'}
for k in list(more.keys()):
    if k in skip:
        del more[k]

count_total = 0
for k, v in more.items():
    if k != v:
        p = re.compile(r'(?<!\w)' + re.escape(k) + r'(?!\w)')
        n_matches = len(p.findall(t))
        if n_matches > 0:
            t = p.sub(v, t)
            count_total += n_matches

open(QMD, 'w', encoding='utf-8').write(t)
print(f'Pass 5 done. {count_total} replacements. size: {len(t)}')
