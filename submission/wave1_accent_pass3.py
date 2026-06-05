"""Pass 3: revert over-eager e→é conversions and add more missing accents."""
import re
QMD = 'TESE_CorePAM_v0.2.qmd'
t = open(QMD, encoding='utf-8').read()

# --- Revert wrong é→e in headings and specific contexts (conjunction usage) ---
reverts = [
    ('intrínsecos é o PAM50', 'intrínsecos e o PAM50'),
    ('expressão gênica é o desafio', 'expressão gênica e o desafio'),
    ('(24) é o ensaio', '(24) e o ensaio'),
    ('OncotypeDX é o MammaPrint', 'OncotypeDX e o MammaPrint'),
    ('proliferação) é o ROR', 'proliferação) e o ROR'),
    ('microarray é, posteriormente', 'microarray e, posteriormente'),
    ('expressão é, por meio', 'expressão e, por meio'),
    ('Análises posteriores do TAILORx (24) e o ensaio', 'Análises posteriores do TAILORx (24) e o ensaio'),  # noop confirmation
    ('desenvolvimento, validação interna, validação externa é modelo final', 'desenvolvimento, validação interna, validação externa e modelo final'),
    ('dependência é uso comercial', 'dependência e uso comercial'),
    ('dependência é necessidade', 'dependência e necessidade'),
    ('custo elevado (US$ 3.000 a 4.000 por amostra em laboratórios de referência nos EUA; indisponível no Sistema Único de Saúde brasileiro), dependência é', 'custo elevado (US$ 3.000 a 4.000 por amostra em laboratórios de referência nos EUA; indisponível no Sistema Único de Saúde brasileiro), dependência e'),
    ('desfechos a longo prazo discordantes', 'desfechos a longo prazo discordantes'),
    # Very common pattern: "X é Y" where Y starts with lowercase noun — LIKELY é
    # But "X e Y" where Y starts with uppercase OR is another noun — LIKELY conjunction
]
for a, b in reverts:
    t = t.replace(a, b)

# --- More accents ---
more = {
    'distincao': 'distinção', 'superexpressao': 'superexpressão',
    'contaminacao': 'contaminação', 'alteracoes': 'alterações',
    'saturacao': 'saturação', 'variacao': 'variação',
    'impoe': 'impõe', 'capitulo': 'capítulo',
    'consorcios': 'consórcios', 'academica': 'acadêmica',
    'asiatica': 'asiática', 'academico': 'acadêmico',
    'histologico': 'histológico', 'histologica': 'histológica',
    'bioligico': 'biológico',
    'dinamica': 'dinâmica', 'dinamico': 'dinâmico',
    'biobanked': 'biobanked',  # OK
    'analoga': 'análoga', 'analogo': 'análogo',
    'indisponivel': 'indisponível',
    'flexibilidade analitica': 'flexibilidade analítica',
    'perfil imuno-histoquimico': 'perfil imuno-histoquímico',
    'plataformas microarray e, posteriormente, RT-qPCR e': 'plataformas microarray e, posteriormente, RT-qPCR e',  # noop
    'hospitais universitarios': 'hospitais universitários',
    'laboratorio de implantacao': 'laboratório de implantação',
    'assinatura genetica': 'assinatura genética',
    'genetica': 'genética', 'genetico': 'genético',
    'geneticos': 'genéticos', 'geneticas': 'genéticas',
    'pre-processamento': 'pré-processamento',
    'prediçao': 'predição', 'predicao': 'predição',
    'estabilizacao': 'estabilização',
    'transformacao': 'transformação',
    'transformacoes': 'transformações',
    'Atribuicao': 'Atribuição',
    'atribuicao': 'atribuição',
    'reamostra': 'reamostra',  # OK
    'reamostras': 'reamostras',  # OK
    'arbitrariedade': 'arbitrariedade',  # OK
    'colineidaade': 'colinearidade',  # typo
    'Regularizacao': 'Regularização',
    'regularizacao': 'regularização',
    'Pre-processamento': 'Pré-processamento',
    'aplicacao': 'aplicação',
    'aplicacoes': 'aplicações',
    'implementacoes': 'implementações',
    'implementacao': 'implementação',
    'organizacao': 'organização',
    'conjuncao': 'conjunção',
    'tecnologica': 'tecnológica', 'tecnologico': 'tecnológico',
    'tres níveis': 'três níveis',
    'tres componentes': 'três componentes',
    'tres propriedades': 'três propriedades',
    'tres objetivos': 'três objetivos',
    'tres etapas': 'três etapas',
    'tres análises': 'três análises',
    'tres etapas explicitas': 'três etapas explícitas',
    'alta vazao': 'alta vazão',
    'vazao': 'vazão',
    'cobrertura': 'cobertura',
    'estabiliza-': 'estabiliza-',  # noop
    'estabilidade': 'estabilidade',  # OK
    'estatistica': 'estatística',
    'hierarquica': 'hierárquica',
    'principios': 'princípios',
    'principio': 'princípio',
    'relação é, por meio': 'relação e, por meio',
}
for k, v in more.items():
    if k != v:
        p = re.compile(r'(?<!\w)' + re.escape(k) + r'(?!\w)')
        t = p.sub(v, t)

# Additional strict safety: specific known-wrong é→e reverts
extra_reverts = [
    # "X e Y" patterns where Y starts with "o"/"a" but is conjunction context
    (' é o desafio da transportabilidade', ' e o desafio da transportabilidade'),
    # capacidade de detectar novos transcritos, maior sensibilidade — OK (no issue)
]
for a, b in extra_reverts:
    t = t.replace(a, b)

open(QMD, 'w', encoding='utf-8').write(t)
print('Pass 3 done. size:', len(t))
