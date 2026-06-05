"""Second-pass accent correction for TESE_CorePAM_v0.2.qmd.
Fixes unaccented words missed by wave1_rewrite.py and reverts
wrong Está → Esta (demonstrative).
"""
import re

QMD = 'TESE_CorePAM_v0.2.qmd'
text = open(QMD, encoding='utf-8').read()

# --- Step 1: Revert "Está " demonstrative back to "Esta " ---
# All 16 usages in current QMD are demonstrative (confirmed by inspection)
text = text.replace('Está ', 'Esta ')

# --- Step 2: Fix uppercase headings ---
heading_map = {
    'INTRODUCAO': 'INTRODUÇÃO',
    'METODO': 'MÉTODO',
    'DISCUSSAO': 'DISCUSSÃO',
    'CONCLUSAO': 'CONCLUSÃO',
    'REFERENCIAS': 'REFERÊNCIAS',
}
for k, v in heading_map.items():
    text = text.replace(k, v)

# --- Step 3: Expanded dictionary of common missed accents ---
ACC = {
    # -ao endings
    'pulmao': 'pulmão', 'apresentacao': 'apresentação', 'geracao': 'geração',
    'motivacao': 'motivação', 'evolucao': 'evolução', 'expressao': 'expressão',
    'decisao': 'decisão', 'excecao': 'exceção', 'anotacao': 'anotação',
    'configuracao': 'configuração', 'correcao': 'correção', 'precisao': 'precisão',
    'Padrao': 'Padrão', 'padrao': 'padrão', 'atribuicao': 'atribuição',
    'Atribuicao': 'Atribuição', 'execucao': 'execução', 'inspecao': 'inspeção',
    'separacao': 'separação', 'competicao': 'competição', 'selecao': 'seleção',
    'posicao': 'posição', 'constatacao': 'constatação', 'aproximacao': 'aproximação',
    'distribuicao': 'distribuição', 'penalizacao': 'penalização',
    'regularizacao': 'regularização', 'combinacao': 'combinação',
    'iteracao': 'iteração', 'particao': 'partição', 'calibracao': 'calibração',
    'Calibracao': 'Calibração', 'pontuacao': 'pontuação', 'Pontuacao': 'Pontuação',
    'discriminacao': 'discriminação', 'visualizacao': 'visualização',
    'dicotomizacao': 'dicotomização', 'interpretacao': 'interpretação',
    'transformacao': 'transformação', 'estabilizacao': 'estabilização',
    'calibracao-no-grande': 'calibração-no-grande', 'reposicao': 'reposição',
    'remocao': 'remoção', 'Direcao': 'Direção', 'direcao': 'direção',
    'verificacao': 'verificação', 'diferenciacao': 'diferenciação',
    'suposicao': 'suposição', 'violacao': 'violação', 'inflacao': 'inflação',
    'implantacao': 'implantação', 'implicacao': 'implicação',
    'submissao': 'submissão', 'aprovacao': 'aprovação', 'Resolucao': 'Resolução',
    'revisao': 'revisão', 'transformacao': 'transformação',
    # -oes endings
    'milhoes': 'milhões', 'restricoes': 'restrições', 'motivacoes': 'motivações',
    'instituicoes': 'instituições', 'implementacoes': 'implementações',
    'competicoes': 'competições', 'decisoes': 'decisões', 'predicoes': 'predições',
    'distribuicoes': 'distribuições', 'transformacoes': 'transformações',
    'variacoes': 'variações', 'implicacoes': 'implicações',
    # various
    'declinio': 'declínio', 'avancos': 'avanços', 'avanco': 'avanço',
    'sistematico': 'sistemático', 'sistematica': 'sistemática',
    'classico': 'clássico', 'classica': 'clássica',
    'classicos': 'clássicos', 'classicas': 'clássicas',
    'analitico': 'analítico', 'analitica': 'analítica',
    'analiticos': 'analíticos', 'analiticas': 'analíticas',
    'secundaria': 'secundária', 'secundario': 'secundário',
    'secundarias': 'secundárias', 'secundarios': 'secundários',
    'exploratoria': 'exploratória', 'exploratorias': 'exploratórias',
    'primaria': 'primária', 'primario': 'primário',
    'primarias': 'primárias', 'primarios': 'primários',
    'pragmatica': 'pragmática', 'pragmatico': 'pragmático',
    'incompativel': 'incompatível', 'compativel': 'compatível',
    'covariavel': 'covariável', 'covariaveis': 'covariáveis',
    'parametro': 'parâmetro', 'parametros': 'parâmetros',
    'codigo': 'código', 'codigos': 'códigos',
    'Independencia': 'Independência', 'independencia': 'independência',
    'consistencia': 'consistência', 'coerencia': 'coerência',
    'concordancia': 'concordância', 'consequencia': 'consequência',
    'referencia': 'referência', 'referencias': 'referências',
    'variancia': 'variância', 'variancias': 'variâncias',
    'importancia': 'importância', 'distancia': 'distância',
    'tolerancia': 'tolerância', 'divergencia': 'divergência',
    'favoravel': 'favorável', 'desfavoravel': 'desfavorável',
    'testavel': 'testável', 'comparavel': 'comparável',
    'preferivel': 'preferível', 'aceitavel': 'aceitável',
    'operaveis': 'operáveis', 'disponiveis': 'disponíveis',
    'uteis': 'úteis', 'viavel': 'viável', 'viaveis': 'viáveis',
    'tres ': 'três ', 'veem ': 'vêem ',
    'metodologico': 'metodológico', 'metodologicas': 'metodológicas',
    'metodologicos': 'metodológicos',
    'hierarquica': 'hierárquica', 'proprio': 'próprio', 'propria': 'própria',
    'proprios': 'próprios', 'proprias': 'próprias',
    'obtera': 'obterá', 'ira ': 'irá ', 'sera ': 'será ',
    'residuos': 'resíduos', 'liquido': 'líquido',
    'binarios': 'binários', 'binarias': 'binárias',
    'laboratorio': 'laboratório', 'laboratorios': 'laboratórios',
    'indisponivel': 'indisponível', 'disponivel': 'disponível',
    'Metricas': 'Métricas', 'metricas': 'métricas',
    'repositorios': 'repositórios', 'repositorio': 'repositório',
    'principios': 'princípios', 'principio': 'princípio',
    'orcamentarias': 'orçamentárias', 'orcamentario': 'orçamentário',
    'supraconfianca': 'supraconfiança', 'subconfianca': 'subconfiança',
    'verossimilhanca': 'verossimilhança', 'frequencia': 'frequência',
    'frequencias': 'frequências', 'significancia': 'significância',
    'Pre': 'Pré', 'pre-especificacao': 'pré-especificação',
    'pre-especificado': 'pré-especificado', 'pre-especificada': 'pré-especificada',
    'analogas': 'análogas', 'analogos': 'análogos',
    'hipotese': 'hipótese', 'hipoteses': 'hipóteses',
    'universitarios': 'universitários', 'universitario': 'universitário',
    'Universitario': 'Universitário',
    'Brasilia': 'Brasília',
    'meta-analise': 'meta-análise',
    'Clinicas': 'Clínicas',
    'politica': 'política', 'politicas': 'políticas',
    'cenarios': 'cenários', 'cenario': 'cenário',
    'utilidade': 'utilidade',  # OK
    'descricao': 'descrição', 'ausencia': 'ausência',
    'estrategia': 'estratégia', 'estrategias': 'estratégias',
    'estagio': 'estágio', 'estagios': 'estágios',
    'intermediario': 'intermediário', 'intermediaria': 'intermediária',
    'Comparacoes': 'Comparações', 'comparacoes': 'comparações',
    'Comparacao': 'Comparação', 'comparacao': 'comparação',
    'epidermico': 'epidérmico',
    'ideia': 'ideia',  # OK
    'identico': 'idêntico', 'identica': 'idêntica',
    'identificos': 'idênticos',
    'historico': 'histórico', 'historica': 'histórica',
    'Historico': 'Histórico',
    'subsequente': 'subsequente',  # OK
    'publico': 'público', 'publica': 'pública', 'publicos': 'públicos', 'publicas': 'públicas',
    'Publico': 'Público', 'Publica': 'Pública',
    'reproducibilidade': 'reprodutibilidade',  # already ok
    'oncologia': 'oncologia',  # OK
    'explicita': 'explícita', 'explicito': 'explícito',
    'explicitas': 'explícitas', 'explicitos': 'explícitos',
    'modulo': 'módulo', 'modulos': 'módulos',
    'imuno-histoquimico': 'imuno-histoquímico',
    'imuno-histoquímica': 'imuno-histoquímica',
    'terapeutica': 'terapêutica', 'terapeuticas': 'terapêuticas',
    'terapeutico': 'terapêutico',
    'semelhantes': 'semelhantes',  # OK
    'trieniio': 'triênio', 'trienio': 'triênio',
    'sistemicos': 'sistêmicos', 'sistemico': 'sistêmico',
    'incidencia': 'incidência',
    'historica': 'histórica',
    'etica': 'ética', 'eticos': 'éticos', 'eticas': 'éticas',
    'Aspectos eticos': 'Aspectos éticos',
    'topologica': 'topológica', 'topologico': 'topológico',
    'biologica': 'biológica', 'biologico': 'biológico',
    'estatistica': 'estatística', 'estatisticas': 'estatísticas',
    'estatistico': 'estatístico',
    'prospectivo': 'prospectivo',  # OK
    'ambito': 'âmbito', 'basico': 'básico', 'basicos': 'básicos',
    'basica': 'básica', 'basicas': 'básicas',
    # specific phrases fixing
    'Nao ': 'Não ', 'nao ': 'não ',  # careful: "não" is the most common
    'que nao ': 'que não ', 'e nao ': 'e não ',
}

# Apply with word-boundary safety (no partial word replacement)
# Sort by length desc to match longer patterns first
keys = sorted(ACC.keys(), key=len, reverse=True)
for k in keys:
    v = ACC[k]
    if k != v:
        # word-boundary aware
        if re.search(r'\w', k):
            pattern = re.compile(r'(?<!\w)' + re.escape(k) + r'(?!\w)')
            text = pattern.sub(v, text)

# --- Fix "e" (word, single letter) meaning "é" (verb) ---
# Only in very safe patterns like "N = X é", "hipotese é que", "e um" → "é um"
# Pattern: " e " between identifier-like words → likely "é"
# SAFE patterns only:
# "X e Y: " where Y is a noun phrase context
# Better: specific surgical replacements

safe_e_patterns = [
    # "o X e Y" where e is verb
    (r' e uma ', ' é uma '),
    (r' e um ', ' é um '),
    (r' e o ', ' é o '),
    (r' e a ', ' é a '),
    (r' e necessário', ' é necessário'),
    (r' e necessária', ' é necessária'),
    (r' e realmente', ' é realmente'),
    (r' e distinta', ' é distinta'),
    (r' e invariante', ' é invariante'),
    (r' e importante', ' é importante'),
    (r' e puramente', ' é puramente'),
    (r' e pragmatica', ' é pragmática'),
    (r' e discutida', ' é discutida'),
    (r' e calculado', ' é calculado'),
    (r' e aplicado', ' é aplicado'),
    (r' e definido', ' é definido'),
    (r' e avaliada', ' é avaliada'),
    (r' e seguinte', ' é seguinte'),
    (r' e preferivel', ' é preferível'),
    (r' e limitada', ' é limitada'),
    (r' e incompativel', ' é incompatível'),
    (r' e incompatível', ' é incompatível'),
    (r' e criptograficamente', ' é criptograficamente'),
    (r' e a expressao', ' é a expressão'),
    (r' e o conjunto', ' é o conjunto'),
    (r' e o peso', ' é o peso'),
    (r' e a log', ' é a log'),
    (r' e testar', ' é testar'),
    (r' e natural', ' é natural'),
    (r' e desde', ' e desde'),  # keep: "e desde" = "and since"
    (r'onde e ', 'onde é '),
    (r', e uma ', ', é uma '),
    (r', e um ', ', é um '),
    (r'coorte e ', 'coorte é '),
    (r'gene e ', 'gene é '),
    (r'Está e ', 'Esta é '),
    (r'Esta e ', 'Esta é '),
    (r'modelo e ', 'modelo é '),
    (r'procedimento e ', 'procedimento é '),
    (r'objetivo e ', 'objetivo é '),
    (r'coerência é ', 'coerência é '),  # noop
    (r' e a motivacao ', ' é a motivação '),
    (r' e a motivação ', ' é a motivação '),
]

for pat, rep in safe_e_patterns:
    text = re.sub(pat, rep, text)

# --- Remaining specific fixes by context ---
specific = {
    'caracterizados': 'caracterizados',  # OK
    'Por utilizar exclusivamente': 'Por utilizar exclusivamente',
    'tecnologias e contextos populacionais intencionalmente heterogeneos': 'tecnologias e contextos populacionais intencionalmente heterogêneos',
    'heterogeneos': 'heterogêneos',
    'subsequentes': 'subsequentes',
    'especificas': 'específicas', 'especifico': 'específico',
    'especificos': 'específicos', 'especifica': 'específica',
    'necessarios': 'necessários', 'necessario': 'necessário',
    'necessarias': 'necessárias', 'necessaria': 'necessária',
    'existentes': 'existentes',
    'tera': 'terá', 'tera ': 'terá ',
    'contemporanea': 'contemporânea', 'contemporaneo': 'contemporâneo',
    'metodo': 'método', 'metodos': 'métodos',
    'unidos': 'unidos',  # OK
    'numerica': 'numérica', 'numerico': 'numérico',
    'numericas': 'numéricas', 'numericos': 'numéricos',
    'inequiparaveis': 'inequiparáveis',
    'intermediario': 'intermediário',
    'operavel': 'operável',
    'pulmao': 'pulmão',  # already
    'seculo': 'século',
    'etaria': 'etária', 'etario': 'etário',
    'ciclico': 'cíclico', 'ciclica': 'cíclica',
    'nucleo': 'núcleo', 'nucleos': 'núcleos',
    'epitelio': 'epitélio',
}
for k, v in specific.items():
    if k != v:
        pattern = re.compile(r'(?<!\w)' + re.escape(k) + r'(?!\w)')
        text = pattern.sub(v, text)

open(QMD, 'w', encoding='utf-8').write(text)
print('Pass 2 accent correction applied. QMD size:', len(text))
