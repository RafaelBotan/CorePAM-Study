"""Pass 6: final accent fixes (Interacao, notavel, quimiossensiveis, etc.)."""
import re
QMD = 'TESE_CorePAM_v0.2.qmd'
t = open(QMD, encoding='utf-8').read()

more = {
    'Interacao': 'Interação',
    'interacao': 'interação',
    'Interacoes': 'Interações',
    'interacoes': 'interações',
    'notavel': 'notável',
    'Notavel': 'Notável',
    'notaveis': 'notáveis',
    'Notaveis': 'Notáveis',
    'quimiossensiveis': 'quimiossensíveis',
    'quimiossensivel': 'quimiossensível',
    'quimiorresistente': 'quimiorresistente',  # OK
    'agressivos': 'agressivos',  # OK
    'consistente': 'consistente',  # OK
    'sistemaicamente': 'sistematicamente',
    'responsiveness': 'responsiveness',  # EN
    # check common typos
    'Dicotomia': 'Dicotomia',  # OK
    'dicotomia': 'dicotomia',  # OK
    'efeitos': 'efeitos',  # OK
    'magnitudes': 'magnitudes',  # OK
    'concentrado': 'concentrado',  # OK
    'ausencia': 'ausência',
    'Ausencia': 'Ausência',
    'ausencias': 'ausências',
    'fortemente': 'fortemente',  # OK
    'disponibilidade': 'disponibilidade',  # OK
    'similaridade': 'similaridade',  # OK
    'necessariamente': 'necessariamente',  # OK
    # Common adverbs/adjectives with accent
    'possivelmente': 'possivelmente',  # OK
    'provavelmente': 'provavelmente',  # OK
    'obviamente': 'obviamente',  # OK
    'seguramente': 'seguramente',  # OK
    'normalmente': 'normalmente',  # OK
    'rapidamente': 'rapidamente',  # OK
    'relevante': 'relevante',  # OK
    'seguinte': 'seguinte',  # OK
    # important residuals
    'orgao': 'órgão',
    'orgaos': 'órgãos',
    'fisio': 'fisio',  # OK prefix
    'biologia': 'biologia',  # OK
    'biologico': 'biológico',
    'biologica': 'biológica',
    'biologicos': 'biológicos',
    'biologicas': 'biológicas',
    'tecnologia': 'tecnologia',  # OK
    'automatica': 'automática',
    'automatico': 'automático',
    'automaticas': 'automáticas',
    'automaticos': 'automáticos',
    'automaticamente': 'automaticamente',  # OK
    'publica': 'pública',
    'publico': 'público',
    'publicas': 'públicas',
    'publicos': 'públicos',
    # Pronouns / common verbs
    'voce': 'você',
    'voces': 'vocês',
    'esta': 'esta',  # OK (demonstrative; DO NOT touch — would break)
    # Future tenses
    'sera': 'sera',  # SKIP - ambiguous
    'explicitas': 'explícitas',
    'explicita': 'explícita',
    'explicitos': 'explícitos',
    'explicito': 'explícito',
    # Preservar
    'Preservar': 'Preservar',  # OK
    # Recem, ja, etc.
    'ja-': 'já-',  # careful - might be in compound words
    'nao-': 'não-',
    'Nao': 'Não',  # CAREFUL — only standalone, not "naos" etc
    # SKIP nao because could appear in "naoclinico" false positives
}

skip = {'esta', 'sera', 'ja-', 'nao-', 'Nao'}
for k in list(more.keys()):
    if k in skip:
        del more[k]

count_total = 0
for k, v in more.items():
    if k != v:
        p = re.compile(r'(?<!\w)' + re.escape(k) + r'(?!\w)')
        n = len(p.findall(t))
        if n > 0:
            t = p.sub(v, t)
            count_total += n

# Special: handle specific cases that are ambiguous
# "Nao" standalone at start of sentence: only replace if followed by space+lowercase word
# Too risky to do automatically. Skip.

open(QMD, 'w', encoding='utf-8').write(t)
print(f'Pass 6 done. {count_total} replacements.')
