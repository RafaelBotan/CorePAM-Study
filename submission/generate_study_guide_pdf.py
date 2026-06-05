"""
Generate PDF study guide for CorePAM with all figures and tables embedded.
Uses fpdf2 (pure Python, no system dependencies).
"""
import os
import re
from fpdf import FPDF

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(BASE)


# ── Figure registry (key → path, caption) ──────────────────────────────
def fig(relpath):
    """Resolve path relative to ROOT."""
    return os.path.join(ROOT, relpath.replace('/', os.sep))


def sfig(relpath):
    """Resolve path relative to BASE (submission/)."""
    return os.path.join(BASE, relpath.replace('/', os.sep))


FIGURE_BLOCKS = {
    # Main figures
    'fig1': {
        'path': sfig('figures/Fig1_StudyDesign_EN.png'),
        'caption': 'Figure 1. Study design overview. Derivacao em SCAN-B (N=3,069), validacao em 4 coortes externas, analise de pCR em 5 coortes.',
        'label': 'Figure 1 — Study Design',
    },
    'fig2': {
        'path': sfig('figures/Fig2_KM_MultiPanel_EN.png'),
        'caption': 'Figure 2. Kaplan-Meier dicotomizadas na mediana do CorePAM score (5 coortes: SCAN-B, TCGA, METABRIC, GSE20685, GSE1456).',
        'label': 'Figure 2 — KM Multi-Panel',
    },
    'fig3': {
        'path': sfig('figures/Fig4_Meta_Forest_HR_per1SD_CorePAM.png'),
        'caption': 'Figure 3. Forest plot meta-analise. HR por 1-SD CorePAM (K=4, REML). Pooled HR = 1.37 (1.24-1.52), I2 = 38.2%.',
        'label': 'Figure 3 — Forest Meta-Analysis',
    },
    'fig4a': {
        'path': sfig('figures/Fig5_DeltaCindex_COREA_vs_COREAplus_CorePAM.png'),
        'caption': 'Figure 4a. DeltaC-index incremental do CorePAM sobre CORE-A por coorte. Bootstrap CI (B=1000).',
        'label': 'Figure 4a — Delta C-index',
    },
    'fig4b': {
        'path': sfig('figures/Fig5_Calibration_60m_Panels_CorePAM.png'),
        'caption': 'Figure 4b. Calibracao: predicted vs. observed em 60m (SCAN-B, METABRIC, GSE20685) e 24m (TCGA).',
        'label': 'Figure 4b — Calibration',
    },
    'fig5': {
        'path': sfig('figures/Fig_pCR1_Forest_OR_EN.png'),
        'caption': 'Figure 5. Forest plot pCR: OR por 1-SD (4 coortes, N=697). Pooled OR = 1.69 (1.39-2.05), I2 = 0%.',
        'label': 'Figure 5 — pCR Forest',
    },
    # Supplementary figures (AF order from manuscript)
    'figs1': {
        'path': sfig('additional_files/Additional_file_09_FigS9_GeneDropout.png'),
        'caption': 'Figure S1 (AF1). Leave-one-gene-out sensitivity. Mudanca no C-index ao remover cada gene individualmente.',
        'label': 'Figure S1 — Gene Dropout',
    },
    'figs2': {
        'path': sfig('additional_files/Additional_file_05_FigS5_METABRIC_Sensitivity.png'),
        'caption': 'Figure S2 (AF2). METABRIC sensibilidade: A) OS, B) Fine-Gray competing risks, C) Quartile KM (DSS).',
        'label': 'Figure S2 — METABRIC Sensitivity',
    },
    'figs3': {
        'path': sfig('additional_files/Additional_file_03_FigS3_Correlation.png'),
        'caption': 'Figure S3 (AF4). Correlacao off-diagonal do CorePAM score entre coortes (Pearson r).',
        'label': 'Figure S3 — Cross-Cohort Correlation',
    },
    'figs4': {
        'path': sfig('additional_files/Additional_file_04_FigS4_PCA.png'),
        'caption': 'Figure S4 (AF5). PCA forense do METABRIC colorido por ER status.',
        'label': 'Figure S4 — METABRIC PCA',
    },
    'figs5': {
        'path': sfig('additional_files/Additional_file_02_FigS2_Weights.png'),
        'caption': 'Figure S5 (AF6). Pesos dos 24 genes CorePAM (lollipop). Vermelho: risco. Azul: protetor.',
        'label': 'Figure S5 — Gene Weights',
    },
    'figs6': {
        'path': sfig('additional_files/Additional_file_01_FigS1_Pareto.png'),
        'caption': 'Figure S6 (AF8). Fronteira de Pareto: gene count vs. OOF C-index. Linha vermelha: margem DeltaC = 0.010.',
        'label': 'Figure S6 — Pareto Frontier',
    },
    'figs7': {
        'path': sfig('additional_files/Additional_file_06_FigS6_TCGA_24m.png'),
        'caption': 'Figure S7 (AF9). TCGA 24-month sensitivity: KM e HR com horizonte restrito.',
        'label': 'Figure S7 — TCGA 24m Sensitivity',
    },
    'figs8': {
        'path': sfig('additional_files/Additional_file_07_FigS7_DCA_OS.png'),
        'caption': 'Figure S8 (AF10). Decision Curve Analysis (OS): net benefit vs. treat-all e treat-none.',
        'label': 'Figure S8 — DCA OS',
    },
    'figs9': {
        'path': sfig('additional_files/Additional_file_10_FigS10_COREA_Sensitivity.png'),
        'caption': 'Figure S9 (AF11). CORE-A sensitivity: HR CorePAM sob diferentes baselines clinicas.',
        'label': 'Figure S9 — COREA Sensitivity',
    },
    'figs10': {
        'path': sfig('additional_files/Additional_file_08_FigS8_Bootstrap.png'),
        'caption': 'Figure S10 (AF16). Bootstrap stability: frequencia de cada gene PAM50 em 200 refits.',
        'label': 'Figure S10 — Bootstrap Stability',
    },
    # Extra figures (from figures/ directory, for deeper study)
    'head2head': {
        'path': fig('figures/supp/en/png/FigS_HeadToHead_24vs50.png'),
        'caption': 'Figure S (AF18). Head-to-head: CorePAM 24 vs. 49 genes. DeltaC-index por coorte (intra e frozen).',
        'label': 'Head-to-Head 24 vs 49 Genes',
    },
    'table3': {
        'path': fig('figures/supp/en/png/FigT_Table3_Survival_Performance.png'),
        'caption': 'Table 1 (visual). Performance prognostica do CorePAM por coorte.',
        'label': 'Table 1 — Survival Performance',
    },
    'km_scanb': {
        'path': fig('figures/main/en/png/Fig3_KM_SCANB_OS_CorePAM.png'),
        'caption': 'KM SCAN-B (coorte de treino). Referencia visual — nao incluida no manuscrito como figura principal.',
        'label': 'KM SCAN-B (Training)',
    },
    'er_forest': {
        'path': fig('figures/supp/en/png/FigS_ER_Stratified_Forest.png'),
        'caption': 'Forest plot ER-estratificado: HR CorePAM em ER+ vs. ER- por coorte.',
        'label': 'ER-Stratified Forest',
    },
    'score_dist': {
        'path': fig('figures/supp/en/png/FigS_CorePAM_RawScore_Distribution_EN.png'),
        'caption': 'Distribuicao do score CorePAM bruto por coorte.',
        'label': 'Score Distribution by Cohort',
    },
    'score_by_er': {
        'path': fig('figures/supp/en/png/FigS_ScoreByER_ByCohort.png'),
        'caption': 'Distribuicao do CorePAM score por status ER e coorte.',
        'label': 'Score by ER Status',
    },
    'cindex_cohort': {
        'path': fig('figures/supp/en/png/FigS_Cindex_ByCohort.png'),
        'caption': 'C-index por coorte: CORE-A vs. CORE-A + CorePAM.',
        'label': 'C-index by Cohort',
    },
    'validation_forest': {
        'path': fig('figures/supp/en/png/FigS_Forest_HR_ValidationCohorts.png'),
        'caption': 'Forest plot: HR ajustado por coorte de validacao.',
        'label': 'Validation Forest HR',
    },
    'pam50_comparison': {
        'path': fig('figures/supp/en/png/FigS_CorePAM_vs_PAM50full_EN.png'),
        'caption': 'CorePAM (24 genes) vs. PAM50 full (50 genes): comparacao de scores.',
        'label': 'CorePAM vs PAM50 Full',
    },
    'pcr_roc': {
        'path': fig('figures/pcr/en/png/Fig_pCR2_ROC_EN.png'),
        'caption': 'Curvas ROC para predicao de pCR por coorte.',
        'label': 'pCR ROC Curves',
    },
    'pcr_quartile': {
        'path': fig('figures/pcr/en/png/Fig_pCR3_QuartileRate_EN.png'),
        'caption': 'Taxa de pCR por quartil do CorePAM score.',
        'label': 'pCR by Quartile',
    },
    'pcr_scoredist': {
        'path': fig('figures/pcr/en/png/Fig_pCR4_ScoreDist_EN.png'),
        'caption': 'Distribuicao do score em respondedores (pCR) vs. nao-respondedores.',
        'label': 'pCR Score Distribution',
    },
    'dca_pcr': {
        'path': fig('figures/supp/en/png/FigS_DCA_pCR_EN.png'),
        'caption': 'Decision Curve Analysis (pCR): net benefit para predicao de resposta.',
        'label': 'DCA pCR',
    },
}


class StudyGuidePDF(FPDF):
    """Custom PDF class with header/footer and utilities."""

    def __init__(self):
        super().__init__(orientation='P', unit='mm', format='A4')
        self.set_auto_page_break(auto=True, margin=20)
        # Register Unicode fonts (Windows system fonts)
        fonts_dir = 'C:/Windows/Fonts'
        self.add_font('Arial', '', os.path.join(fonts_dir, 'arial.ttf'))
        self.add_font('Arial', 'B', os.path.join(fonts_dir, 'arialbd.ttf'))
        self.add_font('Arial', 'I', os.path.join(fonts_dir, 'ariali.ttf'))
        self.add_font('Arial', 'BI', os.path.join(fonts_dir, 'arialbi.ttf'))

    def header(self):
        if self.page_no() > 1:
            self.set_font('Arial', 'I', 8)
            self.set_text_color(120, 120, 120)
            self.cell(0, 5, 'CorePAM Study Guide — Rafael Botan (UnB)', align='L')
            self.ln(8)

    def footer(self):
        self.set_y(-15)
        self.set_font('Arial', 'I', 8)
        self.set_text_color(120, 120, 120)
        self.cell(0, 10, f'p. {self.page_no()}', align='C')

    def add_cover(self):
        self.add_page()
        self.ln(60)
        self.set_font('Arial', 'B', 28)
        self.set_text_color(26, 58, 92)
        self.multi_cell(0, 12, 'Guia de Estudo Completo\nCorePAM', align='C')
        self.ln(10)
        self.set_font('Arial', '', 12)
        self.set_text_color(100, 100, 100)
        self.multi_cell(0, 7, 'A 24-gene PAM50-derived expression score\nwith cross-platform external validation\nfor breast cancer prognosis', align='C')
        self.ln(20)
        self.set_font('Arial', 'B', 12)
        self.set_text_color(50, 50, 50)
        self.cell(0, 8, 'Rafael de Negreiros Botan', align='C', new_x='LMARGIN', new_y='NEXT')
        self.set_font('Arial', '', 10)
        self.cell(0, 7, 'Programa de Ciencias Medicas — Universidade de Brasilia (UnB)', align='C', new_x='LMARGIN', new_y='NEXT')
        self.ln(15)
        self.set_text_color(140, 140, 140)
        self.set_font('Arial', 'I', 10)
        self.cell(0, 7, 'Preparacao para defesa de tese — Abril 2026', align='C', new_x='LMARGIN', new_y='NEXT')
        self.ln(30)
        self.set_font('Arial', '', 9)
        self.set_text_color(100, 100, 100)
        self.multi_cell(0, 5, 'Este guia contem todas as figuras principais e suplementares\ndo manuscrito na ordem em que aparecem no texto,\nacompanhadas de explicacoes detalhadas para estudo.', align='C')

    def section_title(self, title, level=1):
        """Add a section heading."""
        self.ln(4)
        if level == 1:
            self.set_font('Arial', 'B', 16)
            self.set_text_color(26, 58, 92)
            self.multi_cell(0, 8, title)
            self.set_draw_color(26, 58, 92)
            self.line(self.l_margin, self.get_y(), self.w - self.r_margin, self.get_y())
            self.ln(4)
        elif level == 2:
            self.set_font('Arial', 'B', 13)
            self.set_text_color(42, 90, 140)
            self.multi_cell(0, 7, title)
            self.set_draw_color(200, 200, 200)
            self.line(self.l_margin, self.get_y(), self.w - self.r_margin, self.get_y())
            self.ln(3)
        elif level == 3:
            self.set_font('Arial', 'B', 11)
            self.set_text_color(58, 106, 156)
            self.multi_cell(0, 6, title)
            self.ln(2)
        elif level == 4:
            self.set_font('Arial', 'B', 10)
            self.set_text_color(74, 122, 172)
            self.multi_cell(0, 6, title)
            self.ln(1)

    def body_text(self, text):
        """Add body paragraph."""
        self.set_font('Arial', '', 10)
        self.set_text_color(30, 30, 30)
        self.multi_cell(0, 5, text)
        self.ln(2)

    def bold_text(self, text):
        """Add bold paragraph."""
        self.set_font('Arial', 'B', 10)
        self.set_text_color(30, 30, 30)
        self.multi_cell(0, 5, text)
        self.ln(2)

    def bullet(self, text):
        """Add a bullet point."""
        self.set_font('Arial', '', 10)
        self.set_text_color(30, 30, 30)
        x = self.get_x()
        self.cell(6, 5, chr(8226))
        self.multi_cell(0, 5, text)
        self.ln(1)

    def add_figure(self, key, width_pct=0.90):
        """Add a figure from the registry."""
        info = FIGURE_BLOCKS.get(key)
        if not info:
            return
        path = info['path']
        if not os.path.exists(path):
            self.set_font('Arial', 'I', 9)
            self.set_text_color(180, 0, 0)
            self.cell(0, 6, f'[Figura nao encontrada: {key} — {path}]', new_x='LMARGIN', new_y='NEXT')
            return

        avail_w = self.w - self.l_margin - self.r_margin
        img_w = avail_w * width_pct

        # Check if we need a new page (estimate ~100mm for image + caption)
        if self.get_y() > 200:
            self.add_page()

        # Gray background box
        self.set_fill_color(248, 249, 252)
        self.set_draw_color(210, 215, 225)

        # Center the image
        x_offset = self.l_margin + (avail_w - img_w) / 2

        self.ln(3)
        # Label
        self.set_font('Arial', 'B', 9)
        self.set_text_color(26, 58, 92)
        self.cell(0, 5, info['label'], new_x='LMARGIN', new_y='NEXT')
        self.ln(1)

        try:
            self.image(path, x=x_offset, w=img_w)
        except Exception as e:
            self.set_font('Arial', 'I', 9)
            self.set_text_color(180, 0, 0)
            self.cell(0, 6, f'[Erro ao carregar imagem: {e}]', new_x='LMARGIN', new_y='NEXT')
            return

        self.ln(2)
        # Caption
        self.set_font('Arial', 'I', 8)
        self.set_text_color(90, 90, 90)
        self.multi_cell(0, 4, info['caption'])
        self.ln(4)
        self.set_text_color(30, 30, 30)

    def add_simple_table(self, headers, rows):
        """Add a simple table."""
        self.set_font('Arial', '', 8)
        n_cols = len(headers)
        avail_w = self.w - self.l_margin - self.r_margin
        col_w = avail_w / n_cols

        # Header
        self.set_fill_color(26, 58, 92)
        self.set_text_color(255, 255, 255)
        self.set_font('Arial', 'B', 8)
        for h in headers:
            self.cell(col_w, 6, h, border=1, fill=True, align='C')
        self.ln()

        # Rows
        self.set_text_color(30, 30, 30)
        self.set_font('Arial', '', 8)
        for i, row in enumerate(rows):
            if i % 2 == 0:
                self.set_fill_color(245, 247, 250)
            else:
                self.set_fill_color(255, 255, 255)
            for cell in row:
                self.cell(col_w, 5, str(cell), border=1, fill=True, align='C')
            self.ln()
        self.ln(3)


def build_pdf():
    """Build the complete study guide PDF."""
    pdf = StudyGuidePDF()
    pdf.add_cover()

    # ══════════════════════════════════════════════════════════════
    # PARTE I — O QUE O ARTIGO DIZ
    # ══════════════════════════════════════════════════════════════
    pdf.add_page()
    pdf.section_title('PARTE I — O QUE O ARTIGO DIZ')

    # ── 1. TITULO E CLAIM ──
    pdf.section_title('1. Titulo e Claim Central', 2)
    pdf.bold_text('Titulo: "CorePAM: a 24-gene PAM50-derived expression score with cross-platform external validation for breast cancer prognosis"')
    pdf.body_text('Claim central: CorePAM e a menor reducao data-driven do PAM50 que satisfaz um criterio pre-especificado de nao-inferioridade (DeltaC <= 0.010) relativo a um modelo elastic-net Cox com todos os 50 genes, validada em 4 coortes externas (RNA-seq e microarray).')
    pdf.bold_text('O que o artigo NAO diz:')
    pdf.bullet('NAO diz que CorePAM e equivalente ao Prosigna comercial')
    pdf.bullet('NAO diz que CorePAM pode substituir testes genomicos existentes')
    pdf.bullet('NAO diz que CorePAM tem utilidade clinica demonstrada')
    pdf.bullet('NAO diz que o frozen z-score funciona igualmente em todas as plataformas')

    # ── 2. ABSTRACT ──
    pdf.section_title('2. Abstract — Linha por Linha', 2)

    pdf.section_title('Background', 3)
    pdf.body_text('> "The PAM50 classifier predicts breast cancer prognosis but requires 50 genes and specialised platforms."')
    pdf.bullet('PAM50: 50 genes (Perou 2000, Parker 2009) para subtipos intrinsecos')
    pdf.bullet('Specialised platforms: NanoString nCounter (unica plataforma do Prosigna)')

    pdf.body_text('> "We derived CorePAM, the smallest data-driven PAM50 subset maintaining non-inferior prognostic performance relative to a full 50-gene Cox elastic-net model, without pre-specifying gene count."')
    pdf.bullet('"Data-driven": o numero 24 veio dos dados, nao de decisao humana')
    pdf.bullet('"Non-inferior": dentro de margem pre-especificada (0.010 em C-index)')
    pdf.bullet('"Relative to 50-gene Cox elastic-net": comparador e mesmo framework, NAO o Prosigna')

    pdf.add_figure('fig1')

    pdf.section_title('Methods', 3)
    pdf.bullet('alpha = 0.5: Mixing parameter elastic-net (0=Ridge, 1=LASSO)')
    pdf.bullet('10-fold CV deterministico: SHA-256 do patient ID -> fold assignment')
    pdf.bullet('DeltaC = 0.010: Margem de nao-inferioridade (< 0.5 SE do C-index OOF)')
    pdf.bullet('N = 3,069 SCAN-B: Maior coorte publica de RNA-seq de mama')
    pdf.bullet('4 coortes externas: TCGA (RNA-seq), METABRIC (Illumina), GSE20685 (Affy), GSE1456 (Affy)')

    pdf.section_title('Results', 3)
    pdf.bullet('OOF C-index = 0.670 (24 genes) vs. C-index max = 0.679 (50 genes). Gap = 0.009')
    pdf.bullet('HRs: 1.20 (TCGA) a 1.71 (GSE1456) — todas coortes HR > 1')
    pdf.bullet('Meta-analise: HR = 1.37 (1.24-1.52), p = 1.6e-9, I2 = 38.2%')
    pdf.bullet('pCR: OR = 1.69 (1.39-2.05), I2 = 0%')

    pdf.section_title('Conclusions (v3 final)', 3)
    pdf.bullet('"Met the pre-specified DeltaC non-inferiority criterion"')
    pdf.bullet('"Retained prognostic discrimination" (linguagem conservadora)')
    pdf.bullet('"Was associated with pCR" (nao "predicts")')
    pdf.bullet('"May simplify future assay development" (condicional)')
    pdf.bullet('"Platform-specific analytical validation is required" (caveat explicito)')

    # ── 3. BACKGROUND ──
    pdf.add_page()
    pdf.section_title('3. Background — Contexto e Justificativa', 2)

    pdf.section_title('3.1 O que e o PAM50', 3)
    pdf.bullet('Perou 2000: subtipos moleculares por expressao genica')
    pdf.bullet('Parker 2009: refinamento para 50 genes (LumA, LumB, HER2-e, Basal, Normal)')
    pdf.bullet('Prosigna: implementacao comercial (NanoString), aprovado FDA 2013, usa ROR score')
    pdf.bullet('Utilidade: decisao sobre quimioterapia em HR+/HER2- (similar ao OncotypeDX)')

    pdf.section_title('3.2 O problema', 3)
    pdf.bullet('50 genes e caro: NanoString nCounter requer equipamento dedicado')
    pdf.bullet('Barreira de acesso: custo e logistica limitam uso em paises como o Brasil')
    pdf.bullet('Redundancia biologica: genes co-expressos em clusters de proliferacao/diferenciacao')

    pdf.section_title('3.3 O que CorePAM faz de diferente', 3)
    pdf.bullet('Data-driven: numero de genes emerge dos dados (nao hipotese)')
    pdf.bullet('Elastic-net: lida com correlacao entre genes (grouping effect)')
    pdf.bullet('Criterio pre-especificado: DeltaC = 0.010 definido ANTES da analise')
    pdf.bullet('Validacao cross-platform: RNA-seq + 3 plataformas de microarray')
    pdf.bullet('Preprocessamento independente: cada coorte processada separadamente')

    # ── 4. METHODS ──
    pdf.add_page()
    pdf.section_title('4. Methods — Cada Detalhe', 2)

    pdf.section_title('4.1 Coortes de Sobrevivencia', 3)
    pdf.add_simple_table(
        ['Coorte', 'Papel', 'Plataforma', 'Endpoint', 'N', 'Eventos', 'FU (mo)'],
        [
            ['SCAN-B', 'Treino', 'RNA-seq', 'OS', '3,069', '322', '54.9'],
            ['TCGA', 'Validacao', 'RNA-seq', 'OS', '1,072', '150', '32.0'],
            ['METABRIC', 'Validacao', 'Illumina', 'DSS', '1,978', '646', '159.0'],
            ['GSE20685', 'Validacao', 'Affymetrix', 'OS', '327', '83', '112.8'],
            ['GSE1456', 'Validacao', 'Affymetrix', 'OS', '159', '40', '91.0'],
        ]
    )
    pdf.body_text('SCAN-B: Suecia, Brueffer 2018. TCGA: EUA multicentrico. METABRIC: UK/Canada, Curtis 2012. GSE20685: Taiwan, Kao 2011 (Affy HG-U133 Plus 2.0 / GPL570 — CORRIGIDO). GSE1456: Stockholm, Pawitan 2005 (HGU133A, 21/24 genes, sem covariaveis).')

    pdf.section_title('4.2 Coortes pCR', 3)
    pdf.add_simple_table(
        ['Coorte', 'Plataforma', 'N', 'pCR rate'],
        [
            ['GSE25066', 'HGU133A', '182', '23.1%'],
            ['GSE20194', 'HGU133Plus2', '278', '20.1%'],
            ['GSE32646', 'HGU133Plus2', '115', '23.5%'],
            ['I-SPY1', 'Agilent 44K', '122', '26.2%'],
            ['I-SPY2', 'Agilent 44K', '986', '32.4%'],
        ]
    )

    pdf.section_title('4.3 Analytical Freeze', 3)
    pdf.add_simple_table(
        ['Parametro', 'Valor', 'Justificativa'],
        [
            ['DeltaC', '0.010', '< 0.5 SE C-index OOF'],
            ['alpha', '0.5', 'Elastic-net 50/50'],
            ['K folds', '10', 'Standard'],
            ['Fold assign', 'SHA-256', 'Determinismo total'],
            ['Min genes', '80%', 'Permite genes faltantes'],
            ['Gene count', 'Nao fixado', 'Emerge da Pareto'],
        ]
    )
    pdf.body_text('Todos os parametros congelados ANTES da analise. Equivalente a protocolo de ensaio clinico — elimina data dredging.')

    pdf.section_title('4.4 Derivacao do CorePAM', 3)
    pdf.body_text('Input: SCAN-B 3,069 x 50 PAM50 genes (z-scored). Modelo: Cox elastic-net (alpha=0.5, glmnet). CV: 10-fold estratificado. Regularization path: para cada lambda, computa OOF C-index e df (genes != 0). Fronteira de Pareto: menor df onde C(df) >= C(max) - 0.010.')
    pdf.bold_text('Resultado: df = 24, C-index = 0.670, C-max = 0.679, gap = 0.009')
    pdf.add_figure('figs6')

    pdf.section_title('4.5 Calculo do Score', 3)
    pdf.body_text('score = sum(w_i * z_i) / sum(|w_i|)')
    pdf.bullet('w_i: peso elastic-net (frozen do SCAN-B)')
    pdf.bullet('z_i: z-score intra-coorte (media=0, SD=1 dentro da coorte)')
    pdf.bullet('Denominador: soma |w_i| sobre genes presentes (escala invariante a genes faltantes)')
    pdf.bullet('Score positivo = maior risco (confirmado em todas as 5 coortes)')
    pdf.add_figure('figs5')

    pdf.section_title('4.6 Metodos Estatisticos', 3)
    pdf.add_simple_table(
        ['Metodo', 'Uso', 'O que mede'],
        [
            ['Cox univariado', 'HR por 1-SD', 'Associacao score-sobrevida'],
            ['Cox multivariado', 'HR ajustado CORE-A', 'Independencia clinica'],
            ['Kaplan-Meier', 'Curvas mediana', 'Visualizacao do efeito'],
            ['Harrell C-index', 'Discriminacao', 'Ordenacao de pacientes'],
            ["Uno's C (IPCW)", 'Robusto a censoring', 'Sensibilidade'],
            ['Meta-analise RE', 'HR poolado, I2', 'Sintese de evidencia'],
            ['Hartung-Knapp', 'Meta conservadora', 'Sensibilidade K pequeno'],
            ['Fine-Gray', 'Competing risks', 'METABRIC DSS'],
            ['Calibracao logist', 'Slope/intercept', 'Pred vs observado'],
            ['DCA', 'Net benefit', 'Utilidade exploratoria'],
            ['Bootstrap DeltaC', 'Valor incremental', 'Score alem de clinica'],
        ]
    )

    pdf.section_title('4.7 Frozen Z-Score — Conceito Critico', 3)
    pdf.body_text('Problema: z-score intra-coorte requer a coorte inteira. Na clinica, voce tem UM paciente.')
    pdf.body_text('Solucao: congelar media e SD do SCAN-B e aplicar a cada paciente externo: z_frozen = (expr - media_SCANB) / SD_SCANB')
    pdf.add_simple_table(
        ['Coorte', 'DeltaC frozen', 'Pearson r'],
        [
            ['TCGA (RNA-seq)', '+0.0002', '0.999'],
            ['METABRIC (Illumina)', '+0.0007', '0.958'],
            ['GSE20685 (Affy)', '-0.036', '0.923'],
            ['GSE1456 (Affy)', '-0.040', '0.916'],
        ]
    )
    pdf.bold_text('Frozen scoring e platform-conditional, NAO platform-agnostic.')

    # ── 5. RESULTS ──
    pdf.add_page()
    pdf.section_title('5. Results — Cada Numero Explicado', 2)

    pdf.section_title('5.1 Derivacao', 3)
    pdf.body_text('24 genes selecionados. Gap = 0.009 < 0.010 (margem). Sensibilidade: 0.005 -> 31 genes; 0.015 -> 18; 0.020 -> 16.')
    pdf.add_figure('score_dist')

    pdf.section_title('5.2 Head-to-head 24 vs 49 genes (NOVO em R2)', 3)
    pdf.add_simple_table(
        ['Coorte', 'Modo', 'DeltaC', '95% CI', 'Veredicto'],
        [
            ['TCGA', 'intra', '-0.008', '(-0.038, +0.022)', 'Nao-inferior'],
            ['TCGA', 'frozen', '-0.010', '(-0.039, +0.020)', 'Nao-inferior'],
            ['METABRIC', 'intra', '+0.022', '(+0.011, +0.032)', '24-gene SUPERIOR'],
            ['METABRIC', 'frozen', '+0.017', '(+0.007, +0.027)', '24-gene SUPERIOR'],
            ['GSE20685', 'intra', '+0.055', '(+0.010, +0.101)', '24-gene SUPERIOR'],
            ['GSE20685', 'frozen', '+0.047', '(+0.007, +0.085)', '24-gene SUPERIOR'],
            ['GSE1456', 'intra', '+0.097', '(+0.032, +0.166)', '24-gene SUPERIOR'],
            ['GSE1456', 'frozen', '+0.111', '(+0.045, +0.181)', '24-gene SUPERIOR'],
        ]
    )
    pdf.bold_text('Por que o modelo menor e MELHOR? Principio de regularizacao: modelos parcimoniosos generalizam melhor cross-platform (capturam menos ruido plataforma-especifico).')
    pdf.add_figure('head2head')

    pdf.section_title('5.3 Validacao Prognostica (Table 1)', 3)
    pdf.add_figure('table3')
    pdf.body_text('TCGA HR = 1.20 (mais baixo): follow-up mediano de 32 meses e curto para OS em mama. Sensibilidade 24 meses: HR = 1.64 (1.24-2.16, p=4.7e-4).')
    pdf.body_text('METABRIC usa DSS (nao OS): follow-up de 159 meses, muitas mortes por outras causas. DSS isola efeito do cancer. Fine-Gray como sensibilidade.')
    pdf.add_figure('fig2')
    pdf.add_figure('km_scanb')
    pdf.add_figure('figs2')

    pdf.section_title('5.4 Meta-Analise', 3)
    pdf.add_simple_table(
        ['Analise', 'HR (95% CI)', 'p', 'I2'],
        [
            ['Primary RE (REML)', '1.37 (1.24-1.52)', '1.6e-9', '38.2%'],
            ['Hartung-Knapp', '1.37 (1.15-1.64)', '0.011', '38.2%'],
            ['OS-harmonised', '1.36 (1.13-1.64)', '<0.001', '50.5%'],
        ]
    )
    pdf.body_text('I2 = 38.2% (moderado): fontes identificadas — endpoint DSS vs OS e follow-up 32 vs 159 meses. Achado robusto: direcao positiva em 4/4 coortes.')
    pdf.add_figure('fig3')

    pdf.section_title('5.5 Valor Incremental (CORE-A)', 3)
    pdf.add_simple_table(
        ['Coorte', 'CORE-A', 'DeltaC'],
        [
            ['SCAN-B', 'age + ER', '+0.030'],
            ['TCGA', 'age', '+0.038'],
            ['METABRIC', 'age + ER', '+0.043'],
            ['GSE20685', 'age', '+0.094'],
        ]
    )
    pdf.body_text('CORE-A: modelo clinico minimo (age +/- ER). GSE20685 tem maior DeltaC porque o modelo clinico (so idade) e fraco — mais espaco para score molecular.')
    pdf.add_figure('fig4a')
    pdf.add_figure('cindex_cohort')

    pdf.section_title('5.6 Independencia do Estadiamento', 3)
    pdf.body_text('SCAN-B (6 vars): HR adj = 1.75, LRT p = 4.9e-13, DeltaC = +0.017')
    pdf.body_text('METABRIC (6 vars): HR adj = 1.25, LRT p = 6.1e-6, DeltaC = +0.015')
    pdf.body_text('METABRIC ER+: DeltaC = +0.032 (mais pronunciado)')
    pdf.bold_text('CorePAM captura heterogeneidade molecular NAO acessivel pelo estadiamento clinico.')
    pdf.add_figure('validation_forest')
    pdf.add_figure('figs9')

    pdf.section_title('5.7 pCR (Analise Secundaria)', 3)
    pdf.add_simple_table(
        ['Coorte', 'N', 'pCR%', 'OR (95% CI)', 'AUC'],
        [
            ['GSE25066', '182', '23.1%', '1.44 (0.98-2.09)', '0.576'],
            ['GSE20194', '278', '20.1%', '1.73 (1.25-2.38)', '0.660'],
            ['GSE32646', '115', '23.5%', '2.10 (1.31-3.39)', '0.716'],
            ['I-SPY1', '122', '26.2%', '1.66 (1.05-2.63)', '0.622'],
            ['RE pooled', '697', '', '1.69 (1.39-2.05)', ''],
        ]
    )
    pdf.body_text('I2 = 0%: ausencia total de heterogeneidade (extraordinario). I-SPY2 exploratoria (N=986, OR=1.68) confirma.')
    pdf.add_figure('fig5')
    pdf.add_figure('pcr_roc')
    pdf.add_figure('pcr_quartile')

    pdf.section_title('5.8 Sensibilidades', 3)
    pdf.bold_text('Calibracao formal:')
    pdf.add_simple_table(
        ['Coorte', 'Horizonte', 'alpha', 'beta', 'E/O', 'Brier'],
        [
            ['SCAN-B', '60m', '0.56', '0.96', '0.61', '0.158'],
            ['TCGA', '24m', '5.40', '2.74', '0.74', '0.061'],
            ['METABRIC', '60m', '1.15', '1.77', '0.97', '0.143'],
            ['GSE20685', '60m', '0.93', '1.57', '0.98', '0.130'],
        ]
    )
    pdf.body_text('beta ideal = 1. SCAN-B quase ideal (0.96). Externos: slopes > 1 (sub-confianca). Discriminacao preservada; risco absoluto precisa recalibracao por plataforma.')
    pdf.add_figure('fig4b')

    pdf.bold_text('Subgrupo ER:')
    pdf.bullet('ER+: sinal robusto (HR 1.59-1.92)')
    pdf.bullet('ER-: sinal concentrado no SCAN-B (HR=2.53, MAS instavel: 35 eventos), ausente no METABRIC (HR=0.97, 186 eventos)')
    pdf.bullet('Conclusao: sinal clinico concentrado em ER+, consistente com biologia luminal')
    pdf.add_figure('er_forest')
    pdf.add_figure('score_by_er')

    # ── 6. DISCUSSION ──
    pdf.add_page()
    pdf.section_title('6. Discussion — Estrutura e Logica', 2)
    pdf.section_title('Paragrafo 1: Framework e scope', 3)
    pdf.bullet('Ancora nao-inferioridade ao elastic-net framework')
    pdf.bullet('Head-to-head como suporte (nao manchete)')
    pdf.bullet('Caveat: calibracao e frozen sao platform-dependent')
    pdf.bullet('"No claims of therapeutic equivalence" — frase protetora')

    pdf.section_title('Reproducible derivation', 3)
    pdf.bullet('SHA-256 = fold deterministico')
    pdf.bullet('Gene count (24) emerge da Pareto frontier')

    pdf.section_title('Gene selection stability', 3)
    pdf.bullet('Bootstrap B=200: 15/24 genes em >70% dos refits')
    pdf.bullet('2 genes (ACTR3B, NAT1) em 100%')
    pdf.bullet('Leave-one-gene-out: max |DeltaC| = 0.044')
    pdf.add_figure('figs10')
    pdf.add_figure('figs1')
    pdf.add_figure('pam50_comparison')

    pdf.section_title('Heterogeneidade', 3)
    pdf.bullet('I2 = 38.2% explicado por DSS vs OS + follow-up')
    pdf.bullet('OS-harmonised confirma')
    pdf.bullet('Direcao positiva em 4/4: achado robusto')

    pdf.section_title('Coerencia biologica', 3)
    pdf.bullet('Luminal: ESR1, PGR, BCL2, NAT1')
    pdf.bullet('Proliferacao: MYBL2, PTTG1, EXO1, MYC')
    pdf.bullet('Basal: KRT5, KRT17, FOXC1')
    pdf.bullet('HER2: ERBB2, GRB7')
    pdf.bullet('26 genes excluidos: peso zero no lambda selecionado')

    pdf.section_title('Calibracao e transportabilidade', 3)
    pdf.bullet('DCA exploratoria')
    pdf.bullet('Slopes 1.57-2.74 em externas')
    pdf.add_figure('figs8')
    pdf.add_figure('dca_pcr')

    pdf.section_title('Frozen z-score sensitivity', 3)
    pdf.bullet('Platform-conditional, nao platform-agnostic')
    pdf.bullet('Recalibracao por plataforma necessaria')
    pdf.add_figure('figs7')

    pdf.section_title('Limitacoes', 3)
    pdf.bullet('Study design: retrospectivo, dados publicos, TCGA curto, METABRIC DSS, I2 moderado')
    pdf.bullet('Analytical: CORE-A minimo, z-score intra-coorte, frozen platform-dependent')
    pdf.bullet('Clinical translation: validacao analitica e calibracao fora do escopo')

    # ── 7. 24 GENES ──
    pdf.add_page()
    pdf.section_title('7. Os 24 Genes — Biologia Detalhada', 2)

    genes = [
        ['EXO1', '+0.217', 'Risco', 'Prolif', 'Exonuclease reparo DNA'],
        ['NAT1', '-0.205', 'Protetor', 'Luminal', 'N-acetiltransferase'],
        ['BLVRA', '-0.185', 'Protetor', 'Outro', 'Biliverdina redutase'],
        ['ACTR3B', '-0.141', 'Protetor', 'Outro', 'Dinamica citoesqueleto'],
        ['MIA', '-0.119', 'Protetor', 'HER2', 'Melanoma inhib. activity'],
        ['MYBL2', '-0.118', 'Protetor*', 'Prolif', 'Fator transcricao ciclo cel.'],
        ['PTTG1', '+0.114', 'Risco', 'Prolif', 'Securina, separacao cromos.'],
        ['MDM2', '-0.107', 'Protetor', 'Outro', 'Ubiquitina ligase p53'],
        ['SFRP1', '-0.091', 'Protetor', 'Outro', 'Wnt antagonist'],
        ['GPR160', '-0.084', 'Protetor', 'Outro', 'Receptor acoplado Prot G'],
        ['FOXC1', '+0.072', 'Risco', 'Basal', 'Fator forkhead, agressivo'],
        ['PHGDH', '+0.062', 'Risco', 'Outro', 'Via sintese serina'],
        ['MYC', '-0.042', 'Protetor*', 'Prolif', 'Proto-oncogene'],
        ['ESR1', '+0.040', 'Risco*', 'Luminal', 'Receptor estrogeno alpha'],
        ['KRT5', '-0.028', 'Protetor', 'Basal', 'Citoqueratina 5'],
        ['CXXC5', '+0.027', 'Risco', 'Outro', 'Regulador epigenetico'],
        ['PGR', '-0.026', 'Protetor', 'Luminal', 'Receptor progesterona'],
        ['KRT17', '-0.016', 'Protetor', 'Basal', 'Citoqueratina 17'],
        ['CENPF', '+0.015', 'Risco', 'Prolif', 'Proteina centromerica'],
        ['FGFR4', '+0.011', 'Risco', 'Outro', 'Receptor FGF'],
        ['BCL2', '-0.010', 'Protetor', 'Luminal', 'Anti-apoptose (paradoxo)'],
        ['MLPH', '-0.010', 'Protetor', 'Luminal', 'Melanofilina'],
        ['ERBB2', '-0.008', 'Protetor*', 'HER2', 'HER2/neu'],
        ['GRB7', '-0.005', 'Protetor', 'HER2', 'Co-amplificado ERBB2'],
    ]
    pdf.add_simple_table(['Gene', 'Peso', 'Direcao', 'Cluster', 'Funcao'], genes)
    pdf.body_text('*Pesos contra-intuitivos: em contexto multivariado, o peso pode inverter sinal vs univariado. Ex: MYBL2 protetor porque EXO1 ja captura proliferacao-risco; MYBL2 captura componente residual (diferenciacao). ESR1 positivo captura expressao alta de ER em tumores com falha de resposta hormonal.')

    pdf.add_figure('figs3')
    pdf.add_figure('figs4')

    # ══════════════════════════════════════════════════════════════
    # PARTE II — COMO DEFENDER O ARTIGO
    # ══════════════════════════════════════════════════════════════
    pdf.add_page()
    pdf.section_title('PARTE II — COMO DEFENDER O ARTIGO')

    # ── 9. PERGUNTAS DE BANCA ──
    pdf.section_title('9. Perguntas de Banca', 2)

    pdf.section_title('Metodologia', 3)

    all_qa = [
        ('Metodologia', [
            ('Por que elastic-net e nao LASSO puro?',
             'LASSO (alpha=1) tem selecao instavel com genes correlacionados. Elastic-net (alpha=0.5) distribui pesos (grouping effect, Zou & Hastie 2005). Bootstrap confirma: 15/24 em >70% dos refits.'),
            ('Por que nao Random Forest ou Deep Learning?',
             '(1) Interpretabilidade: coeficientes lineares. (2) Selecao embutida: modelos esparsos. (3) Reproducibilidade: deterministico. (4) N=3,069 com 50 features: regiao onde regularizacao e robusta; DL requer ordens de grandeza mais dados.'),
            ('A margem de 0.010 nao e arbitraria?',
             'Toda margem de NI e escolha a priori (como em ensaios clinicos). 0.010 < 0.5 SE do C-index OOF (SE = 0.026). Sensibilidade: 0.005/0.015/0.020 -> 31/18/16 genes. Selecao nao e brittle.'),
            ('Z-score intra-coorte e um problema?',
             'E o framework padrao para validacao retrospectiva. O frozen z-score (4.7) simula deployment: funciona em RNA-seq/Illumina (DeltaC < 0.001), atenua em Affy (DeltaC ~ -0.04). Recalibracao necessaria.'),
            ('SHA-256 com cross-validation?',
             'Fold assignment deterministico e reprodutivel. Qualquer pessoa com mesmos IDs reproduz EXATAMENTE os mesmos folds. Elimina random seed como driver oculto de selecao.'),
            ('Se adicionar/remover genes?',
             'Pareto mostra plateau apos ~15 genes. Gap de 0.009 para 50 genes. Leave-one-gene-out: max |DeltaC| = 0.044. Score renormalizado pelo denominador.'),
        ]),
        ('Estatistica', [
            ('I2 de 38.2% nao e problema?',
             'Moderado (Higgins 2003). Fontes estruturais: DSS vs OS, follow-up 32 vs 159 meses. OS-harmonised confirma. Direcao consistente 4/4.'),
            ('C-index de 0.62-0.66 e bom?',
             'Prosigna/ROR atinge 0.62-0.70. Valor incremental (+0.030 a +0.094 sobre CORE-A) e mais informativo. METABRIC ER+: +0.032 sobre 6 covariaveis.'),
            ('Calibration slopes de 1.57-2.74?',
             'Sub-confianca: subestima riscos extremos. Esperado cross-platform sem recalibracao. Discriminacao preservada. Analogia: termometro preciso em ranking mas precisa recalibracao.'),
        ]),
        ('Biologia', [
            ('CorePAM so funciona em ER+?',
             'Coorte de derivacao inclui todos subtipos. ER+: HR 1.59-1.92. ER-: SCAN-B HR=2.53 (instavel, 35 eventos), METABRIC HR=0.97 (186 eventos). Consistente com OncotypeDX/Prosigna. Genes retidos sao luminais e de proliferacao.'),
            ('Pesos contra-intuitivos fazem sentido?',
             'Sim, em contexto multivariado. EXO1/PTTG1 capturam proliferacao-risco; MYBL2 captura componente RESIDUAL. ESR1 positivo captura expressao alta com falha hormonal.'),
        ]),
        ('Clinica', [
            ('Vantagem pratica de 24 vs 50 genes?',
             'RT-qPCR multiplex pode medir 24 genes mais barato que NanoString. RNA de FFPE pode ser suficiente para 24 alvos. MAS: artigo demonstra viabilidade da reducao. Validacao analitica fora do escopo.'),
            ('Muda a conduta clinica?',
             'Nao ainda. Necessario: validacao prospectiva, validacao analitica, comparacao direta com Prosigna. CorePAM e candidato para esse pipeline.'),
            ('SCAN-B (Suecia) generaliza para o Brasil?',
             'Validacao inclui EUA, UK/Canada, Taiwan, Suecia. Diversidade de plataformas sugere robustez, mas coorte brasileira seria proximo passo.'),
        ]),
        ('Reproducibilidade e Etica', [
            ('Resultados sao reprodutiveis?',
             'Dados publicos (GEO, GDC, cBioPortal). Codigo GitHub. SHA-256 folds. analysis_freeze.csv. Anti-hardcoding assertions.'),
            ('IA foi usada - como garantir sem vies?',
             'Claude como assistente computacional (organizacao, scripts). NENHUMA decisao analitica delegada. Codigo revisado pelo autor. IA nao viu os dados.'),
            ('Etica - dados publicos precisam de aprovacao?',
             'Nao. Dados publicos de-identificados. Nenhum dado novo de humanos coletado.'),
        ]),
    ]

    for section_name, qa_list in all_qa:
        pdf.section_title(section_name, 3)
        for q, a in qa_list:
            pdf.set_font('Arial', 'B', 10)
            pdf.set_text_color(26, 58, 92)
            pdf.set_x(pdf.l_margin)
            pdf.multi_cell(0, 5, 'P: "' + q + '"')
            pdf.set_font('Arial', '', 10)
            pdf.set_text_color(30, 30, 30)
            pdf.set_x(pdf.l_margin)
            pdf.multi_cell(0, 5, 'R: ' + a)
            pdf.ln(3)

    # ── 10. CONEXAO COM A TESE ──
    pdf.add_page()
    pdf.section_title('10. Conexao com a Tese de Doutorado', 2)

    pdf.section_title('Contribuicoes Originais', 3)
    pdf.bullet('Framework: primeiro estudo com reducao data-driven do PAM50, criterio NI pre-especificado, SHA-256 folds, Pareto frontier')
    pdf.bullet('Evidencia: 24 genes retêm essencia prognostica, head-to-head SUPERIOR em 3/4 coortes')
    pdf.bullet('Dual endpoint: prognostico (OS/DSS) + predicao de resposta (pCR), I2=0% no pCR')
    pdf.bullet('Transparencia: pipeline open-source, parametros congelados, anti-hardcoding')

    pdf.section_title('Estrutura Possivel da Tese', 3)
    chapters = [
        ['Cap 1', 'Introducao', 'Epidemiologia mama, classificacao molecular, testes genomicos, gap de acesso'],
        ['Cap 2', 'Revisao Literatura', 'Subtipos, TAILORx/MINDACT, metodos selecao genes, cenario SUS'],
        ['Cap 3', 'Objetivos', 'Geral: derivar e validar CorePAM. Especificos: 5 objetivos'],
        ['Cap 4', 'Metodos', 'Expandido do artigo: derivacao matematica, frozen z-score, meta-analise'],
        ['Cap 5', 'Resultados', 'Todos do artigo + suplementares em detalhe'],
        ['Cap 6', 'Discussao', 'Cenario brasileiro, comparacao Prosigna/OncotypeDX, custos SUS'],
        ['Cap 7', 'Conclusoes', 'Perspectivas: coorte brasileira, RT-qPCR, custo-efetividade'],
    ]
    pdf.add_simple_table(['Cap', 'Titulo', 'Conteudo'], chapters)

    pdf.section_title('Expandir alem do artigo', 3)
    pdf.bullet('Contexto SUS: Portaria 874/2013, CONITEC, acesso a testes genomicos')
    pdf.bullet('Custo: estimativa RT-qPCR 24 genes vs Prosigna (~$3,500 USD)')
    pdf.bullet('Epidemiologia brasileira: INCA, distribuicao de subtipos no Brasil')
    pdf.bullet('Framework regulatorio: ANVISA, CLIA, padronizacao de ensaios moleculares')

    # ── 11. GLOSSARIO ──
    pdf.add_page()
    pdf.section_title('11. Glossario Completo', 2)

    glossary = [
        ['C-index', 'Proporcao de pares concordantes (0.5=random, 1.0=perfeito)'],
        ['DeltaC', 'Diferenca entre C-indices (valor incremental)'],
        ['HR', 'Hazard ratio (>1 = maior risco, por 1-SD)'],
        ['OR', 'Odds ratio (usado no bloco pCR)'],
        ['I2', '% variabilidade que e heterogeneidade real'],
        ['OOF', 'Out-of-fold (predicao em dados nao usados para treino)'],
        ['CORE-A', 'Modelo clinico baseline (age +/- ER)'],
        ['DSS', 'Disease-specific survival'],
        ['OS', 'Overall survival'],
        ['pCR', 'Pathologic complete response'],
        ['Elastic-net', 'Regularizacao L1+L2 (selecao + grouping)'],
        ['LASSO', 'L1 regularization (selecao agressiva)'],
        ['Alpha', 'Mixing parameter (0=Ridge, 1=LASSO)'],
        ['Lambda', 'Regularization strength'],
        ['df', 'Degrees of freedom (genes com peso != 0)'],
        ['REML', 'Restricted Maximum Likelihood'],
        ['LRT', 'Likelihood Ratio Test'],
        ['IPCW', 'Inverse Probability of Censoring Weights'],
        ['DCA', 'Decision Curve Analysis'],
        ['TMM', 'Trimmed Mean of M-values (RNA-seq)'],
        ['SHA-256', 'Hash para fold deterministico'],
        ['TRIPOD', 'Guideline para modelos prognosticos'],
        ['ROR', 'Risk of Recurrence (PAM50-based)'],
        ['FFPE', 'Formalin-Fixed Paraffin-Embedded'],
    ]
    pdf.add_simple_table(['Termo', 'Definicao'], glossary)

    return pdf


def main():
    print("Building study guide PDF with all figures...")

    # Check for missing figures
    missing = []
    for key, info in FIGURE_BLOCKS.items():
        if not os.path.exists(info['path']):
            missing.append(f"  {key}: {info['path']}")
    if missing:
        print(f"WARNING: {len(missing)} figures not found:")
        for m in missing:
            print(m)

    pdf = build_pdf()

    # Save to submission/
    pdf_path = os.path.join(BASE, 'GUIA_ESTUDO_CorePAM_v3.pdf')
    pdf.output(pdf_path)
    print(f"PDF saved: {pdf_path}")

    # Also save to Desktop
    desktop = os.path.join(os.path.expanduser('~'), 'Desktop', 'GUIA_ESTUDO_CorePAM_v3.pdf')
    pdf.output(desktop)
    print(f"PDF also saved: {desktop}")

    return pdf_path


if __name__ == '__main__':
    main()
