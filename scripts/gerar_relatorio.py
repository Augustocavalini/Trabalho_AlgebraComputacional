"""Gera relatorio/relatorio.pdf a partir dos resultados em results/ e das figuras em figures/.

As tabelas e os números citados no texto são lidos dos CSVs, para que o relatório
acompanhe automaticamente uma nova rodada dos experimentos.

Uso: conda run -n tp-alc python scripts/gerar_relatorio.py
"""
import csv
import re
import sys
from pathlib import Path

import matplotlib
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import cm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    Image,
    KeepTogether,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"

U = 1.1e-16

# As fontes DejaVu que acompanham o matplotlib têm os símbolos usados no texto, como κ, ρ e ∞.
PASTA_FONTES = Path(matplotlib.get_data_path()) / "fonts" / "ttf"
pdfmetrics.registerFont(TTFont("Serif", PASTA_FONTES / "DejaVuSerif.ttf"))
pdfmetrics.registerFont(TTFont("Serif-Bold", PASTA_FONTES / "DejaVuSerif-Bold.ttf"))
pdfmetrics.registerFont(TTFont("Serif-Italic", PASTA_FONTES / "DejaVuSerif-Italic.ttf"))
pdfmetrics.registerFont(TTFont("Serif-BoldItalic", PASTA_FONTES / "DejaVuSerif-BoldItalic.ttf"))
pdfmetrics.registerFont(TTFont("Mono", PASTA_FONTES / "DejaVuSansMono.ttf"))
pdfmetrics.registerFontFamily("Serif", normal="Serif", bold="Serif-Bold", italic="Serif-Italic", boldItalic="Serif-BoldItalic")

CORPO = ParagraphStyle("corpo", fontName="Serif", fontSize=9.8, leading=14.6, alignment=TA_JUSTIFY, spaceAfter=7)
TITULO = ParagraphStyle("titulo", fontName="Serif-Bold", fontSize=13, leading=16, alignment=TA_CENTER, spaceAfter=2)
SUBTITULO = ParagraphStyle("subtitulo", fontName="Serif", fontSize=8.8, leading=11.5, alignment=TA_CENTER, spaceAfter=2)
SECAO = ParagraphStyle("secao", fontName="Serif-Bold", fontSize=11.5, leading=14.5, spaceBefore=15, spaceAfter=6, keepWithNext=1)
SUBSECAO = ParagraphStyle("subsecao", fontName="Serif-Bold", fontSize=10, leading=13, spaceBefore=10, spaceAfter=5, keepWithNext=1)
LEGENDA = ParagraphStyle("legenda", fontName="Serif-Italic", fontSize=8.2, leading=10.8, alignment=TA_CENTER, spaceBefore=6, spaceAfter=13)
LEGENDA_TABELA = ParagraphStyle("legenda_tabela", parent=LEGENDA, spaceBefore=9, spaceAfter=4, keepWithNext=1)
CELULA = ParagraphStyle("celula", fontName="Serif", fontSize=7.1, leading=9.4, alignment=TA_CENTER)
CELULA_ESQ = ParagraphStyle("celula_esq", parent=CELULA, alignment=0)
CABECALHO = ParagraphStyle("cabecalho", parent=CELULA, fontName="Serif-Bold")


def ler_csv(nome):
    with open(RESULTS / nome, newline="") as arquivo:
        return list(csv.DictReader(arquivo))


def num(linha, campo):
    return float(linha[campo])


def fixo(valor, casas=2):
    """Número com casas fixas e vírgula decimal."""
    return f"{valor:.{casas}f}".replace(".", ",")


def inteiro(valor):
    """Inteiro com ponto separando os milhares."""
    return f"{int(valor):,}".replace(",", ".")


def sci(valor, casas=1):
    """Notação científica para os textos e as tabelas."""
    if valor == 0:
        return "0"
    mantissa, expoente = f"{valor:.{casas}e}".split("e")
    return f"{mantissa.replace('.', ',')}×10<super>{int(expoente)}</super>"


def p(texto, estilo=CORPO):
    return Paragraph(texto, estilo)


LARGURA_UTIL = 17.2 * cm  # largura da página menos as margens
PADDING = 2.5


def largura_do_texto(texto, fonte, tamanho):
    """Largura de um texto formatado. Os expoentes, em <super>, são menores que o resto."""
    largura = 0.0
    partes = re.split(r"(<super>.*?</super>|<sub>.*?</sub>)", texto)
    for parte in partes:
        reduzido = parte.startswith("<super>") or parte.startswith("<sub>")
        simples = re.sub(r"<[^>]+>", "", parte).replace("&lt;", "<").replace("&gt;", ">")
        largura += pdfmetrics.stringWidth(simples, fonte, tamanho * (0.7 if reduzido else 1.0))
    return largura


def larguras_automaticas(cabecalho, linhas, tamanho):
    """Cada coluna recebe a largura do seu maior valor. O cabeçalho pode quebrar, mas em no máximo
    duas linhas: a coluna tem pelo menos a largura da maior palavra e metade da largura do título."""
    larguras = []
    for k, titulo in enumerate(cabecalho):
        palavras = re.split(r"\s+|<br/>", titulo)
        maior_palavra = max(largura_do_texto(w, "Serif-Bold", tamanho) for w in palavras if w)
        meio_titulo = largura_do_texto(titulo.replace("<br/>", " "), "Serif-Bold", tamanho) / 2
        necessaria = max(maior_palavra, meio_titulo)
        for linha in linhas:
            necessaria = max(necessaria, largura_do_texto(str(linha[k]), "Serif", tamanho))
        larguras.append(necessaria + 2 * PADDING + 2)
    return larguras


def tabela(cabecalho, linhas, larguras=None, alinhar_esquerda=(0,), separadores=()):
    """Monta uma tabela com texto formatado. separadores são as linhas com um fio acima.

    As larguras das colunas são calculadas a partir do conteúdo (o argumento larguras é ignorado).
    Se a tabela não couber na página, a letra diminui; se ainda assim não couber, as colunas
    são reduzidas na mesma proporção."""
    tamanho = 7.1
    colunas = larguras_automaticas(cabecalho, linhas, tamanho)
    while sum(colunas) > LARGURA_UTIL and tamanho > 6.2:
        tamanho -= 0.3
        colunas = larguras_automaticas(cabecalho, linhas, tamanho)
    # Todas as tabelas ocupam a largura da página: a sobra (ou a falta) é dividida entre as colunas
    # na proporção das suas larguras.
    fator = LARGURA_UTIL / sum(colunas)
    colunas = [c * fator for c in colunas]

    celula = ParagraphStyle("celula_t", parent=CELULA, fontSize=tamanho, leading=tamanho * 1.32)
    celula_esq = ParagraphStyle("celula_esq_t", parent=celula, alignment=0)
    cabecalho_estilo = ParagraphStyle("cabecalho_t", parent=celula, fontName="Serif-Bold")

    dados = [[Paragraph(c, cabecalho_estilo) for c in cabecalho]]
    for linha in linhas:
        celulas = []
        for k, valor in enumerate(linha):
            estilo = celula_esq if k in alinhar_esquerda else celula
            celulas.append(Paragraph(str(valor), estilo))
        dados.append(celulas)

    t = Table(dados, colWidths=colunas, repeatRows=1)
    estilo = [
        ("LINEABOVE", (0, 0), (-1, 0), 0.8, colors.black),
        ("LINEBELOW", (0, 0), (-1, 0), 0.5, colors.black),
        ("LINEBELOW", (0, -1), (-1, -1), 0.8, colors.black),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 2.4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 2.4),
        ("LEFTPADDING", (0, 0), (-1, -1), PADDING),
        ("RIGHTPADDING", (0, 0), (-1, -1), PADDING),
    ]
    for indice in separadores:
        estilo.append(("LINEABOVE", (0, indice + 1), (-1, indice + 1), 0.3, colors.grey))
    t.setStyle(TableStyle(estilo))
    t.spaceAfter = 11
    return t


def figura(nome, largura_cm):
    caminho = FIGURES / nome
    from reportlab.lib.utils import ImageReader

    largura_px, altura_px = ImageReader(str(caminho)).getSize()
    imagem = Image(str(caminho), width=largura_cm * cm, height=largura_cm * cm * altura_px / largura_px)
    imagem.spaceBefore = 7
    imagem.spaceAfter = 3
    return imagem


def rodape(canvas, documento):
    canvas.saveState()
    canvas.setFont("Serif", 7.5)
    canvas.drawCentredString(A4[0] / 2, 1.0 * cm, str(documento.page))
    canvas.restoreState()


def ler_validacao():
    """Lê a saída de testar_metodos gravada em results/parte2_validacao.txt."""
    padrao = re.compile(r"^(.+?)\s+(OK|VERIFICAR)\s+retroação = (\S+)\s+erro direto = (\S+)\s+diferença para a LAPACK = (\S+)")
    linhas = []
    with open(RESULTS / "parte2_validacao.txt") as arquivo:
        for texto in arquivo:
            encontrado = padrao.match(texto)
            if encontrado:
                linhas.append(encontrado.groups())
    if not linhas:
        sys.exit("results/parte2_validacao.txt não tem resultados")
    return linhas


def main():
    parte1 = {l["matriz"]: l for l in ler_csv("parte1.csv")}
    bandas = {(l["matriz"], l["ordem"]): l for l in ler_csv("parte3.csv")}
    tempos3 = {(l["matriz"], l["metodo"], l["ordem"]): l for l in ler_csv("parte3_tempos.csv")}
    parte4 = ler_csv("parte4_matrizes.csv")
    escala = ler_csv("parte4_escala.csv")
    inclinacoes = ler_csv("parte4_inclinacoes.csv")
    validacao = ler_validacao()

    def r4(matriz, metodo, implementacao):
        for linha in parte4:
            if (linha["matriz"], linha["metodo"], linha["implementacao"]) == (matriz, metodo, implementacao):
                return linha
        raise KeyError((matriz, metodo, implementacao))

    def e4(metodo, implementacao, m):
        for linha in escala:
            if (linha["metodo"], linha["implementacao"], int(linha["m"])) == (metodo, implementacao, m):
                return linha
        raise KeyError((metodo, implementacao, m))

    def razoes(metodo):
        return [num(l, "implementado_sobre_lapack") for l in parte4 + escala if l["metodo"] == metodo and l["implementacao"] == "implementado"]

    def faixa(valores, casas=1):
        return f"{fixo(min(valores), casas)} a {fixo(max(valores), casas)}"

    def t3(matriz, metodo, ordem):
        return tempos3[(matriz, metodo, ordem)]

    def ganho_medido(matriz, metodo, ordem):
        l = t3(matriz, metodo, ordem)
        return num(l, "fatoracao_densa_s") / num(l, "fatoracao_banda_s")

    def ganho_formula(matriz, metodo, ordem):
        l = t3(matriz, metodo, ordem)
        return num(l, "flops_densa_formula") / num(l, "flops_banda_formula")

    def ganho_exato(matriz, metodo, ordem):
        l = t3(matriz, metodo, ordem)
        return num(l, "flops_densa_exata") / num(l, "flops_banda_exata")

    def taxa(matriz, metodo, ordem, forma):
        l = t3(matriz, metodo, ordem)
        return num(l, f"flops_{forma}_exata") / num(l, f"fatoracao_{forma}_s") / 1e9

    # Confere as afirmações sobre os erros que aparecem no texto da questão 3.
    medicoes = [l for l in parte4 + escala if not l["falha"]]
    for linha in medicoes:
        kappa = num(linha, "limite_kappa_u") / U
        if num(linha, "direto") > kappa * num(linha, "retroacao"):
            sys.exit(f"o texto da questão 3 não vale mais para {linha['matriz']}, {linha['metodo']}")
    acima_de_kappa_u = {l["matriz"] for l in medicoes if num(l, "direto") > num(l, "limite_kappa_u")}
    if acima_de_kappa_u != {"densa SPD aleatória"}:
        sys.exit(f"o texto da questão 3 precisa de revisão: {acima_de_kappa_u}")

    conteudo = []
    adicionar = conteudo.append

    # As tabelas e as figuras são numeradas na ordem em que aparecem.
    contadores = {"Tabela": 0, "Figura": 0}

    def rotulo(tipo):
        contadores[tipo] += 1
        return f"{tipo} {contadores[tipo]}."

    # Título
    adicionar(p("Universidade Federal de Juiz de Fora", SUBTITULO))
    adicionar(p("Programa de Pós-graduação em Modelagem Computacional · Álgebra Linear Computacional · 2026.3", SUBTITULO))
    adicionar(Spacer(1, 3))
    adicionar(p("Métodos diretos e comparação com a LAPACK", TITULO))
    adicionar(p("Augusto De Filippo Cavalini · Prof. José J. Camata", SUBTITULO))
    adicionar(Spacer(1, 4))

    # 1. Matrizes
    b15 = parte1["bcsstk15"]
    bus = parte1["1138_bus"]
    st2 = parte1["steam2"]
    adicionar(p("1. Matrizes escolhidas", SECAO))
    adicionar(p(
        "As três matrizes vêm da coleção Harwell-Boeing, obtida no repositório Matrix Market do NIST, porque os downloads da "
        "SuiteSparse devolviam erro de acesso na época dos testes. A <b>bcsstk15</b> é a matriz de rigidez do módulo de uma plataforma "
        "offshore (análise estática estrutural). É a maior das três e representa o caso (a), SPD. A <b>steam2</b> vem de um modelo "
        "3D de injeção de vapor em reservatório de petróleo (grade 5×5×6 com 4 graus de liberdade por célula). É não simétrica e "
        "cobre o caso (b). Também cobre o caso (c), pois o RCM reduz a sua banda de 331 para 76. A <b>1138_bus</b> é a matriz de "
        "admitância de uma rede elétrica de 1138 barras, SPD. Ela dá um segundo exemplo do caso (c), agora resolvido com Cholesky. "
        "Todas têm 500 ≤ m ≤ 5000, e as origens diferentes dão padrões de esparsidade e condicionamentos variados."
    ))

    # 2. Metodologia
    adicionar(p("2. Máquina e metodologia", SECAO))
    adicionar(p(
        "<b>Máquina:</b> AMD Ryzen 5 5600 (6 núcleos; caches L1d de 32 KiB e L2 de 512 KiB por núcleo, L3 de 32 MiB compartilhada), "
        "15 GiB de memória, Ubuntu 26.04. <b>Software:</b> g++ 15.2.0 com <font name='Mono'>-O2 -march=native -std=c++17</font>, "
        "OpenBLAS 0.3.32, LAPACKE 3.12.1 e Boost 1.90. "
        "<b>Medição:</b> <font name='Mono'>OPENBLAS_NUM_THREADS=1</font>; <font name='Mono'>std::chrono::steady_clock</font>; mediana "
        "de 5 execuções, alternando a versão implementada e a da LAPACK a cada repetição; cópias de A fora do cronômetro. Em todos os "
        "testes, x = (1, ..., 1)<super>T</super> e b = Ax. A verificação de NaN do LAPACKE foi desligada "
        "(<font name='Mono'>LAPACKE_set_nancheck(0)</font>), porque ela percorre os dados antes de cada chamada e não faz parte das "
        "rotinas da LAPACK. Como os erros são medidos na norma infinito, o limite do erro usa κ<sub>∞</sub>(A), estimado "
        "com <font name='Mono'>dlange</font> + <font name='Mono'>dgetrf</font> + <font name='Mono'>dgecon</font>; u = 1,1×10<super>-16</super>."
    ))
    adicionar(p(
        "<b>Resíduo e erro.</b> O <b>resíduo</b> usado em todo o relatório é o resíduo relativo ‖b − Ax̂‖<sub>∞</sub> / (‖A‖<sub>∞</sub> ‖x̂‖<sub>∞</sub>). "
        "Por um resultado de Rigal e Gaches, ele é exatamente a menor perturbação relativa de A para a qual a solução calculada x̂ seria a solução exata "
        "(o <i>backward error</i>). Um resíduo próximo de u significa que o método resolveu exatamente um sistema muito próximo do original, isto é, que ele é estável. "
        "O <b>erro</b> é ‖x̂ − x‖<sub>∞</sub> / ‖x‖<sub>∞</sub>, a distância até a solução verdadeira. Ele depende também do condicionamento: "
        "em primeira ordem, erro ≤ κ(A) × resíduo."
    ))

    # 3. Implementações
    adicionar(p("3. Implementações", SECAO))
    adicionar(p(
        "<b>Armazenamento denso:</b> um único vetor em ordem coluna-maior, a(i, j) = a[i + j·m], o formato que a LAPACK espera. "
        "<b>Banda:</b> o arranjo da LAPACK, com (b<sub>I</sub> + b<sub>S</sub> + 1)·m posições, em que a(i, j) fica na linha "
        "b<sub>S</sub> + i − j da coluna j. No Cholesky em banda, b<sub>S</sub> = 0 e só o triângulo inferior é guardado, em (b + 1)·m "
        "posições, como no <font name='Mono'>dpbtrf</font>. <b>Thomas:</b> três vetores (inferior, diagonal, superior). "
        "Todas as fatorações sobrescrevem a matriz recebida e ficam separadas das substituições."
    ))
    adicionar(p(
        "<b>Decisões:</b> (i) a LU e o Cholesky usam a forma right-looking, com o laço mais interno ao longo das colunas, a ordem da "
        "memória. Em uma medição pontual com m = 2000, fora dos scripts deste relatório, a troca da ordem por linhas pela ordem por "
        "colunas reduziu a LU de 6,8 s para 1,2 s, sem mudar nenhum coeficiente. (ii) O vetor de pivôs usa a convenção da LAPACK, a partir "
        "de 1. (iii) O teste de pivô do Cholesky é <font name='Mono'>!(a<sub>kk</sub> &gt; 0)</font>, que também rejeita NaN. "
        "(iv) A LU em banda e o Thomas param em pivô nulo, mas não têm defesa contra pivôs pequenos; o resíduo é o que revela esse caso."
    ))
    linhas_validacao = [[nome, situacao, sci(float(valor_r)), sci(float(valor_d)), sci(float(valor_l))] for nome, situacao, valor_r, valor_d, valor_l in validacao]
    adicionar(KeepTogether([
        p(f"{rotulo('Tabela')} Validação em uma matriz pequena (tridiagonal SPD, m = 10), que todos os métodos conseguem resolver.", LEGENDA_TABELA),
        tabela(["Método", "Situação", "Resíduo", "Erro", "Diferença para a solução da LAPACK"], linhas_validacao, [4, 2, 2.6, 2.6, 4.8]),
    ]))

    # 4. Parte 1
    linhas = []
    for nome in ("bcsstk15", "steam2", "1138_bus"):
        l = parte1[nome]
        linhas.append([
            nome, inteiro(l["m"]), inteiro(l["nnz"]), sci(num(l, "densidade"), 2),
            "sim" if l["simetrica"] == "1" else "não", "sim" if l["spd"] == "1" else "não",
            l["bI"], l["bS"], sci(num(l, "norma1"), 2), sci(num(l, "kappa1"), 2), l["rotina"],
        ])
    adicionar(KeepTogether([
        p("4. Parte 1: caracterização", SECAO),
        p(f"{rotulo('Tabela')} Caracterização das matrizes.", LEGENDA_TABELA),
        tabela(["Matriz", "m", "nnz", "nnz/m²", "Simétrica", "SPD", "b<sub>I</sub>", "b<sub>S</sub>", "‖A‖<sub>1</sub>", "κ<sub>1</sub> estimado", "Rotina"],
               linhas, [1.7, 1.05, 1.25, 1.9, 1.65, 0.9, 1, 1, 2, 2, 1.5]),
    ]))
    adicionar(p(
        f"O arquivo da steam2 tem 13.760 entradas, mas 8.100 são zeros guardados explicitamente; por isso nnz = {inteiro(st2['nnz'])}. "
        "As bandas e o grafo do RCM consideram só os valores realmente não nulos, como pede o enunciado. A página do NIST informa "
        "banda 438 para a bcsstk15 porque conta a diagonal; aqui b = 437 não a conta. Os valores de κ são da mesma ordem dos "
        f"informados pelas fontes: {sci(num(b15, 'kappa1'), 1)} contra 6,5×10<super>9</super> na bcsstk15, {sci(num(bus, 'kappa1'), 1)} contra "
        f"8,6×10<super>6</super> na 1138_bus e {sci(num(st2, 'kappa1'), 1)} contra 3,5×10<super>6</super> na steam2."
    ))
    painel = Table([[figura("spy_bcsstk15.png", 5.5), figura("spy_steam2.png", 5.5), figura("spy_1138_bus.png", 5.5)]], colWidths=[5.7 * cm] * 3)
    painel.setStyle(TableStyle([("LEFTPADDING", (0, 0), (-1, -1), 0), ("RIGHTPADDING", (0, 0), (-1, -1), 0)]))
    adicionar(KeepTogether([painel, p(f"{rotulo('Figura')} Padrões de esparsidade na ordem original.", LEGENDA)]))

    # 5. Parte 3
    adicionar(p("5. Parte 3: reordenamento RCM e fatoração em banda", SECAO))
    adicionar(p(
        "A permutação vem de <font name='Mono'>boost::cuthill_mckee_ordering</font>, com iterador reverso, sobre o grafo não direcionado "
        "de A. Com ela, formei Ã = PAP<super>T</super> e resolvi Ã(Px) = Pb. O programa confere que ‖Ã(Px) − Pb‖ é da ordem de "
        "10<super>-16</super> em termos relativos e que o nnz não muda. O caso (c) do enunciado é a steam2; a 1138_bus é um segundo "
        "exemplo; a bcsstk15 entra como contraste."
    ))
    linhas = []
    separadores = []
    for k, nome in enumerate(("steam2", "1138_bus", "bcsstk15")):
        if k > 0:
            separadores.append(len(linhas))
        for ordem in ("original", "RCM"):
            l = bandas[(nome, ordem)]
            cholesky = inteiro(l["posicoes_cholesky"]) if l["simetrica"] == "1" else "não simétrica"
            linhas.append([nome, ordem, l["bI"], l["bS"], inteiro(l["posicoes_lu"]), cholesky, inteiro(l["posicoes_densa"])])
    adicionar(KeepTogether([
        p(f"{rotulo('Tabela')} Larguras de banda e armazenamento, em número de posições, antes e depois do RCM.", LEGENDA_TABELA),
        tabela(["Matriz", "Ordem", "b<sub>I</sub>", "b<sub>S</sub>", "LU em banda: (b<sub>I</sub>+b<sub>S</sub>+1)m", "Cholesky em banda: (b+1)m", "Densa: m²"],
               linhas, [2, 1.8, 1.2, 1.2, 3.8, 3.8, 2.4], separadores=separadores),
    ]))
    adicionar(figura("rcm_steam2.png", 11.5))
    adicionar(KeepTogether([
        figura("rcm_1138_bus.png", 11.5),
        figura("rcm_bcsstk15.png", 11.5),
        p(f"{rotulo('Figura')} Padrões antes e depois do RCM. Na steam2 (caso c) e na 1138_bus, a banda cai bastante; na bcsstk15, ela aumenta.", LEGENDA),
    ]))

    # Linhas das tabelas de comparação: LU em todas as matrizes e Cholesky nas SPD, sem e com RCM.
    casos = []
    for nome in ("steam2", "1138_bus", "bcsstk15"):
        for metodo in ("LU", "Cholesky"):
            for ordem in ("original", "RCM"):
                if (nome, metodo, ordem) in tempos3:
                    casos.append((nome, metodo, ordem))
    falhas = [c for c in casos if t3(*c)["falha_banda"]]
    if falhas:
        sys.exit(f"a LU em banda falhou em {falhas}: o texto da Parte 3 precisa de revisão")

    separadores_casos = []
    for k in range(1, len(casos)):
        if casos[k][0] != casos[k - 1][0]:
            separadores_casos.append(k)

    def nome_ordem(ordem):
        return "sem RCM" if ordem == "original" else "com RCM"

    primeira_tabela = contadores["Tabela"] + 1
    adicionar(p(
        "Para cada matriz, resolvi Ax = b com LU, que vale para todas, e com Cholesky, que vale para as duas SPD. Cada método foi "
        "executado nas formas densa e em banda e nas duas ordens: <b>sem RCM</b>, com a matriz A na numeração original, e "
        "<b>com RCM</b>, com Ã = PAP<super>T</super>. <b>Densa</b> é a fatoração com a matriz guardada inteira, em m² posições, "
        "que percorre toda a submatriz restante em cada etapa. <b>Em banda</b> é a mesma fatoração com o armazenamento da Seção 3, "
        "que só percorre a faixa entre b<sub>I</sub> e b<sub>S</sub>. Na LU, a densa usa pivoteamento parcial e a em banda não. "
        f"Todas as fatorações das Tabelas {primeira_tabela} a {primeira_tabela + 2} são as implementadas neste trabalho; "
        "a LAPACK só aparece na Parte 4. Os erros são sempre medidos no sistema original Ax = b: na ordem com RCM, x é recuperado "
        "de Px antes do cálculo."
    ))

    linhas = []
    for nome, metodo, ordem in casos:
        l = t3(nome, metodo, ordem)
        linhas.append([
            nome, metodo, nome_ordem(ordem), l["bI"], l["bS"],
            sci(num(l, "fatoracao_densa_s"), 2), sci(num(l, "fatoracao_banda_s"), 2),
            sci(num(l, "solucao_densa_s"), 2), sci(num(l, "solucao_banda_s"), 2),
        ])
    adicionar(KeepTogether([
        p(f"{rotulo('Tabela')} Tempos de fatoração e de solução, em segundos (mediana de 5 execuções).", LEGENDA_TABELA),
        tabela(["Matriz", "Método", "Ordem", "b<sub>I</sub>", "b<sub>S</sub>", "Fatoração densa (s)", "Fatoração em banda (s)", "Solução densa (s)", "Solução em banda (s)"],
               linhas, [1.6, 1.5, 1.6, 1.1, 1.1, 2.4, 2.4, 2.4, 2.4], alinhar_esquerda=(0, 1, 2), separadores=separadores_casos),
    ]))

    linhas = []
    for caso in casos:
        nome, metodo, ordem = caso
        linhas.append([
            nome, metodo, nome_ordem(ordem), fixo(ganho_medido(*caso), 2), fixo(ganho_formula(*caso), 2),
            fixo(taxa(*caso, "densa")), fixo(taxa(*caso, "banda")),
        ])
    adicionar(KeepTogether([
        p(f"{rotulo('Tabela')} Ganho da fatoração em banda sobre a densa (quantas vezes ela é mais rápida) e taxas das duas formas. O ganho medido é a razão entre os tempos; "
          "o previsto é a razão entre as contagens de operações da fórmula do enunciado: (2/3)m³ ÷ 2m·b<sub>I</sub>·b<sub>S</sub> na LU e "
          "(1/3)m³ ÷ m·b² no Cholesky. As taxas em GFLOP/s dividem pelo tempo a contagem exata das operações que os laços executam.", LEGENDA_TABELA),
        tabela(["Matriz", "Método", "Ordem", "Ganho medido", "Ganho previsto", "GFLOP/s densa", "GFLOP/s em banda"],
               linhas, alinhar_esquerda=(0, 1, 2), separadores=separadores_casos),
    ]))

    linhas = []
    for nome, metodo, ordem in casos:
        l = t3(nome, metodo, ordem)
        linhas.append([
            nome, metodo, nome_ordem(ordem),
            sci(num(l, "retroacao_densa")), sci(num(l, "retroacao_banda")), sci(num(l, "direto_densa")), sci(num(l, "direto_banda")),
        ])
    adicionar(KeepTogether([
        p(f"{rotulo('Tabela')} Erros no sistema original Ax = b (definições na Seção 2).", LEGENDA_TABELA),
        tabela(["Matriz", "Método", "Ordem", "Resíduo: densa", "Resíduo: em banda", "Erro: densa", "Erro: em banda"],
               linhas, [1.8, 1.6, 1.7, 2.6, 2.6, 2.5, 2.5], alinhar_esquerda=(0, 1, 2), separadores=separadores_casos),
    ]))

    taxa_densa_b15 = taxa("bcsstk15", "Cholesky", "original", "densa")
    taxa_banda_b15 = taxa("bcsstk15", "Cholesky", "original", "banda")
    adicionar(p(
        f"<b>Efeito do RCM.</b> Na steam2, a banda em ordem original já ajuda pouco (ganho de {fixo(ganho_medido('steam2', 'LU', 'original'))}); com RCM, o ganho sobe "
        f"para {fixo(ganho_medido('steam2', 'LU', 'RCM'), 1)}. Na 1138_bus, sem RCM a banda ocupa quase toda a matriz (b/m = 0,91) e a forma em banda leva o mesmo tempo da densa "
        f"(ganho de {fixo(ganho_medido('1138_bus', 'LU', 'original'))} na LU e {fixo(ganho_medido('1138_bus', 'Cholesky', 'original'))} no Cholesky); com RCM, ele vai a "
        f"{fixo(ganho_medido('1138_bus', 'LU', 'RCM'), 1)} e {fixo(ganho_medido('1138_bus', 'Cholesky', 'RCM'), 1)}. Na bcsstk15, a numeração original já tem banda estreita, e o RCM "
        f"a piora: o ganho cai de {fixo(ganho_medido('bcsstk15', 'LU', 'original'), 1)} para {fixo(ganho_medido('bcsstk15', 'LU', 'RCM'), 1)} na LU e de "
        f"{fixo(ganho_medido('bcsstk15', 'Cholesky', 'original'), 1)} para {fixo(ganho_medido('bcsstk15', 'Cholesky', 'RCM'), 1)} no Cholesky. A forma densa leva praticamente o "
        "mesmo tempo nas duas ordens, porque percorre a matriz inteira de qualquer forma. Nas matrizes SPD, o Cholesky em banda é cerca de duas vezes mais rápido que a LU em banda, "
        f"como prevê a razão 2m·b² ÷ m·b² = 2 ({sci(num(t3('1138_bus', 'LU', 'RCM'), 'fatoracao_banda_s'), 2)} s contra {sci(num(t3('1138_bus', 'Cholesky', 'RCM'), 'fatoracao_banda_s'), 2)} s na 1138_bus com RCM)."
    ))
    adicionar(p(
        "<b>Medido e previsto.</b> A fórmula do enunciado só vale para b muito menor que m. Quando b se aproxima de m, ela erra bastante: na 1138_bus sem RCM, prevê "
        f"{fixo(ganho_formula('1138_bus', 'LU', 'original'))} (a banda seria mais lenta), mas o ganho medido é {fixo(ganho_medido('1138_bus', 'LU', 'original'))}, porque os laços em banda "
        "nunca passam do tamanho da matriz. A contagem exata corrige isso e acompanha o ganho medido em até cerca de 10% na steam2 e na 1138_bus, onde as duas formas "
        "têm a mesma taxa, perto de 5 GFLOP/s. Na bcsstk15, o ganho medido supera até a contagem exata: a forma densa ocupa 125 MB, bem mais que os 32 MB da L3, e sua taxa "
        f"cai para {fixo(taxa_densa_b15)} GFLOP/s, enquanto a em banda, com cerca de 14 MB, mantém {fixo(taxa_banda_b15)}. A diferença de taxas explica o resto: "
        f"{fixo(ganho_exato('bcsstk15', 'Cholesky', 'original'), 1)} × {fixo(taxa_banda_b15)} ÷ {fixo(taxa_densa_b15)} ≈ {fixo(ganho_exato('bcsstk15', 'Cholesky', 'original') * taxa_banda_b15 / taxa_densa_b15, 1)}, "
        f"perto do ganho medido de {fixo(ganho_medido('bcsstk15', 'Cholesky', 'original'), 1)}."
    ))
    adicionar(p(
        "<b>Precisão.</b> As formas densa e em banda têm resíduos da mesma ordem, entre 10<super>-16</super> e 10<super>-15</super>, em todos os casos, "
        "inclusive nas LU em banda, que não pivotam, e o RCM não altera esse patamar. O erro varia um pouco mais entre as versões, por exemplo de "
        f"{sci(num(t3('bcsstk15', 'LU', 'original'), 'direto_densa'))} sem RCM para {sci(num(t3('bcsstk15', 'LU', 'RCM'), 'direto_densa'))} com RCM na LU densa da bcsstk15, "
        f"mas fica sempre abaixo de κ<sub>∞</sub>u ({sci(num(r4('bcsstk15', 'LU densa', 'implementado'), 'limite_kappa_u'))} nessa matriz). "
        f"O RCM custou {sci(num(t3('steam2', 'LU', 'RCM'), 'tempo_rcm_s'))} s na steam2, {sci(num(t3('1138_bus', 'LU', 'RCM'), 'tempo_rcm_s'))} s na 1138_bus e "
        f"{sci(num(t3('bcsstk15', 'LU', 'RCM'), 'tempo_rcm_s'))} s na bcsstk15, incluindo a montagem do grafo a partir da forma densa."
    ))

    # 6. Parte 4
    adicionar(p("6. Parte 4: comparação com a LAPACK", SECAO))
    linhas = []
    separadores = []
    anterior = None
    for l in parte4:
        if anterior is not None and l["matriz"] != anterior:
            separadores.append(len(linhas))
        anterior = l["matriz"]
        razao = fixo(num(l, "implementado_sobre_lapack")) if l["implementacao"] == "implementado" else "-"
        rho = fixo(num(l, "rho")) if num(l, "rho") >= 0 else "-"
        linhas.append([
            l["matriz"], l["metodo"], l["implementacao"], sci(num(l, "fatoracao_s"), 2), sci(num(l, "solucao_s"), 2), fixo(num(l, "gflops")),
            razao, sci(num(l, "retroacao")), sci(num(l, "direto")), sci(num(l, "limite_kappa_u")), rho,
        ])
    adicionar(p(f"{rotulo('Tabela')} Resultados por matriz e por método, na ordem original das matrizes. Nas linhas de LU, ρ = max|u<sub>ij</sub>| / max|a<sub>ij</sub>|.", LEGENDA_TABELA))
    adicionar(tabela(["Matriz", "Método", "Implementação", "Fatoração<br/>(s)", "Solução<br/>(s)", "GFLOP/s", "Razão<br/>(÷ LAPACK)", "Resíduo", "Erro", "κ<sub>∞</sub>u", "ρ"],
                     linhas, [1.45, 2.45, 2.25, 1.6, 1.6, 1.3, 1.45, 1.5, 1.35, 1.35, 0.75], alinhar_esquerda=(0, 1, 2), separadores=separadores))
    adicionar(Spacer(1, 3))
    adicionar(p(
        "Na forma em banda com ordem original, a taxa em GFLOP/s precisa de cuidado quando b é próximo de m: na 1138_bus (b/m = 0,91), "
        f"a fórmula 2m·b<sub>I</sub>·b<sub>S</sub> vale cerca de 1,6m³, bem acima do custo real, que não passa de (2/3)m³. Por isso a taxa aparente é de {fixo(num(r4('1138_bus', 'LU em banda', 'implementado'), 'gflops'), 1)} GFLOP/s "
        f"na LU em banda implementada, embora o tempo seja o mesmo da densa, com {fixo(num(r4('1138_bus', 'LU densa', 'implementado'), 'gflops'), 1)} GFLOP/s. Além disso, a contagem da tabela não inclui o pivoteamento "
        "do <font name='Mono'>dgbtrf</font>. Sem pivoteamento, ρ pode diferir: na bcsstk15, foi "
        f"{fixo(num(r4('bcsstk15', 'LU em banda', 'implementado'), 'rho'))} na LU em banda implementada e {fixo(num(r4('bcsstk15', 'LU em banda', 'LAPACK'), 'rho'))} no <font name='Mono'>dgbtrf</font>."
    ))

    linhas = []
    for m in (250, 500, 1000, 2000, 4000):
        lu, dg = e4("LU densa", "implementado", m), e4("LU densa", "LAPACK", m)
        ch, dp = e4("Cholesky densa", "implementado", m), e4("Cholesky densa", "LAPACK", m)
        linhas.append([
            inteiro(m), sci(num(lu, "fatoracao_s"), 2), sci(num(dg, "fatoracao_s"), 2), fixo(num(lu, "implementado_sobre_lapack"), 1),
            fixo(num(lu, "gflops")), fixo(num(dg, "gflops")), fixo(num(lu, "rho"), 1),
            sci(num(ch, "fatoracao_s"), 2), sci(num(dp, "fatoracao_s"), 2), fixo(num(ch, "implementado_sobre_lapack"), 1),
            fixo(num(ch, "gflops")), fixo(num(dp, "gflops")),
        ])
    adicionar(KeepTogether([
        p("Estudo de escala", SUBSECAO),
        p(f"{rotulo('Tabela')} Fatoração de matrizes densas aleatórias. A matriz da LU tem coeficientes uniformes em [−1, 1]; a do Cholesky é simétrica, com m somado à diagonal (SPD, κ<sub>∞</sub> ≈ 2,3).", LEGENDA_TABELA),
        tabela(["m", "LU<br/>implementado (s)", "dgetrf (s)", "Razão", "GFLOP/s<br/>implementado", "GFLOP/s dgetrf", "ρ", "Cholesky<br/>implementado (s)", "dpotrf (s)", "Razão", "GFLOP/s<br/>implementado", "GFLOP/s dpotrf"],
               linhas, [0.9, 1.8, 1.4, 1.0, 1.8, 1.35, 0.95, 1.8, 1.4, 1.0, 1.8, 1.35], alinhar_esquerda=()),
    ]))
    linhas = []
    for m in (1000, 10000, 100000, 1000000):
        th, dt = e4("Thomas", "implementado", m), e4("Thomas", "LAPACK", m)
        total = num(th, "fatoracao_s") + num(th, "solucao_s")
        linhas.append([
            inteiro(m), sci(total, 2), sci(num(dt, "fatoracao_s"), 2), fixo(num(th, "implementado_sobre_lapack")),
            fixo(num(th, "gflops")), fixo(num(dt, "gflops")), sci(num(th, "retroacao")), sci(num(th, "direto")), sci(num(th, "limite_kappa_u")),
        ])
    adicionar(KeepTogether([
        p(f"{rotulo('Tabela')} tridiag(−1, 2, −1): Thomas (fatoração + solução) e <font name='Mono'>dgtsv</font>. A taxa usa a contagem de 8m operações.", LEGENDA_TABELA),
        tabela(["m", "Thomas (s)", "dgtsv (s)", "Razão", "GFLOP/s Thomas", "GFLOP/s dgtsv", "Resíduo", "Erro", "κ<sub>∞</sub>u"],
               linhas, [1.6, 2, 2, 1.1, 1.5, 1.5, 2, 2, 2], alinhar_esquerda=()),
    ]))
    adicionar(figura("parte4_escala.png", 17.2))
    adicionar(p(f"{rotulo('Figura')} Tempos em escala log-log. A cor identifica a implementação; o traço e o marcador identificam o método. A reta pontilhada tem a inclinação prevista.", LEGENDA))

    def linha_inclinacao(rotina, etapa):
        for l in inclinacoes:
            if (l["rotina"], l["etapa"]) == (rotina, etapa):
                return l
        raise KeyError((rotina, etapa))

    linhas = []
    for rotina, etapa in (("LU implementado", "fatoração"), ("dgetrf/dgetrs", "fatoração"), ("Cholesky implementado", "fatoração"), ("dpotrf/dpotrs", "fatoração"),
                          ("LU implementado", "solução"), ("dgetrf/dgetrs", "solução"), ("Cholesky implementado", "solução"), ("dpotrf/dpotrs", "solução"),
                          ("Thomas implementado", "fatoração + solução"), ("dgtsv", "fatoração + solução")):
        l = linha_inclinacao(rotina, etapa)
        linhas.append([rotina, etapa, l["prevista"], fixo(float(l["ajustada"])), fixo(float(l["dois_ultimos_pontos"]))])
    adicionar(KeepTogether([
        p(f"{rotulo('Tabela')} Inclinações em log-log: reta ajustada a todos os pontos e inclinação entre os dois últimos pontos.", LEGENDA_TABELA),
        tabela(["Rotina", "Etapa", "Prevista", "Ajustada", "Dois últimos pontos"], linhas, [3.6, 3.6, 1.6, 1.6, 2.8], alinhar_esquerda=(0, 1), separadores=(4, 8)),
    ]))
    lu_final = linha_inclinacao("LU implementado", "fatoração")
    ch_final = linha_inclinacao("Cholesky implementado", "fatoração")
    dg_final = linha_inclinacao("dgetrf/dgetrs", "fatoração")
    adicionar(p(
        "As inclinações concordam com a teoria: cerca de 3 na fatoração densa, 2 na solução e 1 no Thomas. O desvio mais claro é o das fatorações densas "
        f"implementadas entre m = 2000 e 4000 ({fixo(float(lu_final['dois_ultimos_pontos']))} na LU e {fixo(float(ch_final['dois_ultimos_pontos']))} no Cholesky). A taxa cai de cerca de 5 para cerca de 3 GFLOP/s "
        "justamente quando a matriz passa de 32 MB, o tamanho da L3, para 128 MB, enquanto a LAPACK mantém a taxa "
        f"({fixo(float(dg_final['dois_ultimos_pontos']))} no <font name='Mono'>dgetrf</font>). Nas matrizes menores, a LAPACK tem inclinação ajustada abaixo de 3 porque sua taxa ainda está crescendo "
        f"({fixo(num(e4('LU densa', 'LAPACK', 250), 'gflops'), 1)} GFLOP/s em m = 250 e {fixo(num(e4('LU densa', 'LAPACK', 1000), 'gflops'), 1)} em m = 1000). Os tempos de solução são curtos e mais ruidosos, por isso as "
        f"inclinações deles variam mais. Nas densas aleatórias, ρ cresce de {fixo(num(e4('LU densa', 'implementado', 250), 'rho'), 1)} para {fixo(num(e4('LU densa', 'implementado', 4000), 'rho'), 1)}, muito abaixo do pior caso 2<super>m−1</super>."
    ))

    # 7. Discussão
    adicionar(p("7. Questões para discussão", SECAO))
    adicionar(p("<b>1. Por que a LAPACK é mais rápida que a minha LU densa, se os dois códigos fazem o mesmo número de operações?</b>", CORPO))
    adicionar(p(
        f"A contagem de operações não determina o tempo: o acesso à memória e o aproveitamento do processador pesam tanto quanto ela. Minha LU está de "
        f"{faixa(razoes('LU densa'))} vezes atrás, com cerca de 5 GFLOP/s ou menos, enquanto a LAPACK fica entre 50 e 62 GFLOP/s a partir de m = 1000. A minha aplica cada etapa k a toda a submatriz restante, "
        "passando pela matriz inteira m vezes; quando ela não cabe no cache, quase todo acesso custa uma ida à memória. A LAPACK com OpenBLAS organiza a conta em blocos que cabem no cache e usa núcleos "
        "BLAS nível 3 escritos para o processador, com vetorização e várias operações por ciclo. A queda da minha taxa acima de m = 2000, sem queda equivalente na LAPACK, é compatível com isso. "
        "Só medi separadamente o efeito da ordem dos laços; não isolei a contribuição dos blocos e a dos núcleos otimizados."
    ))
    adicionar(p("<b>2. A diferença de desempenho é a mesma para matrizes densas, em banda e tridiagonais?</b>", CORPO))
    adicionar(p(
        f"Não. A razão entre os tempos é de {faixa(razoes('LU densa'))} na LU densa e de {faixa(razoes('Cholesky densa'))} no Cholesky denso; cai para {faixa(razoes('LU em banda'))} na LU em banda e "
        f"{faixa(razoes('Cholesky em banda'))} no Cholesky em banda; e fica em {faixa(razoes('Thomas'), 2)} no Thomas. A vantagem da LAPACK vem de blocos e de operações sobre muitos dados, e esse espaço diminui com a banda. "
        f"Nas faixas estreitas, os laços são curtos e os dados cabem no cache, então a minha taxa nem sofre: {fixo(num(r4('bcsstk15', 'Cholesky em banda', 'implementado'), 'gflops'), 1)} GFLOP/s no Cholesky em banda da bcsstk15, contra "
        f"{fixo(num(r4('bcsstk15', 'Cholesky densa', 'implementado'), 'gflops'), 1)} no denso. As maiores razões em banda aparecem na 1138_bus, com b/m = 0,91, que na prática é densa. No caso tridiagonal, cada etapa faz poucas contas sobre poucos valores, sem "
        f"nada a agrupar em blocos: as duas rotinas ficam perto de 0,6 GFLOP/s, e o <font name='Mono'>dgtsv</font> leva uma vantagem pequena."
    ))
    adicionar(p("<b>3. Os erros são compatíveis com κ(A)u? Em que matrizes o erro é grande mesmo com resíduo pequeno?</b>", CORPO))
    t_final = e4("Thomas", "implementado", 1000000)
    adicionar(p(
        f"Sim. Nas {len(medicoes)} medições, a relação ‖x̂ − x‖/‖x‖ ≤ κ<sub>∞</sub>·(resíduo) vale sempre, e todas as matrizes reais ficam também abaixo de κ<sub>∞</sub>u. A única exceção a κ<sub>∞</sub>u é a SPD aleatória, com κ<sub>∞</sub> ≈ 2,3: em "
        f"m = 4000, o erro é {sci(num(e4('Cholesky densa', 'implementado', 4000), 'direto'))} contra κ<sub>∞</sub>u = {sci(num(e4('Cholesky densa', 'implementado', 4000), 'limite_kappa_u'))}. Isso não indica defeito: κu ignora as constantes que crescem com m, e com κ tão pequeno essas constantes aparecem. "
        f"O erro é grande, apesar do resíduo pequeno, nas matrizes mal condicionadas: bcsstk15 (LU densa: resíduo {sci(num(r4('bcsstk15', 'LU densa', 'implementado'), 'retroacao'))}, erro {sci(num(r4('bcsstk15', 'LU densa', 'implementado'), 'direto'))}), "
        f"1138_bus ({sci(num(r4('1138_bus', 'LU densa', 'implementado'), 'retroacao'))} e {sci(num(r4('1138_bus', 'LU densa', 'implementado'), 'direto'))}) e, principalmente, tridiag(−1, 2, −1) com m = 10<super>6</super> "
        f"({sci(num(t_final, 'retroacao'))} e {sci(num(t_final, 'direto'))}). Nesses casos, os métodos são estáveis: o resíduo pequeno mostra que eles resolvem exatamente um sistema muito próximo, e o erro maior em x vem do condicionamento, não do algoritmo."
    ))
    adicionar(p("<b>4. Quando vale a pena reordenar e usar a fatoração em banda em vez da densa?</b>", CORPO))
    orig_ch = t3("bcsstk15", "Cholesky", "original")
    rcm_ch = t3("bcsstk15", "Cholesky", "RCM")
    aumento_medido = num(rcm_ch, "fatoracao_banda_s") / num(orig_ch, "fatoracao_banda_s")
    aumento_previsto = num(rcm_ch, "flops_banda_exata") / num(orig_ch, "flops_banda_exata")
    adicionar(p(
        "A banda vale a pena quando b é bem menor que m; o reordenamento, quando reduz b. <b>Memória:</b> o armazenamento em banda só é menor que o denso se b<sub>I</sub> + b<sub>S</sub> + 1 &lt; m. Na 1138_bus e na steam2 em ordem original, ele chega a ser maior "
        f"({inteiro(bandas[('1138_bus', 'original')]['posicoes_lu'])} contra {inteiro(bandas[('1138_bus', 'original')]['posicoes_densa'])} posições na 1138_bus), e o RCM inverte isso ({inteiro(bandas[('1138_bus', 'RCM')]['posicoes_cholesky'])} posições no Cholesky em banda). "
        f"<b>Tempo:</b> com o RCM, a fatoração em banda ficou {fixo(ganho_medido('steam2', 'LU', 'RCM'), 0)} vezes mais rápida que a densa na steam2 e {fixo(ganho_medido('1138_bus', 'Cholesky', 'RCM'), 0)} vezes na 1138_bus (Cholesky). Sem reordenar, a 1138_bus (b/m = 0,91) leva o mesmo tempo da densa. "
        f"O RCM é barato em relação à fatoração densa: de {sci(num(t3('steam2', 'LU', 'RCM'), 'tempo_rcm_s'))} a {sci(num(t3('bcsstk15', 'LU', 'RCM'), 'tempo_rcm_s'))} s nestes testes, mesmo montando o grafo a partir da forma densa. Em relação à fatoração em banda, ele não é desprezível: "
        f"equivale a {fixo(100 * num(t3('steam2', 'LU', 'RCM'), 'tempo_rcm_s') / num(t3('steam2', 'LU', 'RCM'), 'fatoracao_banda_s'), 0)}% dela na steam2 e {fixo(100 * num(t3('1138_bus', 'Cholesky', 'RCM'), 'tempo_rcm_s') / num(t3('1138_bus', 'Cholesky', 'RCM'), 'fatoracao_banda_s'), 0)}% na 1138_bus. "
        f"<b>Mas ele nem sempre ajuda:</b> na bcsstk15, a banda aumenta de 437 para 501, e o Cholesky em banda fica {fixo(aumento_medido)} vezes mais lento ({sci(num(orig_ch, 'fatoracao_banda_s'), 2)} s sem RCM, {sci(num(rcm_ch, 'fatoracao_banda_s'), 2)} s com RCM), "
        f"o que concorda com a razão das contagens exatas, {fixo(aumento_previsto)}. Como essa matriz já tem banda estreita na numeração original, o certo é medir b antes e depois e ficar com a menor. <b>Estabilidade:</b> a LU em banda não pivota, então é preciso conferir o "
        "resíduo; nas matrizes SPD, o Cholesky em banda não tem esse risco. Por fim, o custo da solução também cai com a banda, o que aumenta a vantagem quando há muitos lados direitos."
    ))

    # 8. Conclusões
    adicionar(p("8. Conclusões", SECAO))
    adicionar(p(
        f"As implementações se mostraram corretas: o maior resíduo é {sci(max(num(l, 'retroacao') for l in medicoes if l['implementacao'] == 'implementado'))}, no mesmo patamar da LAPACK. A LAPACK é de {faixa(razoes('LU densa'), 0)} vezes mais rápida na forma densa, mas a vantagem encolhe com a banda e quase desaparece "
        "no caso tridiagonal, o que mostra que contar operações não basta: a organização dos acessos à memória é decisiva, como mostra também a ordem dos laços. Por outro lado, escolher a estrutura certa pesa mais que a qualidade da implementação: "
        f"meu Cholesky em banda fatora a bcsstk15 em {sci(num(orig_ch, 'fatoracao_banda_s'), 2)} s, bem menos que os {sci(num(r4('bcsstk15', 'Cholesky densa', 'LAPACK'), 'fatoracao_s'), 2)} s do <font name='Mono'>dpotrf</font> denso. O RCM costuma ser um passo barato e muito vantajoso, mas o resultado depende da matriz, "
        "como mostra a bcsstk15. O passo a passo para reproduzir cada tabela e figura está no README do código."
    ))

    pasta = ROOT / "relatorio"
    pasta.mkdir(exist_ok=True)
    saida = pasta / "relatorio.pdf"
    documento = SimpleDocTemplate(
        str(saida), pagesize=A4, leftMargin=1.9 * cm, rightMargin=1.9 * cm, topMargin=1.7 * cm, bottomMargin=1.7 * cm,
        title="Métodos diretos e comparação com a LAPACK", author="Augusto De Filippo Cavalini",
    )
    documento.build(conteudo, onFirstPage=rodape, onLaterPages=rodape)
    print(saida)


if __name__ == "__main__":
    main()
