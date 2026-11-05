"""Gráficos log-log do estudo de escala da Parte 4 e inclinações ajustadas.

Lê results/parte4_escala.csv, gravado por "parte4 escala". Grava figures/parte4_escala.png
e results/parte4_inclinacoes.csv, e imprime a tabela de inclinações em Markdown.

Uso: conda run -n tp-alc python scripts/graficos_escala.py
"""
import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # gera arquivos de imagem, sem abrir janela
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter

ROOT = Path(__file__).resolve().parent.parent

# A cor identifica a implementação, em todos os painéis.
COR_NOSSA = "#2a78d6"
COR_LAPACK = "#eb6834"
COR_REFERENCIA = "#898781"
COR_GRADE = "#e1e0d9"
COR_TEXTO = "#52514e"

# O traço e o marcador identificam o método, para não depender só da cor.
ESTILOS = {
    "LU densa": {"linestyle": "-", "marker": "o"},
    "Cholesky densa": {"linestyle": "--", "marker": "s"},
    "Thomas": {"linestyle": "-", "marker": "o"},
}

NOME_DA_ROTINA = {
    ("LU densa", "implementado"): "LU implementado",
    ("LU densa", "LAPACK"): "dgetrf/dgetrs",
    ("Cholesky densa", "implementado"): "Cholesky implementado",
    ("Cholesky densa", "LAPACK"): "dpotrf/dpotrs",
    ("Thomas", "implementado"): "Thomas implementado",
    ("Thomas", "LAPACK"): "dgtsv",
}


def ler_resultados():
    caminho = ROOT / "results" / "parte4_escala.csv"
    if not caminho.exists():
        sys.exit(f"{caminho} não existe: rode antes make parte4_escala")

    with open(caminho, newline="") as arquivo:
        linhas = list(csv.DictReader(arquivo))

    series = {}
    for linha in linhas:
        chave = (linha["metodo"], linha["implementacao"])
        serie = series.setdefault(chave, {"m": [], "fatoracao": [], "solucao": []})
        serie["m"].append(float(linha["m"]))
        serie["fatoracao"].append(float(linha["fatoracao_s"]))
        serie["solucao"].append(float(linha["solucao_s"]))

    for serie in series.values():
        ordem = np.argsort(serie["m"])
        for campo in serie:
            serie[campo] = np.array(serie[campo])[ordem]
    return series


def inclinacao(m, tempo):
    """Inclinação da reta ajustada a log(tempo) em função de log(m)."""
    return float(np.polyfit(np.log10(m), np.log10(tempo), 1)[0])


def inclinacao_final(m, tempo):
    """Inclinação entre os dois últimos pontos, que mostra o efeito das matrizes maiores."""
    return float(np.log10(tempo[-1] / tempo[-2]) / np.log10(m[-1] / m[-2]))


def desenhar_referencia(eixos, m, tempo, expoente):
    """Reta de inclinação prevista, passando pela mediana geométrica da série."""
    x = np.array([m[0], m[-1]])
    centro_x = np.sqrt(m[0] * m[-1])
    centro_y = np.exp(np.mean(np.log(tempo)))
    y = centro_y * (x / centro_x) ** expoente
    eixos.plot(x, y, color=COR_REFERENCIA, linestyle=":", linewidth=1.5,
               label=f"inclinação {expoente} (prevista)", zorder=1)


def desenhar_painel(eixos, titulo, curvas, expoente):
    """Desenha um painel: curvas é uma lista de (método, implementação, m, tempo)."""
    referencia = curvas[0]
    desenhar_referencia(eixos, referencia[2], referencia[3], expoente)

    for metodo, implementacao, m, tempo in curvas:
        cor = COR_NOSSA if implementacao == "implementado" else COR_LAPACK
        estilo = ESTILOS[metodo]
        eixos.plot(m, tempo, color=cor, linewidth=2, markersize=6,
                   markeredgecolor="white", markeredgewidth=1,
                   label=NOME_DA_ROTINA[(metodo, implementacao)], zorder=3, **estilo)

    eixos.set_xscale("log")
    eixos.set_yscale("log")

    # O eixo m mostra só os valores medidos, sem os rótulos automáticos, que se sobrepõem.
    valores_de_m = referencia[2]
    eixos.set_xticks(valores_de_m)
    eixos.set_xticklabels([f"{int(valor)}" for valor in valores_de_m])
    eixos.xaxis.set_minor_formatter(NullFormatter())

    eixos.set_title(titulo, fontsize=10, color="#0b0b0b")
    eixos.set_xlabel("m", fontsize=9, color=COR_TEXTO)
    eixos.set_ylabel("tempo (s)", fontsize=9, color=COR_TEXTO)
    eixos.tick_params(labelsize=8, colors=COR_TEXTO)
    eixos.grid(True, which="major", color=COR_GRADE, linewidth=0.8)
    eixos.set_axisbelow(True)
    for lado in ("top", "right"):
        eixos.spines[lado].set_visible(False)
    for lado in ("left", "bottom"):
        eixos.spines[lado].set_color(COR_GRADE)
    eixos.legend(fontsize=8, frameon=False)


def main():
    series = ler_resultados()

    figura, (densa, solucao, tridiagonal) = plt.subplots(1, 3, figsize=(12, 4.2))

    metodos_densos = [("LU densa", "implementado"), ("LU densa", "LAPACK"), ("Cholesky densa", "implementado"), ("Cholesky densa", "LAPACK")]

    curvas = [(met, imp, series[(met, imp)]["m"], series[(met, imp)]["fatoracao"]) for met, imp in metodos_densos]
    desenhar_painel(densa, "Fatoração densa", curvas, 3)

    curvas = [(met, imp, series[(met, imp)]["m"], series[(met, imp)]["solucao"]) for met, imp in metodos_densos]
    desenhar_painel(solucao, "Solução densa", curvas, 2)

    # O dgtsv só informa o tempo total, então o painel usa fatoração mais solução nos dois.
    curvas = []
    for implementacao in ("implementado", "LAPACK"):
        serie = series[("Thomas", implementacao)]
        curvas.append(("Thomas", implementacao, serie["m"], serie["fatoracao"] + serie["solucao"]))
    desenhar_painel(tridiagonal, "tridiag(-1, 2, -1): fatoração + solução", curvas, 1)

    figura.tight_layout()

    pasta_figuras = ROOT / "figures"
    pasta_figuras.mkdir(exist_ok=True)
    saida = pasta_figuras / "parte4_escala.png"
    figura.savefig(saida, dpi=200)
    plt.close(figura)
    print(saida)

    # Tabela de inclinações.
    linhas = []
    for metodo, implementacao in metodos_densos:
        serie = series[(metodo, implementacao)]
        nome = NOME_DA_ROTINA[(metodo, implementacao)]
        linhas.append((nome, "fatoração", 3, inclinacao(serie["m"], serie["fatoracao"]), inclinacao_final(serie["m"], serie["fatoracao"])))
        linhas.append((nome, "solução", 2, inclinacao(serie["m"], serie["solucao"]), inclinacao_final(serie["m"], serie["solucao"])))
    for implementacao in ("implementado", "LAPACK"):
        serie = series[("Thomas", implementacao)]
        total = serie["fatoracao"] + serie["solucao"]
        linhas.append((NOME_DA_ROTINA[("Thomas", implementacao)], "fatoração + solução", 1, inclinacao(serie["m"], total), inclinacao_final(serie["m"], total)))

    print("\n| Rotina | Etapa | Prevista | Ajustada | Dois últimos pontos |")
    print("|---|---|---:|---:|---:|")
    for nome, etapa, prevista, ajustada, final in linhas:
        print(f"| {nome} | {etapa} | {prevista} | {ajustada:.2f} | {final:.2f} |")

    with open(ROOT / "results" / "parte4_inclinacoes.csv", "w", newline="") as arquivo:
        escritor = csv.writer(arquivo)
        escritor.writerow(["rotina", "etapa", "prevista", "ajustada", "dois_ultimos_pontos"])
        for nome, etapa, prevista, ajustada, final in linhas:
            escritor.writerow([nome, etapa, prevista, f"{ajustada:.4f}", f"{final:.4f}"])


if __name__ == "__main__":
    main()
