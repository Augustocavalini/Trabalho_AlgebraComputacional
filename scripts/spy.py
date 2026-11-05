"""Spy plots (padrão de esparsidade) das matrizes em matrices/*.mtx, gravados em figures/.

Uso: conda run -n tp-alc python scripts/spy.py [arquivo.mtx ...]
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # gera arquivos de imagem, sem abrir janela
import matplotlib.pyplot as plt
import scipy.io

ROOT = Path(__file__).resolve().parent.parent


def tamanho_do_marcador(m):
    """Marcador menor para matrizes maiores, para o padrão não virar um borrão."""
    tamanho = 600 / m
    if tamanho > 2.0:
        tamanho = 2.0
    if tamanho < 0.05:
        tamanho = 0.05
    return tamanho


def spy(caminho, pasta_de_saida):
    """Gera a figura de uma matriz e devolve o caminho da imagem."""
    # O mmread já devolve as matrizes simétricas completas.
    A = scipy.io.mmread(caminho, spmatrix=False)
    A = A.tocsr()
    A.eliminate_zeros()  # zeros guardados explicitamente não contam como não nulos
    m = A.shape[0]

    figura, eixos = plt.subplots(figsize=(4.5, 4.5))
    eixos.spy(A, markersize=tamanho_do_marcador(m), color="black")
    eixos.set_title(f"{caminho.stem}  (m = {m}, nnz = {A.nnz})", fontsize=10)
    eixos.tick_params(labelsize=8)
    figura.tight_layout()

    saida = pasta_de_saida / f"spy_{caminho.stem}.png"
    figura.savefig(saida, dpi=200)
    plt.close(figura)
    return saida


def main():
    # Usa os arquivos da linha de comando ou, se não houver, os de matrices/.
    argumentos = sys.argv[1:]
    if argumentos:
        arquivos = [Path(a) for a in argumentos]
    else:
        arquivos = sorted((ROOT / "matrices").glob("*.mtx"))

    if not arquivos:
        sys.exit("nenhum .mtx encontrado em matrices/")

    pasta_de_saida = ROOT / "figures"
    pasta_de_saida.mkdir(exist_ok=True)

    for arquivo in arquivos:
        print(spy(arquivo, pasta_de_saida))


if __name__ == "__main__":
    main()
