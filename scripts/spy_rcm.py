"""Padrão de esparsidade antes e depois do RCM, lado a lado, gravado em figures/rcm_<matriz>.png.

Usa a permutação gravada pelo programa parte3 em results/rcm_<matriz>.txt.
Uso: conda run -n tp-alc python scripts/spy_rcm.py arquivo.mtx [arquivo.mtx ...]
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # gera arquivos de imagem, sem abrir janela
import matplotlib.pyplot as plt
import numpy as np
import scipy.io

from spy import ROOT, tamanho_do_marcador


def larguras_de_banda(A):
    """Devolve (bI, bS): as maiores distâncias dos não nulos até a diagonal."""
    coo = A.tocoo()
    if coo.nnz == 0:
        return 0, 0
    distancia = coo.row - coo.col
    bI = max(int(distancia.max()), 0)
    bS = max(int(-distancia.min()), 0)
    return bI, bS


def desenhar(eixos, A, titulo):
    m = A.shape[0]
    bI, bS = larguras_de_banda(A)
    eixos.spy(A, markersize=tamanho_do_marcador(m), color="black")
    eixos.set_title(f"{titulo}\n$b_I$ = {bI}, $b_S$ = {bS}", fontsize=10)
    eixos.tick_params(labelsize=8)


def figura_rcm(caminho, pasta_de_saida):
    """Gera a figura de uma matriz e devolve o caminho da imagem."""
    arquivo_perm = ROOT / "results" / f"rcm_{caminho.stem}.txt"
    if not arquivo_perm.exists():
        sys.exit(f"{arquivo_perm} não existe: rode antes ./build/parte3 {caminho}")

    # O mmread já devolve as matrizes simétricas completas.
    A = scipy.io.mmread(caminho, spmatrix=False).tocsr()
    A.eliminate_zeros()  # zeros guardados explicitamente não contam como não nulos

    perm = np.loadtxt(arquivo_perm, dtype=int, ndmin=1)
    if len(perm) != A.shape[0]:
        sys.exit(f"{arquivo_perm} não corresponde a {caminho}")

    A_permutada = A[perm][:, perm]  # PAP^T

    figura, (esquerda, direita) = plt.subplots(1, 2, figsize=(8, 4.4))
    desenhar(esquerda, A, f"{caminho.stem}: ordem original")
    desenhar(direita, A_permutada, f"{caminho.stem}: após RCM")
    figura.tight_layout()

    saida = pasta_de_saida / f"rcm_{caminho.stem}.png"
    figura.savefig(saida, dpi=200)
    plt.close(figura)
    return saida


def main():
    arquivos = [Path(a) for a in sys.argv[1:]]
    if not arquivos:
        sys.exit("uso: spy_rcm.py arquivo.mtx [arquivo.mtx ...]")

    pasta_de_saida = ROOT / "figures"
    pasta_de_saida.mkdir(exist_ok=True)

    for arquivo in arquivos:
        print(figura_rcm(arquivo, pasta_de_saida))


if __name__ == "__main__":
    main()
