#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include <cmath>
#include <stdexcept>

#include "metodos.hpp"


DenseMatrix fatoracaoLULapacke(const DenseMatrix& A, std::vector<lapack_int>& pivos) {
    int m = A.m;

    // Fazemos uma cópia da matriz, porque a fatoração sobrescreve os coeficientes.
    DenseMatrix LU = A;

    // A função LAPACKE_dgetrf devolve 0 se tudo deu certo, e > 0 se encontrou um pivô nulo.
    lapack_int info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, m, m, LU.a.data(), m, pivos.data());
    if (info > 0) {
        throw std::runtime_error("a matriz é singular: U tem um zero na diagonal");
    }

    return LU;
}


int linha_do_pivo(const DenseMatrix& LU, int k) {
    int m = LU.m;

    double max_abs = 0.0;
    int linha_do_pivo = k;
    for (int i = k; i < m; i++) {
        double valor_absoluto = std::abs(LU(i, k));
        if (valor_absoluto > max_abs) {
            max_abs = valor_absoluto;
            linha_do_pivo = i;
        }
    }
    return linha_do_pivo;
}

DenseMatrix fatoracaoLu(const DenseMatrix& A, std::vector<lapack_int>& pivos) {
    int m = A.m;

    // Fazemos uma cópia da matriz, porque a fatoração sobrescreve os coeficientes.
    // L e U ficam nesta única matriz: U no triângulo superior e L abaixo da diagonal.
    DenseMatrix LU = A;

    int i, j, k;

    for (k = 0; k < m; k++) {
        // 1) Escolhe o pivô: o maior valor absoluto na coluna k,
        //    entre as linhas k, k+1, ..., m-1.
        int lin_piv = linha_do_pivo(LU, k);

        // A LAPACK conta as linhas a partir de 1, e o dgetrs espera os pivôs nessa convenção.
        pivos[k] = lin_piv + 1;

        // 2) Troca a linha k com a linha do pivô, inteira. Isso inclui as colunas
        //    anteriores a k, onde já estão os multiplicadores de L.
        if (lin_piv != k) {
            for (j = 0; j < m; j++) {
                double aux = LU(k, j);
                LU(k, j) = LU(lin_piv, j);
                LU(lin_piv, j) = aux;
            }
        }

        if (LU(k, k) == 0.0) {
            throw std::runtime_error("a matriz é singular: pivô nulo");
        }

        for (i = k + 1; i < m; i++) {
            // 3) Calcula o multiplicador.
            double multiplicador = LU(i, k) / LU(k, k);
            LU(i, k) = multiplicador;

            // 4) Atualiza a linha i.
            for (j = k + 1; j < m; j++) {
                LU(i, j) -= multiplicador * LU(k, j);
            }
        }
    }

    return LU;
}