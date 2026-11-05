#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include <cmath>
#include <stdexcept>

#include "metodos.hpp"


void fatoracaoLULapacke(DenseMatrix& A, std::vector<lapack_int>& pivos) {
    int m = A.m;

    // A função LAPACKE_dgetrf sobrescreve A com os fatores.
    // Ela devolve 0 se tudo deu certo, e > 0 se encontrou um pivô nulo.
    lapack_int info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, m, m, A.a.data(), m, pivos.data());
    if (info > 0) {
        throw std::runtime_error("a matriz é singular: U tem um zero na diagonal");
    }
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

void fatoracaoLu(DenseMatrix& A, std::vector<lapack_int>& pivos) {
    int m = A.m;

    int i, j, k;

    for (k = 0; k < m; k++) {
        // 1) Escolhe o pivô: o maior valor absoluto na coluna k,
        //    entre as linhas k, k+1, ..., m-1.
        int lin_piv = linha_do_pivo(A, k);

        // A LAPACK conta as linhas a partir de 1. Usamos a mesma convenção.
        pivos[k] = lin_piv + 1;

        // 2) Troca a linha k com a linha do pivô, inteira. Isso inclui as colunas
        //    anteriores a k, onde já estão os multiplicadores de L.
        if (lin_piv != k) {
            for (j = 0; j < m; j++) {
                double aux = A(k, j);
                A(k, j) = A(lin_piv, j);
                A(lin_piv, j) = aux;
            }
        }

        if (A(k, k) == 0.0) {
            throw std::runtime_error("a matriz é singular: pivô nulo");
        }

        // 3) Calcula os multiplicadores da coluna k e os guarda em A(i, k),
        //    as posições que seriam zeradas.
        for (i = k + 1; i < m; i++) {
            A(i, k) = A(i, k) / A(k, k);
        }

        // 4) Atualiza as colunas à direita da coluna k, uma de cada vez.
        //    O laço interno anda ao longo da coluna j, a ordem em que os
        //    coeficientes estão guardados na memória.
        for (j = k + 1; j < m; j++) {
            double valor_da_linha_k = A(k, j);

            for (i = k + 1; i < m; i++) {
                A(i, j) -= A(i, k) * valor_da_linha_k;
            }
        }
    }
}


void aplicar_pivos(const std::vector<lapack_int>& pivos, std::vector<double>& b) {
    int m = static_cast<int>(b.size());

    // Repete em b as trocas de linha da fatoração, na mesma ordem. O resultado é Pb.
    for (int k = 0; k < m; k++) {
        int lin_piv = pivos[k] - 1;  // pivos conta a partir de 1

        double aux = b[k];
        b[k] = b[lin_piv];
        b[lin_piv] = aux;
    }
}


void substituicao_direta(const DenseMatrix& LU, std::vector<double>& y) {
    int m = LU.m;

    // Resolve Ly = y de cima para baixo. A diagonal de L vale 1, então não há divisão.
    //
    // Quando y[j] fica pronto, descontamos sua parcela l_ij * y[j] das equações abaixo.
    // Isso equivale a y_i = b_i - soma(l_ij * y_j), mas percorre a matriz coluna por
    // coluna, na mesma ordem em que os coeficientes estão guardados na memória.
    for (int j = 0; j < m; j++) {
        for (int i = j + 1; i < m; i++) {
            y[i] -= LU(i, j) * y[j];
        }
    }
}


void substituicao_inversa(const DenseMatrix& LU, std::vector<double>& x) {
    int m = LU.m;

    // Resolve Ux = x de baixo para cima, também coluna por coluna.
    //
    // Quando as parcelas à direita já foram descontadas, basta dividir pela diagonal
    // para obter x[j]. Depois, descontamos u_ij * x[j] das equações acima.
    for (int j = m - 1; j >= 0; j--) {
        x[j] = x[j] / LU(j, j);

        for (int i = 0; i < j; i++) {
            x[i] -= LU(i, j) * x[j];
        }
    }
}


std::vector<double> resolver_sistema_lu(const DenseMatrix& LU, const std::vector<lapack_int>& pivos, const std::vector<double>& b) {
    // Fazemos uma cópia do vetor b, porque cada etapa sobrescreve o vetor recebido.
    std::vector<double> x = b;

    aplicar_pivos(pivos, x);         // x passa a ser Pb
    substituicao_direta(LU, x);      // x passa a ser y, de Ly = Pb
    substituicao_inversa(LU, x);     // x passa a ser a solução, de Ux = y

    return x;
}
