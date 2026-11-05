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


void fatoracaoCholeskyLapacke(DenseMatrix& A) {
    int m = A.m;

    // A função LAPACKE_dpotrf sobrescreve o triângulo inferior de A com L.
    // Ela devolve 0 se tudo deu certo, e k > 0 se o pivô da etapa k não for positivo.
    lapack_int info = LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', m, A.a.data(), m);
    if (info > 0) {
        throw std::runtime_error("a matriz não é SPD: pivô não positivo na etapa " + std::to_string(info));
    }
}


void fatoracaoCholesky(DenseMatrix& A) {
    int m = A.m;

    int i, j, k;

    for (k = 0; k < m; k++) {
        // 1) Confere o pivô. Em uma matriz SPD, ele é sempre positivo.
        //    A condição está escrita assim para que um valor NaN também seja rejeitado.
        if (!(A(k, k) > 0.0)) {
            throw std::runtime_error("a matriz não é SPD: pivô não positivo na etapa " + std::to_string(k + 1));
        }

        // 2) A diagonal de L é a raiz quadrada do pivô.
        A(k, k) = std::sqrt(A(k, k));

        // 3) Calcula o restante da coluna k de L.
        for (i = k + 1; i < m; i++) {
            A(i, k) = A(i, k) / A(k, k);
        }

        // 4) Atualiza o triângulo inferior das colunas à direita: a_ij = a_ij - l_ik * l_jk.
        //    Como A é simétrica, cada coluna j só precisa das linhas i >= j.
        //    O laço interno anda ao longo da coluna j, a ordem em que os
        //    coeficientes estão guardados na memória.
        for (j = k + 1; j < m; j++) {
            double l_jk = A(j, k);

            for (i = j; i < m; i++) {
                A(i, j) -= A(i, k) * l_jk;
            }
        }
    }
}


void substituicao_direta_cholesky(const DenseMatrix& L, std::vector<double>& y) {
    int m = L.m;

    // Resolve Ly = y de cima para baixo. Aqui a diagonal de L não vale 1, então é preciso dividir.
    //
    // Quando as parcelas à esquerda já foram descontadas, a divisão fornece y[j].
    // Depois, descontamos l_ij * y[j] das equações abaixo.
    for (int j = 0; j < m; j++) {
        y[j] = y[j] / L(j, j);

        for (int i = j + 1; i < m; i++) {
            y[i] -= L(i, j) * y[j];
        }
    }
}


void substituicao_inversa_cholesky(const DenseMatrix& L, std::vector<double>& x) {
    int m = L.m;

    // Resolve L^T x = x de baixo para cima, usando apenas o triângulo inferior.
    //
    // A linha j de L^T é a coluna j de L. Por isso, somar as parcelas à direita da
    // diagonal significa percorrer a coluna j de L, na ordem da memória.
    for (int j = m - 1; j >= 0; j--) {
        double soma = 0.0;

        for (int i = j + 1; i < m; i++) {
            soma += L(i, j) * x[i];
        }

        x[j] = (x[j] - soma) / L(j, j);
    }
}


std::vector<double> resolver_sistema_cholesky(const DenseMatrix& L, const std::vector<double>& b) {
    // Fazemos uma cópia do vetor b, porque cada etapa sobrescreve o vetor recebido.
    std::vector<double> x = b;

    substituicao_direta_cholesky(L, x);      // x passa a ser y, de Ly = b
    substituicao_inversa_cholesky(L, x);     // x passa a ser a solução, de L^T x = y

    return x;
}
