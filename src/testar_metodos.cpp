// Testa todos os métodos da Parte 2 em uma matriz, comparando cada um com a LAPACK.
//
// Uso: testar_metodos [arquivo.mtx]
// Sem argumentos, usa matrices/teste_10x10.mtx.
//
// A solução exata é x = (1, ..., 1), e b é calculado como b = A*x.
// Um método que não se aplica à matriz aparece na tabela com o motivo.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

#include "parte2.hpp"

typedef std::vector<double> Vetor;

// O teste passa se o erro de retroação for no máximo 100 * m * u, com u = 1,1e-16.
const double FATOR_DE_TOLERANCIA = 100.0;
const double U = 1.1e-16;

// Se a matriz for pequena, as soluções são impressas por extenso.
const int MAIOR_M_PARA_IMPRIMIR_X = 20;

// Calcula o produto y = A*x, percorrendo a matriz coluna por coluna.
static Vetor multiplicar(const DenseMatrix& A, const Vetor& x) {
    int m = A.m;
    Vetor y(m, 0.0);
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < m; i++) {
            y[i] += A(i, j) * x[j];
        }
    }
    return y;
}

// ||v||_inf: o maior valor absoluto do vetor.
static double norma_inf_vetor(const Vetor& v) {
    double maior = 0.0;
    for (size_t i = 0; i < v.size(); i++) {
        maior = std::max(maior, std::abs(v[i]));
    }
    return maior;
}

// ||A||_inf: a maior soma dos valores absolutos de uma linha.
static double norma_inf_matriz(const DenseMatrix& A) {
    int m = A.m;
    Vetor soma_da_linha(m, 0.0);
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < m; i++) {
            soma_da_linha[i] += std::abs(A(i, j));
        }
    }
    return norma_inf_vetor(soma_da_linha);
}

static bool eh_simetrica(const DenseMatrix& A) {
    int m = A.m;
    for (int j = 0; j < m; j++) {
        for (int i = j + 1; i < m; i++) {
            if (A(i, j) != A(j, i)) {
                return false;
            }
        }
    }
    return true;
}

// Cada teste abaixo calcula a solução com o nosso método e com a rotina da LAPACK.
// Se o método não se aplicar à matriz, ele lança uma exceção com o motivo.

static void testar_lu_densa(const DenseMatrix& A, const Vetor& b, Vetor& x_nosso, Vetor& x_lapack) {
    int m = A.m;
    std::vector<lapack_int> pivos(m);

    DenseMatrix LU = A;
    fatoracaoLu(LU, pivos);
    x_nosso = resolver_sistema_lu(LU, pivos, b);

    DenseMatrix LU_lapack = A;
    fatoracaoLULapacke(LU_lapack, pivos);
    x_lapack = b;
    LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', m, 1, LU_lapack.a.data(), m, pivos.data(), x_lapack.data(), m);
}

static void testar_cholesky_densa(const DenseMatrix& A, const Vetor& b, Vetor& x_nosso, Vetor& x_lapack) {
    int m = A.m;

    // O Cholesky só lê um triângulo, então não perceberia sozinho uma matriz não simétrica.
    if (!eh_simetrica(A)) {
        throw std::runtime_error("a matriz não é simétrica");
    }

    DenseMatrix L = A;
    fatoracaoCholesky(L);
    x_nosso = resolver_sistema_cholesky(L, b);

    DenseMatrix L_lapack = A;
    fatoracaoCholeskyLapacke(L_lapack);
    x_lapack = b;
    LAPACKE_dpotrs(LAPACK_COL_MAJOR, 'L', m, 1, L_lapack.a.data(), m, x_lapack.data(), m);
}

static void testar_lu_banda(const DenseMatrix& A, const Vetor& b, Vetor& x_nosso, Vetor& x_lapack) {
    int m = A.m;
    int bI = 0;
    int bS = 0;
    calcular_bandas(A, bI, bS);

    BandMatrix LU = montar_banda(A, bI, bS);
    fatoracaoLuBanda(LU);
    x_nosso = resolver_sistema_lu_banda(LU, b);

    // O dgbtrf pivota e por isso precisa de bI posições extras acima da banda.
    BandMatrix LU_lapack = montar_banda(A, bI, bI + bS);
    std::vector<lapack_int> pivos(m);
    fatoracaoLuBandaLapacke(LU_lapack, pivos);
    x_lapack = b;
    LAPACKE_dgbtrs(LAPACK_COL_MAJOR, 'N', m, bI, bS, 1, LU_lapack.a.data(), LU_lapack.linhas, pivos.data(), x_lapack.data(), m);
}

static void testar_cholesky_banda(const DenseMatrix& A, const Vetor& b, Vetor& x_nosso, Vetor& x_lapack) {
    int m = A.m;
    int bI = 0;
    int bS = 0;
    calcular_bandas(A, bI, bS);

    if (!eh_simetrica(A)) {
        throw std::runtime_error("a matriz não é simétrica");
    }

    BandMatrix L = montar_banda(A, bI, 0);
    fatoracaoCholeskyBanda(L);
    x_nosso = resolver_sistema_cholesky_banda(L, b);

    BandMatrix L_lapack = montar_banda(A, bI, 0);
    fatoracaoCholeskyBandaLapacke(L_lapack);
    x_lapack = b;
    LAPACKE_dpbtrs(LAPACK_COL_MAJOR, 'L', m, bI, 1, L_lapack.a.data(), L_lapack.linhas, x_lapack.data(), m);
}

static void testar_thomas(const DenseMatrix& A, const Vetor& b, Vetor& x_nosso, Vetor& x_lapack) {
    int m = A.m;
    int bI = 0;
    int bS = 0;
    calcular_bandas(A, bI, bS);

    if (bI > 1 || bS > 1) {
        throw std::runtime_error("a matriz não é tridiagonal");
    }

    // Separa a matriz nos três vetores.
    Vetor diagonal(m);
    Vetor inferior(m - 1);
    Vetor superior(m - 1);
    for (int k = 0; k < m; k++) {
        diagonal[k] = A(k, k);
    }
    for (int k = 0; k < m - 1; k++) {
        inferior[k] = A(k + 1, k);
        superior[k] = A(k, k + 1);
    }

    // O dgtsv sobrescreve os vetores, então ele recebe cópias.
    Vetor inferior_lapack = inferior;
    Vetor diagonal_lapack = diagonal;
    Vetor superior_lapack = superior;
    x_lapack = resolver_tridiagonal_lapacke(inferior_lapack, diagonal_lapack, superior_lapack, b);

    fatoracaoThomas(inferior, diagonal, superior);
    x_nosso = resolver_sistema_thomas(inferior, diagonal, superior, b);
}

// O tipo de uma função de teste: todas as funções acima têm esta mesma forma.
typedef void (*FuncaoTeste)(const DenseMatrix&, const Vetor&, Vetor&, Vetor&);

static void imprimir_vetor(const char* nome, const Vetor& v) {
    std::printf("  %-8s", nome);
    for (size_t i = 0; i < v.size(); i++) {
        std::printf(" %.15g", v[i]);
    }
    std::printf("\n");
}

// Roda um teste e imprime o resultado. Devolve false se um método aplicável falhar.
static bool rodar(const char* nome, FuncaoTeste testar, const DenseMatrix& A, const Vetor& b) {
    int m = A.m;
    Vetor x_nosso;
    Vetor x_lapack;

    try {
        testar(A, b, x_nosso, x_lapack);
    } catch (const std::exception& erro) {
        std::printf("%-20s não se aplica: %s\n\n", nome, erro.what());
        return true;
    }

    // Erro de retroação: ||b - A*x||_inf / (||A||_inf * ||x||_inf).
    Vetor Ax = multiplicar(A, x_nosso);
    Vetor residuo(m);
    for (int i = 0; i < m; i++) {
        residuo[i] = b[i] - Ax[i];
    }
    double erro_retroacao = norma_inf_vetor(residuo) / (norma_inf_matriz(A) * norma_inf_vetor(x_nosso));

    // Erro direto: a solução exata é x = (1, ..., 1), que tem ||x||_inf = 1.
    Vetor erro_de_x(m);
    Vetor diferenca_para_lapack(m);
    for (int i = 0; i < m; i++) {
        erro_de_x[i] = x_nosso[i] - 1.0;
        diferenca_para_lapack[i] = x_nosso[i] - x_lapack[i];
    }
    double erro_direto = norma_inf_vetor(erro_de_x);
    double diferenca = norma_inf_vetor(diferenca_para_lapack);

    double tolerancia = FATOR_DE_TOLERANCIA * m * U;
    bool passou = (erro_retroacao <= tolerancia);

    const char* situacao = "VERIFICAR";
    if (passou) {
        situacao = "OK";
    }

    std::printf("%-20s %-9s retroação = %.1e   erro direto = %.1e   diferença para a LAPACK = %.1e\n",
                nome, situacao, erro_retroacao, erro_direto, diferenca);

    if (m <= MAIOR_M_PARA_IMPRIMIR_X) {
        imprimir_vetor("nosso:", x_nosso);
        imprimir_vetor("LAPACK:", x_lapack);
    }
    std::printf("\n");

    return passou;
}

int main(int argc, char** argv) {
    std::string caminho = "matrices/teste_10x10.mtx";
    if (argc > 1) {
        caminho = argv[1];
    }

    try {
        MtxFile arquivo = read_mtx(caminho);
        const DenseMatrix& A = arquivo.A;

        Vetor x_exato(A.m, 1.0);
        Vetor b = multiplicar(A, x_exato);

        std::printf("Matriz %s: m = %d, solução exata x = (1, ..., 1)\n\n", arquivo.name.c_str(), A.m);

        bool tudo_ok = true;
        tudo_ok = rodar("LU densa", testar_lu_densa, A, b) && tudo_ok;
        tudo_ok = rodar("Cholesky densa", testar_cholesky_densa, A, b) && tudo_ok;
        tudo_ok = rodar("LU em banda", testar_lu_banda, A, b) && tudo_ok;
        tudo_ok = rodar("Cholesky em banda", testar_cholesky_banda, A, b) && tudo_ok;
        tudo_ok = rodar("Thomas", testar_thomas, A, b) && tudo_ok;

        if (!tudo_ok) {
            std::printf("Algum método ficou com erro de retroação acima da tolerância.\n");
            return 1;
        }
        std::printf("Todos os métodos aplicáveis passaram.\n");
    } catch (const std::exception& erro) {
        std::fprintf(stderr, "erro: %s\n", erro.what());
        return 1;
    }

    return 0;
}
