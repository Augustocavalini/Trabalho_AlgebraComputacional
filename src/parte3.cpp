// Parte 3, itens 1 e 2: reordenamento RCM e efeito sobre a largura de banda.
//
// Uso: parte3 arquivo.mtx [arquivo.mtx ...]
//
// Para cada matriz:
// - calcula a permutação RCM e forma PAP^T;
// - confere que o sistema reordenado (PAP^T)(Px) = Pb é equivalente ao original;
// - imprime as bandas e o armazenamento em banda antes e depois;
// - grava a permutação em results/rcm_<matriz>.txt, usada por scripts/spy_rcm.py.
//
// A tabela também é gravada em results/parte3.csv.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

#include "matrix.hpp"
#include "rcm.hpp"

// Uma ordem (antes ou depois do RCM) de uma matriz.
struct LinhaTabela {
    std::string nome;
    std::string ordem;
    int m = 0;
    int bI = 0;
    int bS = 0;
    bool simetrica = false;
    long posicoes_lu = 0;        // (bI + bS + 1) * m
    long posicoes_cholesky = 0;  // (b + 1) * m, apenas se a matriz for simétrica
    long posicoes_densa = 0;     // m * m
};

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

static long contar_nao_nulos(const DenseMatrix& A) {
    long total = 0;
    for (size_t p = 0; p < A.a.size(); p++) {
        if (A.a[p] != 0.0) {
            total = total + 1;
        }
    }
    return total;
}

static LinhaTabela medir(const std::string& nome, const std::string& ordem, const DenseMatrix& A) {
    LinhaTabela l;
    l.nome = nome;
    l.ordem = ordem;
    l.m = A.m;
    calcular_bandas(A, l.bI, l.bS);
    l.simetrica = eh_simetrica(A);

    long m = A.m;
    long linhas_lu = l.bI + l.bS + 1;
    l.posicoes_lu = linhas_lu * m;
    l.posicoes_densa = m * m;
    if (l.simetrica) {
        long linhas_cholesky = l.bI + 1;
        l.posicoes_cholesky = linhas_cholesky * m;
    }
    return l;
}

// Confere que perm contém cada índice de 0 a m - 1 exatamente uma vez.
static bool eh_permutacao(const std::vector<int>& perm) {
    int m = static_cast<int>(perm.size());
    std::vector<bool> visto(m, false);

    for (int k = 0; k < m; k++) {
        int indice = perm[k];
        if (indice < 0 || indice >= m || visto[indice]) {
            return false;
        }
        visto[indice] = true;
    }
    return true;
}

// Calcula o produto y = A*x, percorrendo a matriz coluna por coluna.
static std::vector<double> multiplicar(const DenseMatrix& A, const std::vector<double>& x) {
    int m = A.m;
    std::vector<double> y(m, 0.0);
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < m; i++) {
            y[i] += A(i, j) * x[j];
        }
    }
    return y;
}

static double norma_inf_vetor(const std::vector<double>& v) {
    double maior = 0.0;
    for (size_t i = 0; i < v.size(); i++) {
        maior = std::max(maior, std::abs(v[i]));
    }
    return maior;
}

// Mede o quanto (PAP^T)(Px) difere de Pb, relativamente a ||b|| + 1, para x = (1, ..., m).
// Em aritmética exata o resultado é 0. Um x com valores diferentes revela permutações erradas,
// que passariam despercebidas com x = (1, ..., 1).
static double residuo_do_sistema_reordenado(const DenseMatrix& A, const DenseMatrix& permutada, const std::vector<int>& perm) {
    int m = A.m;

    std::vector<double> x(m);
    for (int i = 0; i < m; i++) {
        x[i] = i + 1;
    }

    std::vector<double> b = multiplicar(A, x);
    std::vector<double> Pb = permutar_vetor(b, perm);
    std::vector<double> Px = permutar_vetor(x, perm);
    std::vector<double> lado_esquerdo = multiplicar(permutada, Px);

    std::vector<double> diferenca(m);
    for (int i = 0; i < m; i++) {
        diferenca[i] = lado_esquerdo[i] - Pb[i];
    }
    return norma_inf_vetor(diferenca) / (norma_inf_vetor(b) + 1.0);
}

static bool gravar_permutacao(const std::string& nome, const std::vector<int>& perm) {
    std::string caminho = "results/rcm_" + nome + ".txt";

    FILE* arquivo = std::fopen(caminho.c_str(), "w");
    if (arquivo == nullptr) {
        return false;
    }
    for (size_t k = 0; k < perm.size(); k++) {
        std::fprintf(arquivo, "%d\n", perm[k]);
    }
    std::fclose(arquivo);
    return true;
}

static void imprimir_tabela(const std::vector<LinhaTabela>& tabela) {
    std::printf("\n| Matriz | Ordem | m | b_I | b_S | LU em banda: (b_I+b_S+1)m | Cholesky em banda: (b+1)m | Densa: m² |\n");
    std::printf("|---|---|---:|---:|---:|---:|---:|---:|\n");

    for (const LinhaTabela& l : tabela) {
        std::string cholesky = "não simétrica";
        if (l.simetrica) {
            cholesky = std::to_string(l.posicoes_cholesky);
        }

        std::printf("| %s | %s | %d | %d | %d | %ld | %s | %ld |\n",
                    l.nome.c_str(),
                    l.ordem.c_str(),
                    l.m,
                    l.bI,
                    l.bS,
                    l.posicoes_lu,
                    cholesky.c_str(),
                    l.posicoes_densa);
    }
}

static bool gravar_csv(const std::vector<LinhaTabela>& tabela) {
    FILE* csv = std::fopen("results/parte3.csv", "w");
    if (csv == nullptr) {
        return false;
    }

    std::fprintf(csv, "matriz,ordem,m,bI,bS,simetrica,posicoes_lu,posicoes_cholesky,posicoes_densa\n");
    for (const LinhaTabela& l : tabela) {
        std::fprintf(csv, "%s,%s,%d,%d,%d,%d,%ld,%ld,%ld\n",
                     l.nome.c_str(),
                     l.ordem.c_str(),
                     l.m,
                     l.bI,
                     l.bS,
                     l.simetrica,
                     l.posicoes_lu,
                     l.posicoes_cholesky,
                     l.posicoes_densa);
    }
    std::fclose(csv);
    return true;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "uso: parte3 arquivo.mtx [arquivo.mtx ...]\n");
        return 1;
    }

    std::filesystem::create_directories("results");

    std::vector<LinhaTabela> tabela;

    try {
        for (int k = 1; k < argc; k++) {
            MtxFile arquivo = read_mtx(argv[k]);
            const DenseMatrix& A = arquivo.A;

            std::vector<int> perm = calcular_rcm(A);
            if (!eh_permutacao(perm)) {
                throw std::runtime_error(arquivo.name + ": o RCM não devolveu uma permutação válida");
            }

            DenseMatrix permutada = permutar_matriz(A, perm);

            LinhaTabela antes = medir(arquivo.name, "original", A);
            LinhaTabela depois = medir(arquivo.name, "RCM", permutada);
            tabela.push_back(antes);
            tabela.push_back(depois);

            if (contar_nao_nulos(A) != contar_nao_nulos(permutada)) {
                throw std::runtime_error(arquivo.name + ": o número de não nulos mudou ao permutar");
            }

            double residuo = residuo_do_sistema_reordenado(A, permutada, perm);
            std::fprintf(stderr, "%s: ||(PAP^T)(Px) - Pb|| relativo = %.1e\n", arquivo.name.c_str(), residuo);

            if (!gravar_permutacao(arquivo.name, perm)) {
                throw std::runtime_error(arquivo.name + ": não foi possível gravar a permutação");
            }
        }
    } catch (const std::exception& erro) {
        std::fprintf(stderr, "erro: %s\n", erro.what());
        return 1;
    }

    imprimir_tabela(tabela);

    if (!gravar_csv(tabela)) {
        std::fprintf(stderr, "erro: não foi possível gravar results/parte3.csv\n");
        return 1;
    }
    std::fprintf(stderr, "\ntabela gravada em results/parte3.csv\n");

    return 0;
}
