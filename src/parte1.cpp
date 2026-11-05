// Parte 1: leitura dos .mtx, montagem da forma densa e caracterização das matrizes.
//
// Uso: parte1 [arquivo.mtx ...]
// Sem argumentos, usa todos os arquivos matrices/*.mtx.
// Imprime a tabela em Markdown e grava results/parte1.csv.
#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <limits>
#include <string>
#include <vector>

#include <lapacke.h>

#include "matrix.hpp"

// Uma linha da tabela da Parte 1.
struct Caracterizacao {
    std::string nome;
    int m = 0;
    int nnz = 0;
    double densidade = 0;
    bool simetrica = false;
    bool spd = false;
    int bI = 0;               // largura de banda inferior
    int bS = 0;               // largura de banda superior
    double norma1 = 0;        // ||A||_1, calculada com dlange
    double rcond = 0;         // estimativa de 1/kappa_1, calculada com dgecon ou dpocon
    double kappa1 = 0;        // estimativa de kappa_1 = 1/rcond
    std::string rotina = "";  // rotina usada na estimativa: dpocon ou dgecon
};

// Percorre todos os coeficientes para calcular nnz, densidade, bandas e simetria.
static void medir_estrutura(const std::vector<std::vector<double>>& matriz, Caracterizacao& c) {
    int m = c.m;
    c.simetrica = true;  // vale até encontrarmos um par diferente

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            double valor = matriz[i][j];

            if (valor != 0.0) {
                c.nnz = c.nnz + 1;

                // Abaixo da diagonal: a distância até ela é i - j.
                if (i > j) {
                    c.bI = std::max(c.bI, i - j);
                }
                // Acima da diagonal: a distância até ela é j - i.
                if (j > i) {
                    c.bS = std::max(c.bS, j - i);
                }
            }

            // Basta olhar o triângulo inferior: cada par (i,j), (j,i) é comparado uma vez.
            if (i > j) {
                double espelhado = matriz[j][i];
                if (valor != espelhado) {
                    c.simetrica = false;
                }
            }
        }
    }

    double total_de_coeficientes = static_cast<double>(m) * static_cast<double>(m);
    c.densidade = static_cast<double>(c.nnz) / total_de_coeficientes;
}

// Estima kappa_1(A) = ||A||_1 * ||A^-1||_1. O dpocon e o dgecon não recebem a matriz
// original: eles precisam dos fatores de Cholesky ou de LU e de ||A||_1.
//
// Usa o dpocon se a matriz for SPD e o dgecon nos outros casos.
static void estimar_condicionamento(const DenseMatrix& A, Caracterizacao& c) {
    int m = A.m;

    // As fatorações sobrescrevem o vetor recebido, então trabalhamos em uma cópia.
    std::vector<double> fatores = A.a;

    c.norma1 = LAPACKE_dlange(LAPACK_COL_MAJOR, '1', m, m, fatores.data(), m);

    // Uma matriz simétrica é SPD se e somente se o Cholesky não encontra pivô não positivo.
    if (c.simetrica) {
        lapack_int info_cholesky = LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', m, fatores.data(), m);

        if (info_cholesky == 0) {
            c.spd = true;
            c.rotina = "dpocon";
            LAPACKE_dpocon(LAPACK_COL_MAJOR, 'L', m, fatores.data(), m, c.norma1, &c.rcond);
        } else {
            // O Cholesky falhou no meio e deixou a cópia alterada: voltamos ao original.
            fatores = A.a;
        }
    }

    if (!c.spd) {
        c.rotina = "dgecon";

        std::vector<lapack_int> pivos(m);
        lapack_int info_lu = LAPACKE_dgetrf(LAPACK_COL_MAJOR, m, m, fatores.data(), m, pivos.data());

        if (info_lu == 0) {
            LAPACKE_dgecon(LAPACK_COL_MAJOR, '1', m, fatores.data(), m, c.norma1, &c.rcond);
        } else {
            // info_lu > 0 significa que U tem um zero na diagonal: a matriz é singular.
            c.rcond = 0.0;
        }
    }

    if (c.rcond > 0.0) {
        c.kappa1 = 1.0 / c.rcond;
    } else {
        c.kappa1 = std::numeric_limits<double>::infinity();
    }
}

// Monta a linha da tabela para uma matriz já lida.
static Caracterizacao caracterizar(const MtxFile& arquivo) {
    Caracterizacao c;
    c.nome = arquivo.name;
    c.m = arquivo.A.m;

    medir_estrutura(arquivo.matriz, c);
    estimar_condicionamento(arquivo.A, c);

    return c;
}

// Devolve os arquivos a processar: os da linha de comando ou, se não houver, os de matrices/.
static std::vector<std::string> listar_arquivos(int argc, char** argv) {
    std::vector<std::string> arquivos;

    // argv[0] é o nome do próprio programa; os argumentos começam em argv[1].
    for (int k = 1; k < argc; k++) {
        arquivos.push_back(argv[k]);
    }

    if (arquivos.empty() && std::filesystem::is_directory("matrices")) {
        for (const std::filesystem::directory_entry& item : std::filesystem::directory_iterator("matrices")) {
            if (item.path().extension() == ".mtx") {
                arquivos.push_back(item.path().string());
            }
        }
        std::sort(arquivos.begin(), arquivos.end());
    }

    return arquivos;
}

static const char* sim_ou_nao(bool condicao) {
    if (condicao) {
        return "sim";
    }
    return "não";
}

static void imprimir_tabela(const std::vector<Caracterizacao>& tabela) {
    std::printf("\n| Matriz | m | nnz | nnz/m² | Simétrica | SPD | b_I | b_S | ‖A‖₁ | κ₁(A) estimado | Rotina |\n");
    std::printf("|---|---:|---:|---:|:---:|:---:|---:|---:|---:|---:|---|\n");

    for (const Caracterizacao& c : tabela) {
        std::printf("| %s | %d | %d | %.3e | %s | %s | %d | %d | %.3e | %.2e | %s |\n",
                    c.nome.c_str(),
                    c.m,
                    c.nnz,
                    c.densidade,
                    sim_ou_nao(c.simetrica),
                    sim_ou_nao(c.spd),
                    c.bI,
                    c.bS,
                    c.norma1,
                    c.kappa1,
                    c.rotina.c_str());
    }
}

// Devolve false se não conseguir criar o arquivo.
static bool gravar_csv(const std::vector<Caracterizacao>& tabela) {
    std::filesystem::create_directories("results");

    FILE* csv = std::fopen("results/parte1.csv", "w");
    if (csv == nullptr) {
        return false;
    }

    std::fprintf(csv, "matriz,m,nnz,densidade,simetrica,spd,bI,bS,norma1,rcond,kappa1,rotina\n");

    for (const Caracterizacao& c : tabela) {
        std::fprintf(csv, "%s,%d,%d,%.6e,%d,%d,%d,%d,%.6e,%.6e,%.6e,%s\n",
                     c.nome.c_str(),
                     c.m,
                     c.nnz,
                     c.densidade,
                     c.simetrica,
                     c.spd,
                     c.bI,
                     c.bS,
                     c.norma1,
                     c.rcond,
                     c.kappa1,
                     c.rotina.c_str());
    }

    std::fclose(csv);
    return true;
}

int main(int argc, char** argv) {
    std::vector<std::string> arquivos = listar_arquivos(argc, argv);
    if (arquivos.empty()) {
        std::fprintf(stderr, "nenhum .mtx encontrado; passe os arquivos ou coloque-os em matrices/\n");
        return 1;
    }

    std::vector<Caracterizacao> tabela;

    for (const std::string& caminho : arquivos) {
        // O read_mtx lança uma exceção se houver erro no arquivo; ela é capturada no catch.
        try {
            MtxFile arquivo = read_mtx(caminho);

            std::fprintf(stderr, "%s: %s %s, m=%d, %d entradas armazenadas\n",
                         caminho.c_str(),
                         arquivo.field.c_str(),
                         arquivo.symmetry.c_str(),
                         arquivo.A.m,
                         arquivo.stored_entries);

            Caracterizacao linha = caracterizar(arquivo);
            tabela.push_back(linha);
        } catch (const std::exception& erro) {
            std::fprintf(stderr, "erro: %s\n", erro.what());
            return 1;
        }
    }

    imprimir_tabela(tabela);

    bool csv_gravado = gravar_csv(tabela);
    if (!csv_gravado) {
        std::fprintf(stderr, "erro: não foi possível gravar results/parte1.csv\n");
        return 1;
    }
    std::fprintf(stderr, "\ntabela gravada em results/parte1.csv\n");

    return 0;
}
