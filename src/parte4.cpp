// Parte 4: comparação de cada implementação com a rotina correspondente da LAPACK.
//
// Uso: parte4 matrizes arquivo.mtx [arquivo.mtx ...]   grava results/parte4_matrizes.csv
//      parte4 escala                                   grava results/parte4_escala.csv
//
// Em todos os testes, a solução exata é x = (1, ..., 1) e b = A*x.
// Cada tempo é a mediana de 5 execuções, e as cópias de A ficam fora do cronômetro.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "matrix.hpp"
#include "parte2.hpp"

// Função da OpenBLAS que informa quantas threads a BLAS está usando.
extern "C" int openblas_get_num_threads(void);

typedef std::vector<double> Vetor;

// O enunciado pede a mediana de pelo menos 5 execuções.
const int REPETICOES = 5;

const double U = 1.1e-16;

// Uma matriz já preparada para os testes.
struct Problema {
    std::string nome;
    DenseMatrix A;
    Vetor b;
    int bI = 0;
    int bS = 0;
    double maior_coeficiente = 0;  // max |a_ij|, usado no fator de crescimento
    bool spd = false;
    bool tridiagonal = false;
    double kappa_inf = 0;          // estimativa de kappa_inf(A)
};

// O que uma execução de um método devolve.
struct Execucao {
    double tempo_fatoracao = 0;
    double tempo_solucao = 0;
    Vetor x;
    double rho = -1;  // fator de crescimento; negativo se não se aplicar
};

// O tipo de uma função que executa um método uma vez.
typedef Execucao (*FuncaoExecucao)(const Problema&);

// Uma linha da tabela de resultados.
struct Registro {
    std::string matriz;
    std::string metodo;
    std::string implementacao;  // implementado ou LAPACK
    std::string falha;          // motivo, se o método não funcionou nesta matriz
    int m = 0;
    int bI = 0;
    int bS = 0;
    double tempo_fatoracao = 0;
    double tempo_solucao = 0;
    double flops = 0;
    double gflops = 0;
    double nossa_sobre_lapack = 0;  // só nas linhas da implementação própria
    double retroacao = 0;
    double direto = 0;
    double limite_kappa_u = 0;
    double rho = -1;
};

// O instante atual, em segundos.
static double agora() {
    std::chrono::duration<double> tempo = std::chrono::steady_clock::now().time_since_epoch();
    return tempo.count();
}

static double mediana(std::vector<double> valores) {
    std::sort(valores.begin(), valores.end());
    int meio = static_cast<int>(valores.size()) / 2;
    return valores[meio];  // o número de repetições é ímpar
}

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

static double norma_inf_vetor(const Vetor& v) {
    double maior = 0.0;
    for (size_t i = 0; i < v.size(); i++) {
        maior = std::max(maior, std::abs(v[i]));
    }
    return maior;
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

// Estima kappa_inf(A) = ||A||_inf * ||A^-1||_inf com dlange, dgetrf e dgecon.
// Usamos a norma infinito porque é nela que os erros são medidos.
static double estimar_kappa_inf(const DenseMatrix& A) {
    int m = A.m;
    DenseMatrix fatores = A;
    std::vector<lapack_int> pivos(m);

    double norma = LAPACKE_dlange(LAPACK_COL_MAJOR, 'I', m, m, fatores.a.data(), m);

    lapack_int info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, m, m, fatores.a.data(), m, pivos.data());
    if (info != 0) {
        return std::numeric_limits<double>::infinity();  // matriz singular
    }

    double rcond = 0.0;
    LAPACKE_dgecon(LAPACK_COL_MAJOR, 'I', m, fatores.a.data(), m, norma, &rcond);
    if (rcond > 0.0) {
        return 1.0 / rcond;
    }
    return std::numeric_limits<double>::infinity();
}

static Problema preparar_problema(const std::string& nome, const DenseMatrix& A) {
    Problema p;
    p.nome = nome;
    p.A = A;

    Vetor x_exato(A.m, 1.0);
    p.b = multiplicar(A, x_exato);

    calcular_bandas(A, p.bI, p.bS);
    p.tridiagonal = (p.bI <= 1 && p.bS <= 1);
    p.kappa_inf = estimar_kappa_inf(A);

    for (size_t k = 0; k < A.a.size(); k++) {
        p.maior_coeficiente = std::max(p.maior_coeficiente, std::abs(A.a[k]));
    }

    // Uma matriz simétrica é SPD se o nosso Cholesky não encontrar pivô não positivo.
    if (eh_simetrica(A)) {
        DenseMatrix copia = A;
        try {
            fatoracaoCholesky(copia);
            p.spd = true;
        } catch (const std::exception&) {
            p.spd = false;
        }
    }

    return p;
}

// rho = max |u_ij| / max |a_ij|, com U no triângulo superior da matriz fatorada.
static double crescimento_densa(const DenseMatrix& LU, double maior_coeficiente) {
    int m = LU.m;
    double maior_u = 0.0;
    for (int j = 0; j < m; j++) {
        for (int i = 0; i <= j; i++) {
            maior_u = std::max(maior_u, std::abs(LU(i, j)));
        }
    }
    return maior_u / maior_coeficiente;
}

// O mesmo cálculo em banda: U vai da diagonal até a banda superior do arranjo.
// No dgbtrf, essa banda já inclui o crescimento causado pelo pivoteamento.
static double crescimento_banda(const BandMatrix& LU, double maior_coeficiente) {
    int m = LU.m;
    double maior_u = 0.0;
    for (int j = 0; j < m; j++) {
        int primeira_linha = std::max(0, j - LU.bS);
        for (int i = primeira_linha; i <= j; i++) {
            maior_u = std::max(maior_u, std::abs(LU(i, j)));
        }
    }
    return maior_u / maior_coeficiente;
}

// As funções abaixo executam um método uma vez. Nas rotinas da LAPACK, a cópia de b entra no
// tempo da solução, como acontece dentro das nossas funções resolver_sistema_*.

static Execucao lu_densa_nossa(const Problema& p) {
    Execucao e;
    DenseMatrix LU = p.A;
    std::vector<lapack_int> pivos(p.A.m);

    double inicio = agora();
    fatoracaoLu(LU, pivos);
    double meio = agora();
    e.x = resolver_sistema_lu(LU, pivos, p.b);
    double fim = agora();

    e.tempo_fatoracao = meio - inicio;
    e.tempo_solucao = fim - meio;
    e.rho = crescimento_densa(LU, p.maior_coeficiente);
    return e;
}

static Execucao lu_densa_lapack(const Problema& p) {
    Execucao e;
    int m = p.A.m;
    DenseMatrix LU = p.A;
    std::vector<lapack_int> pivos(m);

    double inicio = agora();
    fatoracaoLULapacke(LU, pivos);
    double meio = agora();
    e.x = p.b;
    LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', m, 1, LU.a.data(), m, pivos.data(), e.x.data(), m);
    double fim = agora();

    e.tempo_fatoracao = meio - inicio;
    e.tempo_solucao = fim - meio;
    e.rho = crescimento_densa(LU, p.maior_coeficiente);
    return e;
}

static Execucao cholesky_densa_nossa(const Problema& p) {
    Execucao e;
    DenseMatrix L = p.A;

    double inicio = agora();
    fatoracaoCholesky(L);
    double meio = agora();
    e.x = resolver_sistema_cholesky(L, p.b);
    double fim = agora();

    e.tempo_fatoracao = meio - inicio;
    e.tempo_solucao = fim - meio;
    return e;
}

static Execucao cholesky_densa_lapack(const Problema& p) {
    Execucao e;
    int m = p.A.m;
    DenseMatrix L = p.A;

    double inicio = agora();
    fatoracaoCholeskyLapacke(L);
    double meio = agora();
    e.x = p.b;
    LAPACKE_dpotrs(LAPACK_COL_MAJOR, 'L', m, 1, L.a.data(), m, e.x.data(), m);
    double fim = agora();

    e.tempo_fatoracao = meio - inicio;
    e.tempo_solucao = fim - meio;
    return e;
}

static Execucao lu_banda_nossa(const Problema& p) {
    Execucao e;
    BandMatrix LU = montar_banda(p.A, p.bI, p.bS);

    double inicio = agora();
    fatoracaoLuBanda(LU);
    double meio = agora();
    e.x = resolver_sistema_lu_banda(LU, p.b);
    double fim = agora();

    e.tempo_fatoracao = meio - inicio;
    e.tempo_solucao = fim - meio;
    e.rho = crescimento_banda(LU, p.maior_coeficiente);
    return e;
}

static Execucao lu_banda_lapack(const Problema& p) {
    Execucao e;
    int m = p.A.m;

    // O dgbtrf pivota e por isso precisa de bI posições extras: 2bI + bS + 1 linhas.
    BandMatrix LU = montar_banda(p.A, p.bI, p.bI + p.bS);
    std::vector<lapack_int> pivos(m);

    double inicio = agora();
    fatoracaoLuBandaLapacke(LU, pivos);
    double meio = agora();
    e.x = p.b;
    LAPACKE_dgbtrs(LAPACK_COL_MAJOR, 'N', m, p.bI, p.bS, 1, LU.a.data(), LU.linhas, pivos.data(), e.x.data(), m);
    double fim = agora();

    e.tempo_fatoracao = meio - inicio;
    e.tempo_solucao = fim - meio;
    e.rho = crescimento_banda(LU, p.maior_coeficiente);
    return e;
}

static Execucao cholesky_banda_nossa(const Problema& p) {
    Execucao e;
    BandMatrix L = montar_banda(p.A, p.bI, 0);

    double inicio = agora();
    fatoracaoCholeskyBanda(L);
    double meio = agora();
    e.x = resolver_sistema_cholesky_banda(L, p.b);
    double fim = agora();

    e.tempo_fatoracao = meio - inicio;
    e.tempo_solucao = fim - meio;
    return e;
}

static Execucao cholesky_banda_lapack(const Problema& p) {
    Execucao e;
    int m = p.A.m;
    BandMatrix L = montar_banda(p.A, p.bI, 0);

    double inicio = agora();
    fatoracaoCholeskyBandaLapacke(L);
    double meio = agora();
    e.x = p.b;
    LAPACKE_dpbtrs(LAPACK_COL_MAJOR, 'L', m, p.bI, 1, L.a.data(), L.linhas, e.x.data(), m);
    double fim = agora();

    e.tempo_fatoracao = meio - inicio;
    e.tempo_solucao = fim - meio;
    return e;
}

// Separa uma matriz tridiagonal nos três vetores.
static void extrair_tridiagonal(const DenseMatrix& A, Vetor& inferior, Vetor& diagonal, Vetor& superior) {
    int m = A.m;
    inferior.assign(m - 1, 0.0);
    diagonal.assign(m, 0.0);
    superior.assign(m - 1, 0.0);

    for (int k = 0; k < m; k++) {
        diagonal[k] = A(k, k);
    }
    for (int k = 0; k < m - 1; k++) {
        inferior[k] = A(k + 1, k);
        superior[k] = A(k, k + 1);
    }
}

static Execucao thomas_nossa(const Problema& p) {
    Execucao e;
    Vetor inferior;
    Vetor diagonal;
    Vetor superior;
    extrair_tridiagonal(p.A, inferior, diagonal, superior);

    double inicio = agora();
    fatoracaoThomas(inferior, diagonal, superior);
    double meio = agora();
    e.x = resolver_sistema_thomas(inferior, diagonal, superior, p.b);
    double fim = agora();

    e.tempo_fatoracao = meio - inicio;
    e.tempo_solucao = fim - meio;
    return e;
}

// O dgtsv faz fatoração e solução juntas, então todo o tempo fica na fatoração.
static Execucao thomas_lapack(const Problema& p) {
    Execucao e;
    Vetor inferior;
    Vetor diagonal;
    Vetor superior;
    extrair_tridiagonal(p.A, inferior, diagonal, superior);

    double inicio = agora();
    e.x = resolver_tridiagonal_lapacke(inferior, diagonal, superior, p.b);
    double fim = agora();

    e.tempo_fatoracao = fim - inicio;
    e.tempo_solucao = 0.0;
    return e;
}

// Preenche um registro com as medianas e os erros de uma série de execuções.
static void preencher_registro(Registro& r, const Problema& p, const std::vector<double>& fatoracoes, const std::vector<double>& solucoes, const Execucao& ultima, bool flops_inclui_solucao) {
    int m = p.A.m;

    r.tempo_fatoracao = mediana(fatoracoes);
    r.tempo_solucao = mediana(solucoes);
    r.rho = ultima.rho;

    double tempo_da_taxa = r.tempo_fatoracao;
    if (flops_inclui_solucao) {
        tempo_da_taxa = r.tempo_fatoracao + r.tempo_solucao;
    }
    r.gflops = r.flops / tempo_da_taxa / 1e9;

    // Erro de retroação: ||b - A*x||_inf / (||A||_inf * ||x||_inf).
    Vetor Ax = multiplicar(p.A, ultima.x);
    Vetor residuo(m);
    Vetor erro_de_x(m);
    Vetor soma_da_linha(m, 0.0);
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < m; i++) {
            soma_da_linha[i] += std::abs(p.A(i, j));
        }
    }
    for (int i = 0; i < m; i++) {
        residuo[i] = p.b[i] - Ax[i];
        erro_de_x[i] = ultima.x[i] - 1.0;
    }
    r.retroacao = norma_inf_vetor(residuo) / (norma_inf_vetor(soma_da_linha) * norma_inf_vetor(ultima.x));

    // Erro direto: como x = (1, ..., 1), temos ||x||_inf = 1.
    r.direto = norma_inf_vetor(erro_de_x);
    r.limite_kappa_u = p.kappa_inf * U;
}

// Mede um método nas duas implementações, alternando-as a cada repetição para que
// variações da máquina ao longo do tempo afetem as duas da mesma forma.
static void medir_metodo(const Problema& p, const std::string& metodo, double flops, bool flops_inclui_solucao,
                         FuncaoExecucao nossa, FuncaoExecucao lapack, std::vector<Registro>& tabela) {
    Registro r_nossa;
    Registro r_lapack;
    r_nossa.implementacao = "implementado";
    r_lapack.implementacao = "LAPACK";

    Registro* registros[2] = {&r_nossa, &r_lapack};
    for (int k = 0; k < 2; k++) {
        registros[k]->matriz = p.nome;
        registros[k]->metodo = metodo;
        registros[k]->m = p.A.m;
        registros[k]->bI = p.bI;
        registros[k]->bS = p.bS;
        registros[k]->flops = flops;
    }

    std::vector<double> fatoracoes_nossa;
    std::vector<double> solucoes_nossa;
    std::vector<double> fatoracoes_lapack;
    std::vector<double> solucoes_lapack;
    Execucao ultima_nossa;
    Execucao ultima_lapack;

    for (int repeticao = 0; repeticao < REPETICOES; repeticao++) {
        // Uma implementação pode falhar sem que a outra falhe: a LU em banda
        // não pivota e pode encontrar um pivô nulo que o dgbtrf evitaria.
        if (r_nossa.falha.empty()) {
            try {
                ultima_nossa = nossa(p);
                fatoracoes_nossa.push_back(ultima_nossa.tempo_fatoracao);
                solucoes_nossa.push_back(ultima_nossa.tempo_solucao);
            } catch (const std::exception& erro) {
                r_nossa.falha = erro.what();
            }
        }
        if (r_lapack.falha.empty()) {
            try {
                ultima_lapack = lapack(p);
                fatoracoes_lapack.push_back(ultima_lapack.tempo_fatoracao);
                solucoes_lapack.push_back(ultima_lapack.tempo_solucao);
            } catch (const std::exception& erro) {
                r_lapack.falha = erro.what();
            }
        }
    }

    if (r_nossa.falha.empty()) {
        preencher_registro(r_nossa, p, fatoracoes_nossa, solucoes_nossa, ultima_nossa, flops_inclui_solucao);
    }
    if (r_lapack.falha.empty()) {
        preencher_registro(r_lapack, p, fatoracoes_lapack, solucoes_lapack, ultima_lapack, flops_inclui_solucao);
    }

    // No Thomas, o dgtsv só informa o tempo total, então a razão usa o total nos dois lados.
    if (r_nossa.falha.empty() && r_lapack.falha.empty()) {
        double tempo_nossa = r_nossa.tempo_fatoracao;
        double tempo_lapack = r_lapack.tempo_fatoracao;
        if (flops_inclui_solucao) {
            tempo_nossa = r_nossa.tempo_fatoracao + r_nossa.tempo_solucao;
            tempo_lapack = r_lapack.tempo_fatoracao + r_lapack.tempo_solucao;
        }
        r_nossa.nossa_sobre_lapack = tempo_nossa / tempo_lapack;
    }

    tabela.push_back(r_nossa);
    tabela.push_back(r_lapack);
}

static void imprimir_cabecalho() {
    std::printf("\n| Matriz | Método | Impl. | m | b_I | b_S | Fatoração (s) | Solução (s) | GFLOP/s | Implementado/LAPACK | Retroação | Direto | κ∞(A)u | Direto ≤ κu | ρ |\n");
    std::printf("|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---:|---:|\n");
}

static void imprimir_registro(const Registro& r) {
    if (!r.falha.empty()) {
        std::printf("| %s | %s | %s | %d | %d | %d | falhou: %s |\n",
                    r.matriz.c_str(), r.metodo.c_str(), r.implementacao.c_str(), r.m, r.bI, r.bS, r.falha.c_str());
        std::fflush(stdout);
        return;
    }

    char texto_razao[32] = "-";
    if (r.nossa_sobre_lapack > 0.0) {
        std::snprintf(texto_razao, sizeof(texto_razao), "%.2f", r.nossa_sobre_lapack);
    }

    char texto_rho[32] = "-";
    if (r.rho >= 0.0) {
        std::snprintf(texto_rho, sizeof(texto_rho), "%.2f", r.rho);
    }

    // O dgtsv faz fatoração e solução juntas, e isso aparece na coluna da solução.
    char texto_solucao[32] = "junto";
    if (r.tempo_solucao > 0.0) {
        std::snprintf(texto_solucao, sizeof(texto_solucao), "%.3e", r.tempo_solucao);
    }

    char texto_gflops[32] = "-";
    if (r.flops > 0.0) {
        std::snprintf(texto_gflops, sizeof(texto_gflops), "%.2f", r.gflops);
    }

    const char* dentro_do_limite = "não";
    if (r.direto <= r.limite_kappa_u) {
        dentro_do_limite = "sim";
    }

    std::printf("| %s | %s | %s | %d | %d | %d | %.3e | %s | %s | %s | %.1e | %.1e | %.1e | %s | %s |\n",
                r.matriz.c_str(),
                r.metodo.c_str(),
                r.implementacao.c_str(),
                r.m,
                r.bI,
                r.bS,
                r.tempo_fatoracao,
                texto_solucao,
                texto_gflops,
                texto_razao,
                r.retroacao,
                r.direto,
                r.limite_kappa_u,
                dentro_do_limite,
                texto_rho);
    std::fflush(stdout);  // mostra cada linha assim que ela fica pronta
}

static bool gravar_csv(const std::string& caminho, const std::vector<Registro>& tabela) {
    std::filesystem::create_directories("results");

    FILE* csv = std::fopen(caminho.c_str(), "w");
    if (csv == nullptr) {
        return false;
    }

    std::fprintf(csv, "matriz,metodo,implementacao,m,bI,bS,fatoracao_s,solucao_s,flops,gflops,implementado_sobre_lapack,retroacao,direto,limite_kappa_u,rho,falha\n");
    for (const Registro& r : tabela) {
        // O nome da matriz fica entre aspas porque pode ter vírgulas, como em tridiag(-1,2,-1).
        std::fprintf(csv, "\"%s\",%s,%s,%d,%d,%d,%.6e,%.6e,%.6e,%.4f,%.4f,%.6e,%.6e,%.6e,%.4f,%s\n",
                     r.matriz.c_str(),
                     r.metodo.c_str(),
                     r.implementacao.c_str(),
                     r.m,
                     r.bI,
                     r.bS,
                     r.tempo_fatoracao,
                     r.tempo_solucao,
                     r.flops,
                     r.gflops,
                     r.nossa_sobre_lapack,
                     r.retroacao,
                     r.direto,
                     r.limite_kappa_u,
                     r.rho,
                     r.falha.c_str());
    }
    std::fclose(csv);
    return true;
}

// Imprime as linhas novas da tabela, a partir da posição indicada.
static void imprimir_novos(const std::vector<Registro>& tabela, size_t primeira) {
    for (size_t k = primeira; k < tabela.size(); k++) {
        imprimir_registro(tabela[k]);
    }
}

// Roda todos os métodos que se aplicam a uma matriz, com as contagens de operações do enunciado.
static void comparar_matriz(const Problema& p, std::vector<Registro>& tabela) {
    double m = p.A.m;
    double bI = p.bI;
    double bS = p.bS;
    size_t primeira = tabela.size();

    medir_metodo(p, "LU densa", (2.0 / 3.0) * m * m * m, false, lu_densa_nossa, lu_densa_lapack, tabela);
    medir_metodo(p, "LU em banda", 2.0 * m * bI * bS, false, lu_banda_nossa, lu_banda_lapack, tabela);

    if (p.spd) {
        medir_metodo(p, "Cholesky densa", (1.0 / 3.0) * m * m * m, false, cholesky_densa_nossa, cholesky_densa_lapack, tabela);
        medir_metodo(p, "Cholesky em banda", m * bI * bI, false, cholesky_banda_nossa, cholesky_banda_lapack, tabela);
    }

    if (p.tridiagonal) {
        medir_metodo(p, "Thomas", 8.0 * m, true, thomas_nossa, thomas_lapack, tabela);
    }

    imprimir_novos(tabela, primeira);
}

// Modo "matrizes": todos os métodos aplicáveis em cada arquivo.
static int rodar_matrizes(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "uso: parte4 matrizes arquivo.mtx [arquivo.mtx ...]\n");
        return 1;
    }

    std::vector<Registro> tabela;
    imprimir_cabecalho();

    for (int k = 2; k < argc; k++) {
        MtxFile arquivo = read_mtx(argv[k]);
        Problema p = preparar_problema(arquivo.name, arquivo.A);
        comparar_matriz(p, tabela);
    }

    if (!gravar_csv("results/parte4_matrizes.csv", tabela)) {
        std::fprintf(stderr, "erro: não foi possível gravar results/parte4_matrizes.csv\n");
        return 1;
    }
    std::fprintf(stderr, "\ntabela gravada em results/parte4_matrizes.csv\n");
    return 0;
}

// Matriz densa com coeficientes uniformes entre -1 e 1. A semente fixa torna o estudo reproduzível.
static DenseMatrix densa_aleatoria(int m) {
    std::mt19937 gerador(12345 + m);
    std::uniform_real_distribution<double> distribuicao(-1.0, 1.0);

    DenseMatrix A(m);
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < m; i++) {
            A(i, j) = distribuicao(gerador);
        }
    }
    return A;
}

// Matriz densa SPD: simétrica e aleatória, com m somado à diagonal. Como cada linha soma
// menos de m em valor absoluto fora da diagonal, ela é estritamente dominante e, portanto, SPD.
static DenseMatrix densa_spd_aleatoria(int m) {
    std::mt19937 gerador(54321 + m);
    std::uniform_real_distribution<double> distribuicao(-1.0, 1.0);

    DenseMatrix A(m);
    for (int j = 0; j < m; j++) {
        for (int i = j; i < m; i++) {
            double valor = distribuicao(gerador);
            A(i, j) = valor;
            A(j, i) = valor;
        }
    }
    for (int k = 0; k < m; k++) {
        A(k, k) += m;
    }
    return A;
}

// Estudo de escala da tridiag(-1, 2, -1). Com m = 10^6 a forma densa não caberia na memória,
// por isso esta medição só usa os três vetores.
static void medir_tridiagonal(int m, std::vector<Registro>& tabela) {
    Vetor inferior_original(m - 1, -1.0);
    Vetor diagonal_original(m, 2.0);
    Vetor superior_original(m - 1, -1.0);

    // b = A*x com x = (1, ..., 1): só a primeira e a última equações não se cancelam.
    Vetor b(m, 0.0);
    b[0] = 1.0;
    b[m - 1] = 1.0;

    // kappa_inf estimado com dgttrf e dgtcon. A norma infinito desta matriz é 4.
    double norma_inf = 4.0;
    Vetor inf_cond = inferior_original;
    Vetor diag_cond = diagonal_original;
    Vetor sup_cond = superior_original;
    Vetor sup2_cond(m - 2);
    std::vector<lapack_int> pivos(m);
    LAPACKE_dgttrf(m, inf_cond.data(), diag_cond.data(), sup_cond.data(), sup2_cond.data(), pivos.data());
    double rcond = 0.0;
    LAPACKE_dgtcon('I', m, inf_cond.data(), diag_cond.data(), sup_cond.data(), sup2_cond.data(), pivos.data(), norma_inf, &rcond);

    std::vector<double> fatoracoes_nossa;
    std::vector<double> solucoes_nossa;
    std::vector<double> totais_lapack;
    Vetor x_nossa;
    Vetor x_lapack;

    for (int repeticao = 0; repeticao < REPETICOES; repeticao++) {
        Vetor inferior = inferior_original;
        Vetor diagonal = diagonal_original;
        Vetor superior = superior_original;

        double inicio = agora();
        fatoracaoThomas(inferior, diagonal, superior);
        double meio = agora();
        x_nossa = resolver_sistema_thomas(inferior, diagonal, superior, b);
        double fim = agora();
        fatoracoes_nossa.push_back(meio - inicio);
        solucoes_nossa.push_back(fim - meio);

        inferior = inferior_original;
        diagonal = diagonal_original;
        superior = superior_original;

        inicio = agora();
        x_lapack = resolver_tridiagonal_lapacke(inferior, diagonal, superior, b);
        fim = agora();
        totais_lapack.push_back(fim - inicio);
    }

    const Vetor* solucoes[2] = {&x_nossa, &x_lapack};
    const char* implementacoes[2] = {"implementado", "LAPACK"};
    double tempos_fatoracao[2] = {mediana(fatoracoes_nossa), mediana(totais_lapack)};
    double tempos_solucao[2] = {mediana(solucoes_nossa), 0.0};

    for (int k = 0; k < 2; k++) {
        const Vetor& x = *solucoes[k];

        // Multiplica pela matriz tridiagonal e compara com b.
        double maior_residuo = 0.0;
        for (int i = 0; i < m; i++) {
            double Ax_i = 2.0 * x[i];
            if (i > 0) {
                Ax_i -= x[i - 1];
            }
            if (i < m - 1) {
                Ax_i -= x[i + 1];
            }
            maior_residuo = std::max(maior_residuo, std::abs(b[i] - Ax_i));
        }

        double maior_erro = 0.0;
        for (int i = 0; i < m; i++) {
            maior_erro = std::max(maior_erro, std::abs(x[i] - 1.0));
        }

        Registro r;
        r.matriz = "tridiag(-1,2,-1)";
        r.metodo = "Thomas";
        r.implementacao = implementacoes[k];
        r.m = m;
        r.bI = 1;
        r.bS = 1;
        r.tempo_fatoracao = tempos_fatoracao[k];
        r.tempo_solucao = tempos_solucao[k];
        r.flops = 8.0 * m;
        r.gflops = r.flops / (r.tempo_fatoracao + r.tempo_solucao) / 1e9;
        r.retroacao = maior_residuo / (norma_inf * norma_inf_vetor(x));
        r.direto = maior_erro;
        r.limite_kappa_u = U / rcond;
        if (k == 0) {
            r.nossa_sobre_lapack = (tempos_fatoracao[0] + tempos_solucao[0]) / tempos_fatoracao[1];
        }
        tabela.push_back(r);
    }
}

// Modo "escala": tridiagonal até m = 10^6 e densas aleatórias de m = 250 a 4000.
static int rodar_escala() {
    std::vector<Registro> tabela;
    imprimir_cabecalho();

    std::vector<int> dimensoes_tridiagonal = {1000, 10000, 100000, 1000000};
    for (size_t k = 0; k < dimensoes_tridiagonal.size(); k++) {
        size_t primeira = tabela.size();
        medir_tridiagonal(dimensoes_tridiagonal[k], tabela);
        imprimir_novos(tabela, primeira);
    }

    std::vector<int> dimensoes_densas = {250, 500, 1000, 2000, 4000};
    for (size_t k = 0; k < dimensoes_densas.size(); k++) {
        int inteiro_m = dimensoes_densas[k];
        double m = inteiro_m;
        size_t primeira = tabela.size();

        Problema aleatoria = preparar_problema("densa aleatória", densa_aleatoria(inteiro_m));
        medir_metodo(aleatoria, "LU densa", (2.0 / 3.0) * m * m * m, false, lu_densa_nossa, lu_densa_lapack, tabela);

        Problema spd = preparar_problema("densa SPD aleatória", densa_spd_aleatoria(inteiro_m));
        medir_metodo(spd, "Cholesky densa", (1.0 / 3.0) * m * m * m, false, cholesky_densa_nossa, cholesky_densa_lapack, tabela);

        imprimir_novos(tabela, primeira);
    }

    if (!gravar_csv("results/parte4_escala.csv", tabela)) {
        std::fprintf(stderr, "erro: não foi possível gravar results/parte4_escala.csv\n");
        return 1;
    }
    std::fprintf(stderr, "\ntabela gravada em results/parte4_escala.csv\n");
    return 0;
}

int main(int argc, char** argv) {
    int threads = openblas_get_num_threads();
    std::fprintf(stderr, "BLAS com %d thread(s)\n", threads);
    if (threads != 1) {
        std::fprintf(stderr, "aviso: rode com OPENBLAS_NUM_THREADS=1 para comparar códigos seriais\n");
    }

    // Por padrão, o LAPACKE percorre toda a matriz procurando NaN antes de cada chamada. Isso não
    // faz parte da rotina da LAPACK e pesaria nas rotinas baratas, como o dgtsv, então desligamos.
    LAPACKE_set_nancheck(0);

    std::string modo;
    if (argc > 1) {
        modo = argv[1];
    }

    try {
        if (modo == "matrizes") {
            return rodar_matrizes(argc, argv);
        }
        if (modo == "escala") {
            return rodar_escala();
        }
    } catch (const std::exception& erro) {
        std::fprintf(stderr, "erro: %s\n", erro.what());
        return 1;
    }

    std::fprintf(stderr, "uso: parte4 matrizes arquivo.mtx [arquivo.mtx ...]\n");
    std::fprintf(stderr, "     parte4 escala\n");
    return 1;
}
