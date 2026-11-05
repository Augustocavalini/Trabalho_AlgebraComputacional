// Parte 3: reordenamento RCM e fatoração em banda.
//
// Uso: parte3 arquivo.mtx [arquivo.mtx ...]
//
// Para cada matriz:
// - calcula a permutação RCM e forma PAP^T;
// - confere que o sistema reordenado (PAP^T)(Px) = Pb é equivalente ao original;
// - imprime as bandas e o armazenamento em banda antes e depois;
// - grava a permutação em results/rcm_<matriz>.txt, usada por scripts/spy_rcm.py;
// - resolve o sistema com LU (e com Cholesky, se a matriz for SPD) nas formas densa e em banda,
//   sem e com RCM, e compara os tempos entre si e com a previsão das contagens de operações.
//
// As tabelas também são gravadas em results/parte3.csv e results/parte3_tempos.csv.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

#include "matrix.hpp"
#include "parte2.hpp"
#include "rcm.hpp"

// O enunciado pede a mediana de pelo menos 5 execuções.
const int REPETICOES = 5;

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

// Uma linha da comparação entre as formas densa e em banda: um método (LU ou Cholesky)
// aplicado a uma matriz em uma ordem (original ou RCM). As duas formas são implementações próprias.
struct ComparacaoTempos {
    std::string nome;
    std::string metodo;       // LU ou Cholesky
    std::string ordem;        // original (A) ou RCM (PAP^T)
    std::string falha_banda;  // motivo, se a fatoração em banda não foi possível
    int m = 0;
    int bI = 0;
    int bS = 0;
    double tempo_rcm = 0;
    double fatoracao_densa = 0;
    double fatoracao_banda = 0;
    double solucao_densa = 0;
    double solucao_banda = 0;
    double flops_densa_formula = 0;  // contagens do enunciado
    double flops_banda_formula = 0;
    double flops_densa_exata = 0;    // operações que os nossos laços realmente fazem
    double flops_banda_exata = 0;
    double retroacao_densa = 0;
    double retroacao_banda = 0;
    double direto_densa = 0;
    double direto_banda = 0;
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

// ||A||_inf: a maior soma dos valores absolutos de uma linha.
static double norma_inf_matriz(const DenseMatrix& A) {
    int m = A.m;
    std::vector<double> soma_da_linha(m, 0.0);
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < m; i++) {
            soma_da_linha[i] += std::abs(A(i, j));
        }
    }
    return norma_inf_vetor(soma_da_linha);
}

// Calcula os erros de uma solução x do sistema original Ax = b, com solução exata (1, ..., 1).
static void calcular_erros(const DenseMatrix& A, const std::vector<double>& b, const std::vector<double>& x, double& retroacao, double& direto) {
    int m = A.m;
    std::vector<double> Ax = multiplicar(A, x);

    std::vector<double> residuo(m);
    std::vector<double> erro_de_x(m);
    for (int i = 0; i < m; i++) {
        residuo[i] = b[i] - Ax[i];
        erro_de_x[i] = x[i] - 1.0;
    }

    retroacao = norma_inf_vetor(residuo) / (norma_inf_matriz(A) * norma_inf_vetor(x));
    direto = norma_inf_vetor(erro_de_x);
}

// Uma matriz simétrica é SPD se o nosso Cholesky não encontrar pivô não positivo.
static bool eh_spd(const DenseMatrix& A) {
    if (!eh_simetrica(A)) {
        return false;
    }

    DenseMatrix copia = A;
    try {
        fatoracaoCholesky(copia);
    } catch (const std::exception&) {
        return false;
    }
    return true;
}

// Operações que a nossa LU faz com larguras bI e bS (a densa usa bI = bS = m - 1).
// Na etapa k, há r multiplicadores (uma divisão cada) e r * c atualizações
// (uma multiplicação e uma subtração cada). A busca do pivô não é contada.
static double contar_operacoes_lu(int m, int bI, int bS) {
    double total = 0.0;
    for (int k = 0; k < m; k++) {
        double r = std::min(bI, m - 1 - k);
        double c = std::min(bS, m - 1 - k);
        total += r + 2.0 * r * c;
    }
    return total;
}

// Operações que o nosso Cholesky faz com largura b (o denso usa b = m - 1).
// Na etapa k, há uma raiz, r divisões e r(r + 1)/2 atualizações de duas operações cada.
static double contar_operacoes_cholesky(int m, int b) {
    double total = 0.0;
    for (int k = 0; k < m; k++) {
        double r = std::min(b, m - 1 - k);
        total += 1.0 + r + r * (r + 1.0);
    }
    return total;
}

// Mede a LU densa (com pivoteamento) e a LU em banda (sem pivoteamento) no sistema A_sistema x = b_sistema.
static void medir_lu(const DenseMatrix& A_sistema, const std::vector<double>& b_sistema, ComparacaoTempos& c,
                     std::vector<double>& x_densa, std::vector<double>& x_banda) {
    BandMatrix banda_original = montar_banda(A_sistema, c.bI, c.bS);
    std::vector<lapack_int> pivos(c.m);

    std::vector<double> tempos_fatoracao_densa;
    std::vector<double> tempos_solucao_densa;
    std::vector<double> tempos_fatoracao_banda;
    std::vector<double> tempos_solucao_banda;

    // As formas se alternam a cada repetição, para que variações da máquina
    // ao longo do tempo afetem as duas da mesma forma. As cópias ficam fora do cronômetro.
    for (int repeticao = 0; repeticao < REPETICOES; repeticao++) {
        DenseMatrix LU_densa = A_sistema;
        double inicio = agora();
        fatoracaoLu(LU_densa, pivos);
        double meio = agora();
        x_densa = resolver_sistema_lu(LU_densa, pivos, b_sistema);
        double fim = agora();
        tempos_fatoracao_densa.push_back(meio - inicio);
        tempos_solucao_densa.push_back(fim - meio);

        // A LU em banda não pivota e pode encontrar um pivô nulo que a densa evitaria.
        if (c.falha_banda.empty()) {
            try {
                BandMatrix LU_banda = banda_original;
                inicio = agora();
                fatoracaoLuBanda(LU_banda);
                meio = agora();
                x_banda = resolver_sistema_lu_banda(LU_banda, b_sistema);
                fim = agora();
                tempos_fatoracao_banda.push_back(meio - inicio);
                tempos_solucao_banda.push_back(fim - meio);
            } catch (const std::exception& erro) {
                c.falha_banda = erro.what();
            }
        }
    }

    c.fatoracao_densa = mediana(tempos_fatoracao_densa);
    c.solucao_densa = mediana(tempos_solucao_densa);
    if (c.falha_banda.empty()) {
        c.fatoracao_banda = mediana(tempos_fatoracao_banda);
        c.solucao_banda = mediana(tempos_solucao_banda);
    }

    double m = c.m;
    double bI = c.bI;
    double bS = c.bS;
    c.flops_densa_formula = (2.0 / 3.0) * m * m * m;
    c.flops_banda_formula = 2.0 * m * bI * bS;
    c.flops_densa_exata = contar_operacoes_lu(c.m, c.m - 1, c.m - 1);
    c.flops_banda_exata = contar_operacoes_lu(c.m, c.bI, c.bS);
}

// Mede o Cholesky denso e em banda no sistema A_sistema x = b_sistema.
static void medir_cholesky(const DenseMatrix& A_sistema, const std::vector<double>& b_sistema, ComparacaoTempos& c,
                           std::vector<double>& x_densa, std::vector<double>& x_banda) {
    BandMatrix banda_original = montar_banda(A_sistema, c.bI, 0);

    std::vector<double> tempos_fatoracao_densa;
    std::vector<double> tempos_solucao_densa;
    std::vector<double> tempos_fatoracao_banda;
    std::vector<double> tempos_solucao_banda;

    for (int repeticao = 0; repeticao < REPETICOES; repeticao++) {
        DenseMatrix L_densa = A_sistema;
        double inicio = agora();
        fatoracaoCholesky(L_densa);
        double meio = agora();
        x_densa = resolver_sistema_cholesky(L_densa, b_sistema);
        double fim = agora();
        tempos_fatoracao_densa.push_back(meio - inicio);
        tempos_solucao_densa.push_back(fim - meio);

        BandMatrix L_banda = banda_original;
        inicio = agora();
        fatoracaoCholeskyBanda(L_banda);
        meio = agora();
        x_banda = resolver_sistema_cholesky_banda(L_banda, b_sistema);
        fim = agora();
        tempos_fatoracao_banda.push_back(meio - inicio);
        tempos_solucao_banda.push_back(fim - meio);
    }

    c.fatoracao_densa = mediana(tempos_fatoracao_densa);
    c.solucao_densa = mediana(tempos_solucao_densa);
    c.fatoracao_banda = mediana(tempos_fatoracao_banda);
    c.solucao_banda = mediana(tempos_solucao_banda);

    double m = c.m;
    double b = c.bI;
    c.flops_densa_formula = (1.0 / 3.0) * m * m * m;
    c.flops_banda_formula = m * b * b;
    c.flops_densa_exata = contar_operacoes_cholesky(c.m, c.m - 1);
    c.flops_banda_exata = contar_operacoes_cholesky(c.m, c.bI);
}

// Mede um método em uma ordem. Na ordem RCM, o sistema resolvido é (PAP^T)(Px) = Pb, e x é
// recuperado com P^T antes de medir os erros, que são sempre calculados no sistema original Ax = b.
static ComparacaoTempos comparar(const std::string& nome, const std::string& metodo, const std::string& ordem,
                                 const DenseMatrix& A, const std::vector<double>& b,
                                 const DenseMatrix& A_sistema, const std::vector<double>& b_sistema,
                                 const std::vector<int>& perm, double tempo_rcm) {
    ComparacaoTempos c;
    c.nome = nome;
    c.metodo = metodo;
    c.ordem = ordem;
    c.m = A.m;
    c.tempo_rcm = tempo_rcm;
    calcular_bandas(A_sistema, c.bI, c.bS);

    std::fprintf(stderr, "medindo %s, %s, ordem %s\n", nome.c_str(), metodo.c_str(), ordem.c_str());

    std::vector<double> x_densa;
    std::vector<double> x_banda;
    if (metodo == "LU") {
        medir_lu(A_sistema, b_sistema, c, x_densa, x_banda);
    } else {
        medir_cholesky(A_sistema, b_sistema, c, x_densa, x_banda);
    }

    if (ordem == "RCM") {
        x_densa = despermutar_vetor(x_densa, perm);
        if (c.falha_banda.empty()) {
            x_banda = despermutar_vetor(x_banda, perm);
        }
    }

    calcular_erros(A, b, x_densa, c.retroacao_densa, c.direto_densa);
    if (c.falha_banda.empty()) {
        calcular_erros(A, b, x_banda, c.retroacao_banda, c.direto_banda);
    }
    return c;
}

// Quantas vezes a forma em banda é mais barata que a densa, ou 0 se a conta não se aplicar.
static double razao(double numerador, double denominador) {
    if (denominador <= 0.0) {
        return 0.0;
    }
    return numerador / denominador;
}

static void imprimir_tempos(const std::vector<ComparacaoTempos>& comparacoes) {
    std::printf("\nTempos (s), mediana de %d execuções\n", REPETICOES);
    std::printf("| Matriz | Método | Ordem | b_I | b_S | Fatoração densa | Fatoração em banda | Solução densa | Solução em banda |\n");
    std::printf("|---|---|---|---:|---:|---:|---:|---:|---:|\n");
    for (const ComparacaoTempos& c : comparacoes) {
        if (!c.falha_banda.empty()) {
            std::printf("| %s | %s | %s | %d | %d | %.3e | falhou: %s | %.3e | - |\n",
                        c.nome.c_str(), c.metodo.c_str(), c.ordem.c_str(), c.bI, c.bS,
                        c.fatoracao_densa, c.falha_banda.c_str(), c.solucao_densa);
            continue;
        }
        std::printf("| %s | %s | %s | %d | %d | %.3e | %.3e | %.3e | %.3e |\n",
                    c.nome.c_str(), c.metodo.c_str(), c.ordem.c_str(), c.bI, c.bS,
                    c.fatoracao_densa, c.fatoracao_banda, c.solucao_densa, c.solucao_banda);
    }

    std::printf("\nGanho da fatoração em banda sobre a densa\n");
    std::printf("| Matriz | Método | Ordem | Medido | Previsto (fórmula do enunciado) | Previsto (contagem exata) |\n");
    std::printf("|---|---|---|---:|---:|---:|\n");
    for (const ComparacaoTempos& c : comparacoes) {
        if (!c.falha_banda.empty()) {
            continue;
        }
        std::printf("| %s | %s | %s | %.2f | %.2f | %.2f |\n",
                    c.nome.c_str(), c.metodo.c_str(), c.ordem.c_str(),
                    razao(c.fatoracao_densa, c.fatoracao_banda),
                    razao(c.flops_densa_formula, c.flops_banda_formula),
                    razao(c.flops_densa_exata, c.flops_banda_exata));
    }

    std::printf("\nTaxa da fatoração em GFLOP/s, com a contagem exata\n");
    std::printf("| Matriz | Método | Ordem | Densa | Em banda |\n");
    std::printf("|---|---|---|---:|---:|\n");
    for (const ComparacaoTempos& c : comparacoes) {
        if (!c.falha_banda.empty()) {
            continue;
        }
        std::printf("| %s | %s | %s | %.2f | %.2f |\n",
                    c.nome.c_str(), c.metodo.c_str(), c.ordem.c_str(),
                    c.flops_densa_exata / c.fatoracao_densa / 1e9,
                    c.flops_banda_exata / c.fatoracao_banda / 1e9);
    }

    std::printf("\nErros no sistema original Ax = b\n");
    std::printf("| Matriz | Método | Ordem | Retroação densa | Retroação em banda | Direto densa | Direto em banda |\n");
    std::printf("|---|---|---|---:|---:|---:|---:|\n");
    for (const ComparacaoTempos& c : comparacoes) {
        if (!c.falha_banda.empty()) {
            continue;
        }
        std::printf("| %s | %s | %s | %.1e | %.1e | %.1e | %.1e |\n",
                    c.nome.c_str(), c.metodo.c_str(), c.ordem.c_str(),
                    c.retroacao_densa, c.retroacao_banda, c.direto_densa, c.direto_banda);
    }
}

static bool gravar_csv_tempos(const std::vector<ComparacaoTempos>& comparacoes) {
    FILE* csv = std::fopen("results/parte3_tempos.csv", "w");
    if (csv == nullptr) {
        return false;
    }

    std::fprintf(csv, "matriz,metodo,ordem,m,bI,bS,fatoracao_densa_s,fatoracao_banda_s,solucao_densa_s,solucao_banda_s,"
                      "flops_densa_formula,flops_banda_formula,flops_densa_exata,flops_banda_exata,"
                      "retroacao_densa,retroacao_banda,direto_densa,direto_banda,tempo_rcm_s,falha_banda\n");
    for (const ComparacaoTempos& c : comparacoes) {
        std::fprintf(csv, "%s,%s,%s,%d,%d,%d,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%s\n",
                     c.nome.c_str(),
                     c.metodo.c_str(),
                     c.ordem.c_str(),
                     c.m,
                     c.bI,
                     c.bS,
                     c.fatoracao_densa,
                     c.fatoracao_banda,
                     c.solucao_densa,
                     c.solucao_banda,
                     c.flops_densa_formula,
                     c.flops_banda_formula,
                     c.flops_densa_exata,
                     c.flops_banda_exata,
                     c.retroacao_densa,
                     c.retroacao_banda,
                     c.direto_densa,
                     c.direto_banda,
                     c.tempo_rcm,
                     c.falha_banda.c_str());
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
    std::vector<ComparacaoTempos> comparacoes;

    try {
        for (int k = 1; k < argc; k++) {
            MtxFile arquivo = read_mtx(argv[k]);
            const DenseMatrix& A = arquivo.A;

            // O RCM também é medido, porque é um custo a mais de quem decide reordenar.
            std::vector<int> perm;
            std::vector<double> tempos_rcm;
            for (int repeticao = 0; repeticao < REPETICOES; repeticao++) {
                double inicio = agora();
                perm = calcular_rcm(A);
                tempos_rcm.push_back(agora() - inicio);
            }

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

            // b = A*x, com a solução exata x = (1, ..., 1). O sistema reordenado usa Pb.
            std::vector<double> x_exato(A.m, 1.0);
            std::vector<double> b = multiplicar(A, x_exato);
            std::vector<double> Pb = permutar_vetor(b, perm);
            double tempo_rcm = mediana(tempos_rcm);

            // A LU vale para qualquer matriz; o Cholesky, só para as SPD.
            // Ser SPD não depende da ordem, porque PAP^T é congruente a A.
            std::vector<std::string> metodos = {"LU"};
            if (eh_spd(A)) {
                metodos.push_back("Cholesky");
            }

            for (const std::string& metodo : metodos) {
                comparacoes.push_back(comparar(arquivo.name, metodo, "original", A, b, A, b, perm, tempo_rcm));
                comparacoes.push_back(comparar(arquivo.name, metodo, "RCM", A, b, permutada, Pb, perm, tempo_rcm));
            }
        }
    } catch (const std::exception& erro) {
        std::fprintf(stderr, "erro: %s\n", erro.what());
        return 1;
    }

    imprimir_tabela(tabela);
    imprimir_tempos(comparacoes);

    if (!gravar_csv(tabela) || !gravar_csv_tempos(comparacoes)) {
        std::fprintf(stderr, "erro: não foi possível gravar os resultados em results/\n");
        return 1;
    }
    std::fprintf(stderr, "\ntabelas gravadas em results/parte3.csv e results/parte3_tempos.csv\n");

    return 0;
}
