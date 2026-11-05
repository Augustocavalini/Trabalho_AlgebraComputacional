#include "matrix.hpp"

#include <cstdio>
#include <stdexcept>

// O mmio é uma biblioteca em C. O extern "C" avisa o compilador de C++ disso,
// para que ele encontre as funções na hora de ligar o programa.
extern "C" {
#include "mmio.h"
}

// Extrai o nome da matriz a partir do caminho do arquivo.
// Exemplo: "matrices/bcsstk14.mtx" vira "bcsstk14".
static std::string base_name(const std::string& path) {
    std::string nome = path;

    // Remove tudo até a última barra, se houver alguma.
    // A busca devolve std::string::npos quando não encontra o que procura.
    size_t posicao_da_barra = nome.find_last_of('/');
    if (posicao_da_barra != std::string::npos) {
        nome = nome.substr(posicao_da_barra + 1);
    }

    // Remove a extensão ".mtx", se houver.
    size_t posicao_da_extensao = nome.rfind(".mtx");
    if (posicao_da_extensao != std::string::npos) {
        nome = nome.substr(0, posicao_da_extensao);
    }

    return nome;
}

// Fecha o arquivo e interrompe a leitura, enviando a mensagem de erro para quem chamou.
static void fechar_e_lancar_erro(FILE* arquivo, const std::string& mensagem) {
    std::fclose(arquivo);
    throw std::runtime_error(mensagem);
}

MtxFile read_mtx(const std::string& path) {
    // 1) Abre o arquivo para leitura.
    FILE* arquivo = std::fopen(path.c_str(), "r");
    if (arquivo == nullptr) {
        throw std::runtime_error("não foi possível abrir " + path);
    }

    // 2) Lê a primeira linha (o banner) com o mmio e confere se o tipo é suportado.
    // O MM_typecode guarda o tipo da matriz em 4 letras: "MCRS" significa
    // matriz, coordenada, real, simétrica.
    MM_typecode tipo;
    int erro_no_banner = mm_read_banner(arquivo, &tipo);
    if (erro_no_banner != 0) {
        fechar_e_lancar_erro(arquivo, path + ": cabeçalho Matrix Market inválido");
    }
    if (!mm_is_matrix(tipo) || !mm_is_coordinate(tipo)) {
        fechar_e_lancar_erro(arquivo, path + ": só o formato coordenado é suportado");
    }
    if (mm_is_complex(tipo) || mm_is_hermitian(tipo)) {
        fechar_e_lancar_erro(arquivo, path + ": matrizes complexas não são suportadas");
    }

    // 3) Lê a linha com as dimensões. O mmio pula os comentários antes dela.
    int linhas = 0;
    int colunas = 0;
    int entradas = 0;
    int erro_nas_dimensoes = mm_read_mtx_crd_size(arquivo, &linhas, &colunas, &entradas);
    if (erro_nas_dimensoes != 0) {
        fechar_e_lancar_erro(arquivo, path + ": erro ao ler as dimensões");
    }
    if (linhas != colunas) {
        fechar_e_lancar_erro(arquivo, path + ": a matriz não é quadrada");
    }

    // 4) Descobre como as entradas devem ser interpretadas.
    bool arquivo_sem_valores = mm_is_pattern(tipo);
    bool simetrica = mm_is_symmetric(tipo);
    bool antissimetrica = mm_is_skew(tipo);

    // 5) Preenche as informações do resultado e cria a matriz zerada.
    MtxFile resultado;
    resultado.name = base_name(path);
    resultado.stored_entries = entradas;
    resultado.A = DenseMatrix(linhas);

    if (simetrica) {
        resultado.symmetry = "symmetric";
    } else if (antissimetrica) {
        resultado.symmetry = "skew-symmetric";
    } else {
        resultado.symmetry = "general";
    }

    if (arquivo_sem_valores) {
        resultado.field = "pattern";
    } else if (mm_is_integer(tipo)) {
        resultado.field = "integer";
    } else {
        resultado.field = "real";
    }

    // 6) Lê as entradas, uma por linha do arquivo, e coloca cada uma na matriz.
    for (int k = 0; k < entradas; k++) {
        int i = 0;
        int j = 0;
        double valor = 1.0;  // valor usado nos arquivos pattern

        // O fscanf devolve quantos itens conseguiu ler da linha.
        bool leitura_ok;
        if (arquivo_sem_valores) {
            int itens_lidos = std::fscanf(arquivo, "%d %d", &i, &j);
            leitura_ok = (itens_lidos == 2);
        } else {
            int itens_lidos = std::fscanf(arquivo, "%d %d %lg", &i, &j, &valor);
            leitura_ok = (itens_lidos == 3);
        }
        if (!leitura_ok) {
            fechar_e_lancar_erro(arquivo, path + ": entrada " + std::to_string(k + 1) + " ilegível");
        }

        // O formato Matrix Market conta os índices a partir de 1, e o C++ a partir de 0.
        i = i - 1;
        j = j - 1;

        bool i_valido = (i >= 0 && i < linhas);
        bool j_valido = (j >= 0 && j < colunas);
        if (!i_valido || !j_valido) {
            fechar_e_lancar_erro(arquivo, path + ": índice fora do intervalo");
        }

        // Soma em vez de atribuir: assim, entradas repetidas se acumulam.
        resultado.A(i, j) += valor;

        // O arquivo guarda só um triângulo das matrizes simétricas, então a entrada
        // é copiada para a posição espelhada. A diagonal não tem espelho.
        bool fora_da_diagonal = (i != j);
        if (fora_da_diagonal && simetrica) {
            resultado.A(j, i) += valor;
        }
        if (fora_da_diagonal && antissimetrica) {
            resultado.A(j, i) -= valor;  // na antissimétrica, a_ji = -a_ij
        }
    }

    std::fclose(arquivo);

    // 7) Guarda também a cópia em formato de matriz.
    resultado.matriz = montar_matriz(resultado.A);

    return resultado;
}

std::vector<std::vector<double>> montar_matriz(const DenseMatrix& A) {
    int m = A.m;

    // Cria m linhas, cada uma com m zeros.
    std::vector<double> linha_zerada(m, 0.0);
    std::vector<std::vector<double>> matriz(m, linha_zerada);

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            matriz[i][j] = A(i, j);
        }
    }

    return matriz;
}
