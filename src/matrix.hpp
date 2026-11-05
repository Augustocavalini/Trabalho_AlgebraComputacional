#pragma once

#include <string>
#include <vector>

// Matriz densa quadrada guardada em ordem coluna-maior, como a LAPACK espera.
//
// Os coeficientes ficam em um único vetor, com as colunas uma depois da outra.
// Exemplo com m = 3:
//
//     | a00 a01 a02 |
//     | a10 a11 a12 |   vira   a = [a00, a10, a20,  a01, a11, a21,  a02, a12, a22]
//     | a20 a21 a22 |              (coluna 0)       (coluna 1)      (coluna 2)
//
// Como m <= 5000 neste trabalho, o vetor tem no máximo 25 milhões de posições,
// e por isso um int é suficiente para guardar qualquer posição.
struct DenseMatrix {
    int m;                  // número de linhas (e de colunas)
    std::vector<double> a;  // os m*m coeficientes

    // Construtor sem argumentos: cria uma matriz vazia, de dimensão 0.
    // Ele é chamado quando escrevemos apenas "DenseMatrix A;".
    DenseMatrix() {
        m = 0;
    }

    // Construtor com a dimensão: cria uma matriz m por m com todos os coeficientes zerados.
    // Ele é chamado quando escrevemos "DenseMatrix A(100);".
    DenseMatrix(int dimensao) {
        m = dimensao;

        int total_de_posicoes = m * m;
        a.assign(total_de_posicoes, 0.0);  // o vetor passa a ter essa quantidade de zeros
    }

    // Converte (linha i, coluna j) na posição correspondente do vetor a.
    //
    // Antes da coluna j existem j colunas completas, cada uma com m coeficientes.
    // Depois, andamos mais i posições dentro da coluna j.
    int posicao(int i, int j) const {
        int coeficientes_antes_da_coluna = j * m;
        int posicao_no_vetor = coeficientes_antes_da_coluna + i;
        return posicao_no_vetor;
    }

    // Permite escrever na matriz com a sintaxe A(i, j), por exemplo: A(2, 1) = 5.0;
    // O "&" no tipo de retorno faz a função devolver o próprio coeficiente guardado
    // no vetor, e não uma cópia. É por isso que a atribuição altera a matriz.
    double& operator()(int i, int j) {
        int p = posicao(i, j);
        return a[p];
    }

    // Permite ler a matriz com a mesma sintaxe, por exemplo: double v = A(2, 1);
    // Esta versão é usada quando a matriz é "const", isto é, quando não pode ser alterada.
    // Ela devolve uma cópia do valor, sem "&", e o "const" no fim promete não mudar a matriz.
    double operator()(int i, int j) const {
        int p = posicao(i, j);
        return a[p];
    }
};

// Conteúdo de um arquivo Matrix Market já montado na forma densa.
struct MtxFile {
    std::string name;      // nome do arquivo, sem diretório e sem extensão
    std::string symmetry;  // cabeçalho: general, symmetric ou skew-symmetric
    std::string field;     // cabeçalho: real, integer ou pattern
    int stored_entries;    // entradas guardadas no arquivo (só um triângulo, se simétrica)
    DenseMatrix A;         // a matriz completa, já espelhada, em um único vetor

    // Os mesmos coeficientes em formato de matriz: matriz[i][j] é a linha i, coluna j.
    // É uma cópia feita ao final da leitura: alterar A depois não altera esta matriz.
    std::vector<std::vector<double>> matriz;
};

// Copia os coeficientes do vetor coluna-maior para uma matriz com m linhas e m colunas.
std::vector<std::vector<double>> montar_matriz(const DenseMatrix& A);

// Lê um arquivo .mtx em formato coordenado e monta a forma densa.
//
// - Em matrizes symmetric e skew-symmetric, espelha os coeficientes fora da diagonal.
// - Em arquivos pattern, que não têm valores, cada entrada recebe o valor 1.
// - Entradas repetidas são somadas.
// - Em caso de erro, lança uma exceção std::runtime_error com a mensagem do problema.
MtxFile read_mtx(const std::string& path);
