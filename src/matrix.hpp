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

// Matriz em banda: guarda só as diagonais entre a banda inferior bI e a banda superior bS.
//
// O arranjo é o mesmo da LAPACK. Cada coluna da matriz vira uma coluna de um arranjo com
// bI + bS + 1 linhas, e cada diagonal da matriz vira uma linha desse arranjo. Exemplo com
// m = 4, bI = 1 e bS = 1 (as posições marcadas com * nunca são usadas):
//
//     | a00 a01  0   0  |                 |  *  a01 a12 a23 |   diagonal superior
//     | a10 a11 a12  0  |   é guardada    | a00 a11 a22 a33 |   diagonal principal
//     |  0  a21 a22 a23 |   como          | a10 a21 a32  *  |   diagonal inferior
//     |  0   0  a32 a33 |
//
// As colunas desse arranjo ficam uma depois da outra no vetor a, que tem
// (bI + bS + 1) * m posições. Só se pode acessar A(i, j) dentro da banda, isto é,
// com i - j <= bI e j - i <= bS.
struct BandMatrix {
    int m;                  // número de linhas (e de colunas) da matriz
    int bI;                 // largura de banda inferior
    int bS;                 // largura de banda superior
    int linhas;             // número de linhas do arranjo: bI + bS + 1
    std::vector<double> a;  // os coeficientes da banda, coluna por coluna

    // Cria uma matriz em banda com todos os coeficientes zerados.
    BandMatrix(int dimensao, int banda_inferior, int banda_superior) {
        m = dimensao;
        bI = banda_inferior;
        bS = banda_superior;
        linhas = bI + bS + 1;

        int total_de_posicoes = linhas * m;
        a.assign(total_de_posicoes, 0.0);
    }

    // Converte (linha i, coluna j) da matriz na posição correspondente do vetor a.
    //
    // Antes da coluna j existem j colunas do arranjo. Dentro da coluna j, a diagonal
    // principal fica na linha bS do arranjo, e cada passo para baixo na matriz
    // é também um passo para baixo no arranjo.
    int posicao(int i, int j) const {
        int coeficientes_antes_da_coluna = j * linhas;
        int linha_no_arranjo = bS + i - j;
        int posicao_no_vetor = coeficientes_antes_da_coluna + linha_no_arranjo;
        return posicao_no_vetor;
    }

    // Permite escrever na banda com a sintaxe A(i, j), por exemplo: A(2, 1) = 5.0;
    double& operator()(int i, int j) {
        int p = posicao(i, j);
        return a[p];
    }

    // Permite ler a banda com a mesma sintaxe, por exemplo: double v = A(2, 1);
    double operator()(int i, int j) const {
        int p = posicao(i, j);
        return a[p];
    }
};

// Calcula as larguras de banda de A: bI = max(i - j) e bS = max(j - i) sobre os não nulos.
void calcular_bandas(const DenseMatrix& A, int& bI, int& bS);

// Copia para o armazenamento em banda os coeficientes de A dentro das larguras pedidas.
// Coeficientes fora dessas larguras não são copiados. Isso é usado de propósito no Cholesky
// em banda, que guarda apenas o triângulo inferior: montar_banda(A, b, 0).
BandMatrix montar_banda(const DenseMatrix& A, int bI, int bS);

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
