#pragma once

#include <vector>

#include "matrix.hpp"

// Reordenamento Cuthill-McKee reverso (RCM).
//
// A permutação é guardada como um vetor perm: perm[k] é a linha ou coluna da matriz
// original que passa a ocupar a posição k. Exemplo com perm = (2, 0, 1): a nova linha 0
// é a antiga linha 2, a nova linha 1 é a antiga linha 0 e a nova linha 2 é a antiga linha 1.
// Isso equivale a uma matriz de permutação P com 1 na posição (k, perm[k]).

// Calcula a permutação RCM com a Boost, a partir do grafo de adjacência de A: um vértice por
// linha e uma aresta para cada a_ij != 0, com i != j. A aresta não tem direção, então em
// uma matriz não simétrica basta que a_ij ou a_ji seja não nulo.
std::vector<int> calcular_rcm(const DenseMatrix& A);

// Forma a matriz reordenada PAP^T: o coeficiente (i, j) dela é A(perm[i], perm[j]).
DenseMatrix permutar_matriz(const DenseMatrix& A, const std::vector<int>& perm);

// Calcula Pv: a posição k do resultado é v[perm[k]].
std::vector<double> permutar_vetor(const std::vector<double>& v, const std::vector<int>& perm);

// Desfaz a permutação, calculando P^T v. Serve para recuperar x a partir de Px.
std::vector<double> despermutar_vetor(const std::vector<double>& v, const std::vector<int>& perm);
