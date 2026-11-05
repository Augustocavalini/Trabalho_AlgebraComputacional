#pragma once

#include <vector>

#include <lapacke.h>

#include "matrix.hpp"

// Fatoração PA = LU com o dgetrf da LAPACK. Usada como referência de comparação.
// Sobrescreve A com os fatores. O vetor pivos deve ter m posições.
void fatoracaoLULapacke(DenseMatrix& A, std::vector<lapack_int>& pivos);

// Nossa fatoração PA = LU com pivoteamento parcial.
//
// L e U sobrescrevem A: U fica no triângulo superior, com a diagonal, e L abaixo da
// diagonal. A diagonal de L vale 1 e não é guardada. O vetor pivos deve ter m posições;
// pivos[k] é a linha trocada com a linha k na etapa k, contada a partir de 1, como na LAPACK.
// Se encontrar um pivô nulo, lança std::runtime_error.
void fatoracaoLu(DenseMatrix& A, std::vector<lapack_int>& pivos);

// Repete em b as trocas de linha da fatoração, transformando b em Pb.
void aplicar_pivos(const std::vector<lapack_int>& pivos, std::vector<double>& b);

// Substituição direta: resolve Ly = b. Recebe b em y e devolve a solução em y.
void substituicao_direta(const DenseMatrix& LU, std::vector<double>& y);

// Substituição inversa: resolve Ux = y. Recebe y em x e devolve a solução em x.
void substituicao_inversa(const DenseMatrix& LU, std::vector<double>& x);

// Resolve Ax = b a partir de uma fatoração já calculada: aplica os pivôs e faz as
// substituições direta e inversa. A fatoração pode vir de fatoracaoLu ou de fatoracaoLULapacke.
std::vector<double> resolver_sistema_lu(const DenseMatrix& LU, const std::vector<lapack_int>& pivos, const std::vector<double>& b);
