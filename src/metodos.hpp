#pragma once

#include <vector>

#include <lapacke.h>

#include "matrix.hpp"

// Resolve Ax = b com o dgetrs da LAPACK, a partir de uma fatoração LU já calculada.
// A fatoração pode vir de fatoracaoLULapacke ou de fatoracaoLu.
std::vector<double> resolver_sistema_lu(const DenseMatrix& LU, const std::vector<lapack_int>& pivos, const std::vector<double>& b);

// Fatoração PA = LU com o dgetrf da LAPACK. O vetor pivos deve ter m posições.
DenseMatrix fatoracaoLULapacke(const DenseMatrix& A, std::vector<lapack_int>& pivos);

// Nossa fatoração PA = LU com pivoteamento parcial. O vetor pivos deve ter m posições
// e recebe as linhas contadas a partir de 1, como na LAPACK.
DenseMatrix fatoracaoLu(const DenseMatrix& A, std::vector<lapack_int>& pivos);
