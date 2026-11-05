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

// Fatoração A = LL^T com o dpotrf da LAPACK. Usada como referência de comparação.
// Sobrescreve o triângulo inferior de A com L.
void fatoracaoCholeskyLapacke(DenseMatrix& A);

// Nossa fatoração de Cholesky, A = LL^T, para matrizes simétricas positivas definidas.
//
// L sobrescreve o triângulo inferior de A, com a diagonal. Só esse triângulo é lido, e o
// triângulo superior não é alterado. Se encontrar um pivô não positivo, a matriz não é SPD
// e a função lança std::runtime_error, informando a etapa, contada a partir de 1.
void fatoracaoCholesky(DenseMatrix& A);

// Substituição direta: resolve Ly = b. Recebe b em y e devolve a solução em y.
void substituicao_direta_cholesky(const DenseMatrix& L, std::vector<double>& y);

// Substituição inversa: resolve L^T x = y. Recebe y em x e devolve a solução em x.
void substituicao_inversa_cholesky(const DenseMatrix& L, std::vector<double>& x);

// Resolve Ax = b a partir de uma fatoração já calculada, com as substituições direta e
// inversa. A fatoração pode vir de fatoracaoCholesky ou de fatoracaoCholeskyLapacke.
std::vector<double> resolver_sistema_cholesky(const DenseMatrix& L, const std::vector<double>& b);

// Nossa LU em banda, sem pivoteamento: A = LU.
//
// L e U sobrescrevem A dentro da própria banda: L fica nas bI diagonais inferiores, com a
// diagonal valendo 1 sem ser guardada, e U nas demais. Se encontrar um pivô nulo, lança
// std::runtime_error. Sem pivoteamento, pivôs pequenos podem gerar erros grandes.
void fatoracaoLuBanda(BandMatrix& A);

// Substituição direta, restrita à banda: resolve Ly = b. Recebe b em y e devolve a solução em y.
void substituicao_direta_banda(const BandMatrix& LU, std::vector<double>& y);

// Substituição inversa, restrita à banda: resolve Ux = y. Recebe y em x e devolve a solução em x.
void substituicao_inversa_banda(const BandMatrix& LU, std::vector<double>& x);

// Resolve Ax = b a partir de uma fatoração de fatoracaoLuBanda.
std::vector<double> resolver_sistema_lu_banda(const BandMatrix& LU, const std::vector<double>& b);

// Fatoração em banda com o dgbtrf da LAPACK, que faz pivoteamento. Usada como referência.
// A deve vir de montar_banda(A_densa, bI, bI + bS), que reserva as posições extras exigidas
// (2bI + bS + 1 linhas). O vetor pivos deve ter m posições. A solução usa LAPACKE_dgbtrs.
void fatoracaoLuBandaLapacke(BandMatrix& A, std::vector<lapack_int>& pivos);

// Nosso Cholesky em banda, A = LL^T, para matrizes SPD.
//
// A deve vir de montar_banda(A_densa, b, 0), que guarda só o triângulo inferior
// (b + 1 linhas). L sobrescreve A. Se encontrar um pivô não positivo, a matriz não é SPD
// e a função lança std::runtime_error, informando a etapa, contada a partir de 1.
void fatoracaoCholeskyBanda(BandMatrix& A);

// Substituição direta, restrita à banda: resolve Ly = b. Recebe b em y e devolve a solução em y.
void substituicao_direta_cholesky_banda(const BandMatrix& L, std::vector<double>& y);

// Substituição inversa, restrita à banda: resolve L^T x = y. Recebe y em x e devolve a solução em x.
void substituicao_inversa_cholesky_banda(const BandMatrix& L, std::vector<double>& x);

// Resolve Ax = b a partir de uma fatoração de fatoracaoCholeskyBanda ou fatoracaoCholeskyBandaLapacke.
std::vector<double> resolver_sistema_cholesky_banda(const BandMatrix& L, const std::vector<double>& b);

// Cholesky em banda com o dpbtrf da LAPACK. Usado como referência. A deve vir de
// montar_banda(A_densa, b, 0), como em fatoracaoCholeskyBanda.
void fatoracaoCholeskyBandaLapacke(BandMatrix& A);

// Algoritmo de Thomas: LU sem pivoteamento para matrizes tridiagonais, usando apenas três vetores.
//
//     | d0 s0  0  0 |      diagonal = (d0, d1, d2, d3)   m posições
//     | i0 d1 s1  0 |      inferior = (i0, i1, i2)       m - 1 posições
//     |  0 i1 d2 s2 |      superior = (s0, s1, s2)       m - 1 posições
//     |  0  0 i2 d3 |
//
// Sobrescreve inferior com os multiplicadores de L e diagonal com a diagonal de U.
// O vetor superior não muda, porque já é a parte de U acima da diagonal.
// Se encontrar um pivô nulo, lança std::runtime_error.
void fatoracaoThomas(std::vector<double>& inferior, std::vector<double>& diagonal, const std::vector<double>& superior);

// Substituição direta: resolve Ly = b. Recebe b em y e devolve a solução em y.
void substituicao_direta_thomas(const std::vector<double>& inferior, std::vector<double>& y);

// Substituição inversa: resolve Ux = y. Recebe y em x e devolve a solução em x.
void substituicao_inversa_thomas(const std::vector<double>& diagonal, const std::vector<double>& superior, std::vector<double>& x);

// Resolve Ax = b a partir dos três vetores já fatorados por fatoracaoThomas.
std::vector<double> resolver_sistema_thomas(const std::vector<double>& inferior, const std::vector<double>& diagonal, const std::vector<double>& superior, const std::vector<double>& b);

// Resolve Ax = b com o dgtsv da LAPACK, que faz fatoração e solução juntas, com pivoteamento.
// Usada como referência. Sobrescreve os três vetores.
std::vector<double> resolver_tridiagonal_lapacke(std::vector<double>& inferior, std::vector<double>& diagonal, std::vector<double>& superior, const std::vector<double>& b);
