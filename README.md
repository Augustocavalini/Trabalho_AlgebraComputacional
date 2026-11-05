* [ ]

# Trabalho_AlgebraComputacional

TP1 de Álgebra Linear Computacional (UFJF): métodos diretos e comparação com a LAPACK.

## Ambiente

Testado em Ubuntu 26.04, g++ 15.2.

### Dependências do sistema

```bash
sudo apt install build-essential liblapacke-dev libopenblas-dev libboost-graph-dev
```

`src/mmio.c` e `src/mmio.h` são do NIST (https://math.nist.gov/MatrixMarket/mmio-c.html) e já estão no repositório.

### Python (só para os gráficos)

```bash
conda create -n tp-alc python=3.12 numpy scipy matplotlib threadpoolctl
conda activate tp-alc
```

### Verificação

```bash
make check
```

Saída esperada: `dgesv` com `info=0` e `x=(1, 1, 1)`, o typecode do mmio e uma permutação RCM do Boost.

Todas as medições de tempo usam a BLAS em uma única thread:

```bash
export OPENBLAS_NUM_THREADS=1
```

## Matrizes

As matrizes não estão no repositório. Baixe-as do repositório Matrix Market do NIST para `matrices/`:

```bash
scripts/download_matrices.sh Harwell-Boeing/bcsstruc2/bcsstk14 Harwell-Boeing/oilgen/orsirr_1 Harwell-Boeing/psadmit/1138_bus
```

Páginas das matrizes:

- https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/bcsstruc2/bcsstk14.html
- https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/oilgen/orsirr_1.html
- https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/psadmit/1138_bus.html

As mesmas matrizes estão na SuiteSparse Matrix Collection (https://sparse.tamu.edu), no grupo `HB`.

## Parte 1 – Leitura e caracterização

```bash
make parte1                        # tabela em Markdown no terminal e em results/parte1.csv
conda run -n tp-alc python scripts/spy.py   # spy plots em figures/spy_<matriz>.png
```

Sem argumentos, os dois usam todas as `matrices/*.mtx`; também aceitam arquivos específicos
(`./build/parte1 matrices/bcsstk14.mtx`).

- `src/matrix.{hpp,cpp}`: leitura do `.mtx` com o mmio e montagem da forma densa em coluna-maior.
  Matrizes `symmetric` (e `skew-symmetric`) têm os coeficientes fora da diagonal espelhados;
  arquivos `pattern` recebem valor 1; entradas repetidas são somadas.
  A matriz fica em dois formatos: `A` (vetor coluna-maior, usado pela LAPACK) e `matriz[i][j]`
  (cópia em formato de matriz, feita ao final da leitura).
- `src/parte1.cpp`: nnz e densidade contados na forma densa (já espelhada), simetria verificada
  numericamente, bandas `b_I = max(i-j)` e `b_S = max(j-i)` sobre os não nulos, e κ₁(A) estimado
  com `dlange` + `dpocon` se `dpotrf` tiver sucesso (matriz SPD) ou `dlange` + `dgetrf` + `dgecon`
  caso contrário.

## Parte 2 – Validação dos métodos

```bash
make testar                                        # usa matrices/teste_10x10.mtx
./build/testar_metodos matrices/bcsstk14.mtx       # usa outra matriz
```

`src/parte2.cpp` resolve o sistema com solução exata `x = (1, ..., 1)` por LU densa, Cholesky denso,
LU em banda, Cholesky em banda e Thomas, e compara cada solução com a da rotina correspondente da LAPACK.
Um método que não se aplica à matriz aparece com o motivo. Os métodos estão em `src/parte2.{hpp,cpp}`.

`matrices/teste_10x10.mtx` é uma matriz tridiagonal SPD, escolhida para funcionar em todos os métodos.
Ela não fica no repositório, como as demais matrizes (ver `.gitignore`). Para recriá-la, grave em
`matrices/teste_10x10.mtx`:

```text
%%MatrixMarket matrix coordinate real symmetric
10 10 19
1 1 4
2 1 -1
2 2 5
3 2 2
3 3 6
4 3 -1
4 4 4
5 4 -2
5 5 5
6 5 1
6 6 6
7 6 -1
7 7 4
8 7 2
8 8 5
9 8 -2
9 9 6
10 9 1
10 10 4
```

## Parte 3 – Reordenamento RCM

```bash
make parte3 MATRIZ=matrices/bcsstk14.mtx    # bandas e armazenamento: results/parte3.csv; tempos: results/parte3_tempos.csv
make spy_rcm MATRIZ=matrices/bcsstk14.mtx   # padrões lado a lado: figures/rcm_bcsstk14.png
```

`src/rcm.{hpp,cpp}` calcula a permutação com `boost::cuthill_mckee_ordering`, a partir do grafo de adjacência de A
(uma aresta sem direção para cada `a_ij != 0`, `i != j`), e forma `PAPᵀ`. Em `perm`, `perm[k]` é a linha da matriz
original que passa a ocupar a posição `k`. O programa `src/parte3.cpp` confere que `(PAPᵀ)(Px) = Pb` e grava a
permutação em `results/rcm_<matriz>.txt`, que é usada pelo `scripts/spy_rcm.py`.

O programa também resolve o sistema reordenado nas formas densa e em banda, com Cholesky se a matriz for SPD e com
LU nos outros casos, e reporta a mediana de 5 execuções. O ganho previsto vem das contagens de operações:
`(1/3)m³ ÷ m·b²` no Cholesky e `(2/3)m³ ÷ 2m·b_I·b_S` na LU. Os erros são medidos no sistema original, depois de
recuperar `x` a partir de `Px`.

## Parte 4 – Comparação com a LAPACK

```bash
make parte4 MATRIZES="matrices/bcsstk14.mtx matrices/orsirr_1.mtx"   # results/parte4_matrizes.csv
make parte4_escala                                                   # results/parte4_escala.csv (alguns minutos)
make graficos_escala                                                 # figures/parte4_escala.png e results/parte4_inclinacoes.csv
```

Sem `MATRIZES=`, o `make parte4` usa todas as `matrices/*.mtx`. Cada método que se aplica à matriz é medido na nossa
implementação e na rotina da LAPACK, com a mediana de 5 execuções, alternando as duas. `src/parte4.cpp` desliga a
verificação de NaN do LAPACKE (`LAPACKE_set_nancheck(0)`), que não faz parte das rotinas da LAPACK. O limite do erro
direto usa κ∞(A), estimado com `dgecon`, porque os erros são medidos na norma infinito.
