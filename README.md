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

As matrizes não estão no repositório. Baixe-as da SuiteSparse Matrix Collection para `matrices/`:

```bash
scripts/download_matrices.sh HB/bcsstk14 HB/orsirr_1 HB/1138_bus
```

Se o download automático falhar, baixe o arquivo `.tar.gz` (formato Matrix Market) pela página de cada matriz e extraia o `.mtx` em `matrices/`:

- https://sparse.tamu.edu/HB/bcsstk14
- https://sparse.tamu.edu/HB/orsirr_1
- https://sparse.tamu.edu/HB/1138_bus

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

`src/testar_metodos.cpp` resolve o sistema com solução exata `x = (1, ..., 1)` por LU densa, Cholesky denso,
LU em banda, Cholesky em banda e Thomas, e compara cada solução com a da rotina correspondente da LAPACK.
Um método que não se aplica à matriz aparece com o motivo. Os métodos estão em `src/metodos.{hpp,cpp}`.

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
