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
pip install reportlab    # só para gerar o relatório em PDF
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

As matrizes não estão no repositório. As três usadas no relatório vêm do repositório Matrix Market do NIST:

- `bcsstk15`: https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/bcsstruc2/bcsstk15.html
- `steam2`: https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/steam/steam2.html
- `1138_bus`: https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/psadmit/1138_bus.html

Baixe o `.mtx.gz` pela página de cada matriz e descompacte-o em `matrices/`, por exemplo:

```bash
gunzip -c ~/Downloads/steam2.mtx.gz > matrices/steam2.mtx
```

As mesmas matrizes estão na SuiteSparse Matrix Collection (https://sparse.tamu.edu), no grupo `HB`.

## Reproduzir o relatório

Com `matrices/` contendo apenas as três matrizes acima (nos comandos sem argumentos, os programas usam todas as `matrices/*.mtx`):

| Item do relatório | Comando | Saída |
|---|---|---|
| Tabela 1 | `OPENBLAS_NUM_THREADS=1 ./build/testar_metodos teste_10x10.mtx > results/parte2_validacao.txt` (a matriz é a da seção Parte 2) | `results/parte2_validacao.txt` |
| Tabela 2 | `make parte1` | `results/parte1.csv` |
| Figura 1 | `make spy` | `figures/spy_<matriz>.png` |
| Tabelas 3 a 6 | `make parte3 MATRIZ="matrices/steam2.mtx matrices/1138_bus.mtx matrices/bcsstk15.mtx"` | `results/parte3.csv`, `results/parte3_tempos.csv` |
| Figura 2 | `make spy_rcm MATRIZ="matrices/steam2.mtx matrices/1138_bus.mtx matrices/bcsstk15.mtx"` | `figures/rcm_<matriz>.png` |
| Tabela 7 | `make parte4` | `results/parte4_matrizes.csv` |
| Tabelas 8 e 9 | `make parte4_escala` (alguns minutos) | `results/parte4_escala.csv` |
| Figura 3 e Tabela 10 | `make graficos_escala` | `figures/parte4_escala.png`, `results/parte4_inclinacoes.csv` |
| Relatório em PDF | `make relatorio` | `relatorio/relatorio.pdf` |

## Parte 1 – Leitura e caracterização

```bash
make parte1                        # tabela em Markdown no terminal e em results/parte1.csv
conda run -n tp-alc python scripts/spy.py   # spy plots em figures/spy_<matriz>.png
```

Sem argumentos, os dois usam todas as `matrices/*.mtx`; também aceitam arquivos específicos
(`./build/parte1 matrices/bcsstk15.mtx`).

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
./build/testar_metodos matrices/bcsstk15.mtx       # usa outra matriz
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
make parte3 MATRIZ=matrices/bcsstk15.mtx    # bandas e armazenamento: results/parte3.csv; tempos: results/parte3_tempos.csv
make spy_rcm MATRIZ=matrices/bcsstk15.mtx   # padrões lado a lado: figures/rcm_bcsstk15.png
```

`src/rcm.{hpp,cpp}` calcula a permutação com `boost::cuthill_mckee_ordering`, a partir do grafo de adjacência de A
(uma aresta sem direção para cada `a_ij != 0`, `i != j`), e forma `PAPᵀ`. Em `perm`, `perm[k]` é a linha da matriz
original que passa a ocupar a posição `k`. O programa `src/parte3.cpp` confere que `(PAPᵀ)(Px) = Pb` e grava a
permutação em `results/rcm_<matriz>.txt`, que é usada pelo `scripts/spy_rcm.py`.

O programa também resolve o sistema com LU (todas as matrizes) e com Cholesky (matrizes SPD), nas formas densa e em
banda, sem RCM (matriz original) e com RCM (`PAPᵀ`), e reporta a mediana de 5 execuções. O ganho previsto vem das
contagens de operações, pela fórmula do enunciado (`(2/3)m³ ÷ 2m·b_I·b_S` na LU, `(1/3)m³ ÷ m·b²` no Cholesky) e pela
contagem exata das operações dos laços. Os erros são medidos no sistema original, depois de recuperar `x` a partir de `Px`.

## Parte 4 – Comparação com a LAPACK

```bash
make parte4 MATRIZES="matrices/bcsstk15.mtx matrices/steam2.mtx"   # results/parte4_matrizes.csv
make parte4_escala                                                   # results/parte4_escala.csv (alguns minutos)
make graficos_escala                                                 # figures/parte4_escala.png e results/parte4_inclinacoes.csv
```

Sem `MATRIZES=`, o `make parte4` usa todas as `matrices/*.mtx`. Cada método que se aplica à matriz é medido na nossa
implementação e na rotina da LAPACK, com a mediana de 5 execuções, alternando as duas. `src/parte4.cpp` desliga a
verificação de NaN do LAPACKE (`LAPACKE_set_nancheck(0)`), que não faz parte das rotinas da LAPACK. O limite do erro
direto usa κ∞(A), estimado com `dgecon`, porque os erros são medidos na norma infinito.
