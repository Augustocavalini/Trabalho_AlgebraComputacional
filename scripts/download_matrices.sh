#!/usr/bin/env bash
# Baixa matrizes da SuiteSparse Matrix Collection (https://sparse.tamu.edu) para matrices/.
# Uso: scripts/download_matrices.sh [Grupo/Nome ...]
# Exemplo: scripts/download_matrices.sh HB/bcsstk14 HB/orsirr_1

# Interrompe o script no primeiro erro.
set -euo pipefail

DESTINO="$(dirname "$0")/../matrices"
mkdir -p "$DESTINO"

# Sem argumentos, baixa as matrizes sugeridas.
if [ "$#" -eq 0 ]; then
    set -- HB/bcsstk14 HB/orsirr_1 HB/1138_bus
fi

for matriz in "$@"; do
    # Remove o grupo do nome: "HB/bcsstk14" vira "bcsstk14".
    nome="${matriz#*/}"

    if [ -f "$DESTINO/$nome.mtx" ]; then
        echo "já existe: $nome.mtx"
        continue
    fi

    echo "baixando $matriz"

    # O curl baixa o .tar.gz e o tar extrai só o .mtx. O --strip-components=1
    # descarta a pasta interna, de modo que o arquivo fica direto em matrices/.
    curl -fsSL "https://sparse.tamu.edu/MM/$matriz.tar.gz" \
        | tar -xz -C "$DESTINO" --strip-components=1 "$nome/$nome.mtx"
done
