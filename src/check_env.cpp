// Teste rápido do ambiente: confere se LAPACKE, mmio e Boost estão funcionando.
#include <cstdio>
#include <vector>

#include <lapacke.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>

extern "C" {
#include "mmio.h"
}

// LAPACKE: resolve um sistema 3x3 com dgesv. A solução esperada é x = (1, 1, 1).
// Devolve o info do dgesv, que é 0 quando tudo dá certo.
static int testar_lapacke() {
    // Matriz em ordem coluna-maior: cada grupo de três números é uma coluna.
    double A[9] = {
        4, 1, 0,
        1, 4, 1,
        0, 1, 4,
    };
    double b[3] = {5, 6, 5};  // na saída, b contém a solução x
    lapack_int pivos[3];

    lapack_int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, 3, 1, A, 3, pivos, b, 3);

    std::printf("LAPACKE dgesv: info=%d x=(%g, %g, %g)\n", static_cast<int>(info), b[0], b[1], b[2]);
    return static_cast<int>(info);
}

// mmio: monta um código de tipo e pede a descrição dele em texto.
static void testar_mmio() {
    MM_typecode tipo;
    mm_initialize_typecode(&tipo);
    mm_set_matrix(&tipo);
    mm_set_coordinate(&tipo);
    mm_set_real(&tipo);
    mm_set_symmetric(&tipo);

    std::printf("mmio: typecode = %s\n", mm_typecode_to_str(tipo));
}

// Boost: calcula a ordem RCM de um caminho 0-3-1-4-2 (cinco vértices em fila).
static void testar_boost() {
    // Tipo do grafo: não direcionado, com as propriedades que o algoritmo exige
    // em cada vértice (uma cor, para marcar visitados, e o grau).
    typedef boost::property<boost::vertex_degree_t, int> PropriedadeGrau;
    typedef boost::property<boost::vertex_color_t, boost::default_color_type, PropriedadeGrau> PropriedadesVertice;
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, PropriedadesVertice> Grafo;

    Grafo grafo(5);
    boost::add_edge(0, 3, grafo);
    boost::add_edge(3, 1, grafo);
    boost::add_edge(1, 4, grafo);
    boost::add_edge(4, 2, grafo);

    std::vector<Grafo::vertex_descriptor> ordem(boost::num_vertices(grafo));

    // Cuthill-McKee escreve a ordem no iterador recebido. Passando rbegin(), o vetor é
    // preenchido de trás para a frente, o que dá a ordem reversa (RCM).
    boost::cuthill_mckee_ordering(grafo,
                                  ordem.rbegin(),
                                  boost::get(boost::vertex_color, grafo),
                                  boost::make_degree_map(grafo));

    std::printf("Boost RCM:");
    for (size_t k = 0; k < ordem.size(); k++) {
        std::printf(" %d", static_cast<int>(ordem[k]));
    }
    std::printf("\n");
}

int main() {
    int info = testar_lapacke();
    testar_mmio();
    testar_boost();

    if (info != 0) {
        return 1;
    }
    return 0;
}
