#include "rcm.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>

std::vector<int> calcular_rcm(const DenseMatrix& A) {
    int m = A.m;

    // Tipo do grafo: não direcionado, com as propriedades que o algoritmo exige
    // em cada vértice (uma cor, para marcar visitados, e o grau).
    typedef boost::property<boost::vertex_degree_t, int> PropriedadeGrau;
    typedef boost::property<boost::vertex_color_t, boost::default_color_type, PropriedadeGrau> PropriedadesVertice;
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, PropriedadesVertice> Grafo;

    Grafo grafo(m);

    // Cada par (i, j) é visitado uma só vez, com i > j. Se a_ij e a_ji fossem tratados
    // separadamente, uma matriz simétrica teria arestas duplicadas e graus errados.
    for (int j = 0; j < m; j++) {
        for (int i = j + 1; i < m; i++) {
            if (A(i, j) != 0.0 || A(j, i) != 0.0) {
                boost::add_edge(i, j, grafo);
            }
        }
    }

    std::vector<Grafo::vertex_descriptor> ordem(m);

    // Cuthill-McKee escreve a ordem no iterador recebido. Passando rbegin(), o vetor é
    // preenchido de trás para a frente, o que dá a ordem reversa (RCM).
    // Se o grafo tiver partes desconectadas, a Boost percorre cada uma delas.
    boost::cuthill_mckee_ordering(grafo,
                                  ordem.rbegin(),
                                  boost::get(boost::vertex_color, grafo),
                                  boost::make_degree_map(grafo));

    std::vector<int> perm(m);
    for (int k = 0; k < m; k++) {
        perm[k] = static_cast<int>(ordem[k]);
    }
    return perm;
}

DenseMatrix permutar_matriz(const DenseMatrix& A, const std::vector<int>& perm) {
    int m = A.m;
    DenseMatrix permutada(m);

    for (int j = 0; j < m; j++) {
        for (int i = 0; i < m; i++) {
            permutada(i, j) = A(perm[i], perm[j]);
        }
    }
    return permutada;
}

std::vector<double> permutar_vetor(const std::vector<double>& v, const std::vector<int>& perm) {
    int m = static_cast<int>(v.size());
    std::vector<double> permutado(m);

    for (int k = 0; k < m; k++) {
        permutado[k] = v[perm[k]];
    }
    return permutado;
}

std::vector<double> despermutar_vetor(const std::vector<double>& v, const std::vector<int>& perm) {
    int m = static_cast<int>(v.size());
    std::vector<double> original(m);

    for (int k = 0; k < m; k++) {
        original[perm[k]] = v[k];
    }
    return original;
}
