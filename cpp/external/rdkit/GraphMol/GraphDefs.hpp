#ifndef RDKIT_GRAPH_DEFS_H
#define RDKIT_GRAPH_DEFS_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include <boost/graph/graph_traits.hpp>

#include "csr_undirected_extension.hpp"

// clang-format off
namespace RDKit {
class Atom;
class Bond;

//! These are the BGL types used to store the topology:
using MolGraph    = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, Atom*, Bond*>;
using CSRMolGraph = boost::compressed_sparse_row_graph<boost::undirectedS, Atom*, Bond*>;


enum class GraphType { MolGraph, CSRMolGraph };

template <typename GraphType>
struct MolGraphTraits {};

template <>
struct MolGraphTraits<MolGraph> {
    using vertex_descriptor   = boost::graph_traits<MolGraph>::vertex_descriptor;
    using edge_descriptor     = boost::graph_traits<MolGraph>::edge_descriptor;
    using vertex_iterator     = boost::graph_traits<MolGraph>::vertex_iterator;
    using edge_iterator       = boost::graph_traits<MolGraph>::edge_iterator;
    using out_edge_iterator   = boost::graph_traits<MolGraph>::out_edge_iterator;
    using adjacency_iterator  = boost::graph_traits<MolGraph>::adjacency_iterator;

    static Atom* getAtom(const MolGraph& graph, vertex_descriptor vertex) {
        return graph[vertex];
    }

    static Bond* getBond(const MolGraph& graph, edge_descriptor edge) {
        return graph[edge];
    }
};

template <>
struct MolGraphTraits<CSRMolGraph> {
    using vertex_descriptor   = boost::graph_traits<CSRMolGraph>::vertex_descriptor;
    using edge_descriptor     = boost::graph_traits<CSRMolGraph>::edge_descriptor;
    using vertex_iterator     = boost::graph_traits<CSRMolGraph>::vertex_iterator;
    using edge_iterator       = boost::graph_traits<CSRMolGraph>::edge_iterator;
    using out_edge_iterator   = boost::graph_traits<CSRMolGraph>::out_edge_iterator;
    using adjacency_iterator  = boost::graph_traits<CSRMolGraph>::adjacency_iterator;

    static Atom* getAtom(const CSRMolGraph& graph, vertex_descriptor vertex) {
        return graph[vertex];
    }

    static Bond* getBond(const CSRMolGraph& graph, edge_descriptor edge) {
        return graph[edge];
    }
};
} // namespace RDKit

#endif // RDKIT_GRAPH_DEFS_H

