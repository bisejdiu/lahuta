#ifndef BOOST_GRAPH_CSR_UNDIRECTED_EXTENSION_HPP
#define BOOST_GRAPH_CSR_UNDIRECTED_EXTENSION_HPP

#include <boost/graph/compressed_sparse_row_graph.hpp>

namespace boost
{

//
// undirectedS specialization: inherits the directed implementation
// and merely flips the trait tag. All edges are duplicated (u,v) & (v,u).
// This is not the most efficient, not even complete, implementation,
// but it is the simplest and works fine.
// Tried to match the styling guide. - Besian, August 2025
//
template < typename VertexProperty, typename EdgeProperty, typename GraphProperty,
    typename Vertex, typename EdgeIndex >
class compressed_sparse_row_graph< undirectedS, VertexProperty, EdgeProperty,
    GraphProperty, Vertex, EdgeIndex >
: public compressed_sparse_row_graph< directedS, VertexProperty, EdgeProperty,
      GraphProperty, Vertex, EdgeIndex >
{
public:
    typedef compressed_sparse_row_graph< directedS, VertexProperty, EdgeProperty,
        GraphProperty, Vertex, EdgeIndex > base_type;

    using base_type::base_type;

    // trait overrides
    typedef undirected_tag directed_category;
    typedef typename base_type::out_edge_iterator in_edge_iterator;

    // public typedefs
    typedef typename base_type::vertex_descriptor vertex_descriptor;
    typedef typename base_type::edge_descriptor edge_descriptor;
    typedef typename base_type::out_edge_iterator out_edge_iterator;
    typedef typename base_type::degree_size_type degree_size_type;

    // in-edge wrappers
    friend inline std::pair< in_edge_iterator, in_edge_iterator > in_edges(
        vertex_descriptor v, const compressed_sparse_row_graph& g)
    {
        return out_edges(v, static_cast< const base_type& >(g));
    }

    friend inline degree_size_type in_degree(
        vertex_descriptor v, const compressed_sparse_row_graph& g)
    {
        return out_degree(v, static_cast< const base_type& >(g));
    }
};

} // end namespace boost

#endif // BOOST_GRAPH_CSR_UNDIRECTED_EXTENSION_HPP
