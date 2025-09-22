#ifndef RDKIT_GRAPH_ITERATORS_H
#define RDKIT_GRAPH_ITERATORS_H

#include <cassert>
#include <iterator>
#include <type_traits>
#include <utility>
#include <variant>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include <boost/graph/graph_traits.hpp>

#include "GraphDefs.hpp"

// clang-format off
namespace RDKit {

namespace detail {

// Behavioral differences between Atom and Bond iterators
template <typename ElementPtr>
struct GraphElementTraits;

template <>
struct GraphElementTraits<Atom*> {
  using ElementType = Atom*;

  // We need to wrap 4 underlying BGL iterator types
  using MolPrimaryIterator    = MolGraphTraits<MolGraph>::vertex_iterator;
  using MolSecondaryIterator  = MolGraphTraits<MolGraph>::adjacency_iterator;
  using CSRPrimaryIterator    = MolGraphTraits<CSRMolGraph>::vertex_iterator;
  using CSRSecondaryIterator  = MolGraphTraits<CSRMolGraph>::adjacency_iterator;

  // Retrieve element from graph given an iterator
  template <typename Graph, typename Iter>
  static ElementType getElement(const Graph& graph, Iter iter) {
    return MolGraphTraits<std::decay_t<Graph>>::getAtom(graph, *iter);
  }
};

template <>
struct GraphElementTraits<Bond*> {
  using ElementType = Bond*;

  using MolPrimaryIterator    = MolGraphTraits<MolGraph>::edge_iterator;
  using MolSecondaryIterator  = MolGraphTraits<MolGraph>::out_edge_iterator;
  using CSRPrimaryIterator    = MolGraphTraits<CSRMolGraph>::edge_iterator;
  using CSRSecondaryIterator  = MolGraphTraits<CSRMolGraph>::out_edge_iterator;

  template <typename Graph, typename Iter>
  static ElementType getElement(const Graph& graph, Iter iter) {
    return MolGraphTraits<std::decay_t<Graph>>::getBond(graph, *iter);
  }
};

// One iterator to rule them all
template <typename ElementPtr>
class GraphElementIterator {
private:
  using Traits = GraphElementTraits<ElementPtr>;

public:
  using MolPrimaryIterator   = typename Traits::MolPrimaryIterator;
  using MolSecondaryIterator = typename Traits::MolSecondaryIterator;
  using CSRPrimaryIterator   = typename Traits::CSRPrimaryIterator;
  using CSRSecondaryIterator = typename Traits::CSRSecondaryIterator;

  using IterVariant = std::variant<std::monostate,
                                   MolPrimaryIterator,
                                   MolSecondaryIterator,
                                   CSRPrimaryIterator,
                                   CSRSecondaryIterator>;
  using GraphVariant = std::variant<std::monostate, const MolGraph*, const CSRMolGraph*>;

  using iterator_category = std::input_iterator_tag;
  using difference_type   = std::ptrdiff_t;
  using value_type        = ElementPtr;
  using pointer           = ElementPtr;
  using reference         = ElementPtr;

  // Conditional noexcept traits
  static constexpr bool NoexceptInc =
      noexcept(++std::declval<MolPrimaryIterator&>())   &&
      noexcept(++std::declval<MolSecondaryIterator&>()) &&
      noexcept(++std::declval<CSRPrimaryIterator&>())   &&
      noexcept(++std::declval<CSRSecondaryIterator&>());
  static constexpr bool NoexceptEq =
      noexcept(std::declval<const GraphVariant&>() == std::declval<const GraphVariant&>()) &&
      noexcept(std::declval<const IterVariant&>()  == std::declval<const IterVariant&>());

  GraphElementIterator() = default;
  GraphElementIterator(const GraphElementIterator &) = default;
  GraphElementIterator(GraphElementIterator &&) noexcept = default;
  GraphElementIterator &operator=(const GraphElementIterator &) = default;
  GraphElementIterator &operator=(GraphElementIterator &&) noexcept = default;

  GraphElementIterator &operator++() noexcept(NoexceptInc) {
    assert(!std::holds_alternative<std::monostate>(m_graph) && "incrementing iterator with null graph");
    assert(!std::holds_alternative<std::monostate>(m_iter) && "invalid iterator");
    std::visit([&](auto &it) {
      using T = std::decay_t<decltype(it)>;
      if constexpr (!std::is_same_v<T, std::monostate>) {
        ++it;
      } else {
        assert(false && "incrementing invalid iterator");
      }
    }, m_iter);
    return *this;
  }

  GraphElementIterator operator++(int) {
    GraphElementIterator tmp(*this);
    ++(*this);
    return tmp;
  }

  ElementPtr operator*() const {
    assert(!std::holds_alternative<std::monostate>(m_graph) && "dereferencing iterator with null graph");
    assert(!std::holds_alternative<std::monostate>(m_iter)  && "invalid iterator");
    ElementPtr result = nullptr;
    std::visit(
        [&](auto gptr, auto &it) {
          using GPtr = std::decay_t<decltype(gptr)>;
          using ItT  = std::decay_t<decltype(it)>;
          if constexpr (std::is_same_v<GPtr, const MolGraph*> &&
                       (std::is_same_v<ItT, MolPrimaryIterator> || std::is_same_v<ItT, MolSecondaryIterator>)) {
            result = Traits::getElement(*gptr, it);
          } else if constexpr (std::is_same_v<GPtr, const CSRMolGraph*> &&
                              (std::is_same_v<ItT, CSRPrimaryIterator> || std::is_same_v<ItT, CSRSecondaryIterator>)) {
            result = Traits::getElement(*gptr, it);
          } else {
            assert(false && "iterator type does not match graph");
          }
        },
        m_graph, m_iter);
    return result;
  }

  ElementPtr operator->() const { return operator*(); }

  bool operator==(const GraphElementIterator &other) const noexcept(NoexceptEq) {
    return m_graph == other.m_graph && m_iter == other.m_iter;
  }

  bool operator!=(const GraphElementIterator &other) const noexcept(NoexceptEq) { return !(*this == other); }

  // Generic factory methods
  static GraphElementIterator molPrimary(const MolGraph *graph, MolPrimaryIterator iter) noexcept {
    GraphElementIterator result;
    result.m_graph = graph;
    result.m_iter  = iter;
    return result;
  }

  static GraphElementIterator molSecondary(const MolGraph *graph, MolSecondaryIterator iter) noexcept {
    GraphElementIterator result;
    result.m_graph = graph;
    result.m_iter  = iter;
    return result;
  }

  static GraphElementIterator csrPrimary(const CSRMolGraph *graph, CSRPrimaryIterator iter) noexcept {
    GraphElementIterator result;
    result.m_graph = graph;
    result.m_iter  = iter;
    return result;
  }

  static GraphElementIterator csrSecondary(const CSRMolGraph *graph, CSRSecondaryIterator iter) noexcept {
    GraphElementIterator result;
    result.m_graph = graph;
    result.m_iter  = iter;
    return result;
  }

  // Backward-compatible aliases for Atom* specialization
  template <typename T = ElementPtr>
  static std::enable_if_t<std::is_same_v<T, Atom*>, GraphElementIterator>
  molVertices(const MolGraph *graph, MolPrimaryIterator iter) {
    return molPrimary(graph, iter);
  }

  template <typename T = ElementPtr>
  static std::enable_if_t<std::is_same_v<T, Atom*>, GraphElementIterator>
  molNeighbors(const MolGraph *graph, MolSecondaryIterator iter) {
    return molSecondary(graph, iter);
  }

  template <typename T = ElementPtr>
  static std::enable_if_t<std::is_same_v<T, Atom*>, GraphElementIterator>
  csrVertices(const CSRMolGraph *graph, CSRPrimaryIterator iter) {
    return csrPrimary(graph, iter);
  }

  template <typename T = ElementPtr>
  static std::enable_if_t<std::is_same_v<T, Atom*>, GraphElementIterator>
  csrNeighbors(const CSRMolGraph *graph, CSRSecondaryIterator iter) {
    return csrSecondary(graph, iter);
  }

  // Backward-compatible aliases for Bond* specialization
  template <typename T = ElementPtr>
  static std::enable_if_t<std::is_same_v<T, Bond*>, GraphElementIterator>
  molEdges(const MolGraph *graph, MolPrimaryIterator iter) {
    return molPrimary(graph, iter);
  }

  template <typename T = ElementPtr>
  static std::enable_if_t<std::is_same_v<T, Bond*>, GraphElementIterator>
  molOutEdges(const MolGraph *graph, MolSecondaryIterator iter) {
    return molSecondary(graph, iter);
  }

  template <typename T = ElementPtr>
  static std::enable_if_t<std::is_same_v<T, Bond*>, GraphElementIterator>
  csrEdges(const CSRMolGraph *graph, CSRPrimaryIterator iter) {
    return csrPrimary(graph, iter);
  }

  template <typename T = ElementPtr>
  static std::enable_if_t<std::is_same_v<T, Bond*>, GraphElementIterator>
  csrOutEdges(const CSRMolGraph *graph, CSRSecondaryIterator iter) {
    return csrSecondary(graph, iter);
  }

private:
  IterVariant  m_iter{};
  GraphVariant m_graph{};
};

} // namespace detail

namespace lt {

using AtomIterator = detail::GraphElementIterator<Atom*>;
using BondIterator = detail::GraphElementIterator<Bond*>;

class AtomRange {
private:
  AtomIterator m_begin;
  AtomIterator m_end;

public:
  AtomRange(AtomIterator begin, AtomIterator end) : m_begin(std::move(begin)), m_end(std::move(end)) {}

  AtomIterator begin() const { return m_begin; }
  AtomIterator end()   const { return m_end;   }
};

class BondRange {
private:
  BondIterator m_begin;
  BondIterator m_end;

public:
  BondRange(BondIterator begin, BondIterator end) : m_begin(std::move(begin)), m_end(std::move(end)) {}

  BondIterator begin() const { return m_begin; }
  BondIterator end()   const { return m_end;   }
};

template <typename GraphType>
std::pair<AtomIterator, AtomIterator> createAtomIterators(const GraphType &graph) {
  auto [begin, end] = boost::vertices(graph);
  using Decayed = std::decay_t<GraphType>;

  if constexpr (std::is_same_v<Decayed, MolGraph>) {
    return {
        AtomIterator::molVertices(&graph, begin),
        AtomIterator::molVertices(&graph, end)};
  } else if constexpr (std::is_same_v<Decayed, CSRMolGraph>) {
    return {
        AtomIterator::csrVertices(&graph, begin),
        AtomIterator::csrVertices(&graph, end)};
  } else {
    static_assert(sizeof(GraphType) == 0, "Unsupported graph type for AtomIterator");
  }
}

template <typename GraphType>
std::pair<BondIterator, BondIterator> createBondIterators(const GraphType &graph) {
  auto [begin, end] = boost::edges(graph);
  using Decayed = std::decay_t<GraphType>;

  if constexpr (std::is_same_v<Decayed, MolGraph>) {
    return {
        BondIterator::molEdges(&graph, begin),
        BondIterator::molEdges(&graph, end)};
  } else if constexpr (std::is_same_v<Decayed, CSRMolGraph>) {
    return {
        BondIterator::csrEdges(&graph, begin),
        BondIterator::csrEdges(&graph, end)};
  } else {
    static_assert(sizeof(GraphType) == 0, "Unsupported graph type for BondIterator");
  }
}

} // namespace lt

} // namespace RDKit

#endif // RDKIT_GRAPH_ITERATORS_H
