#ifndef RDKIT_GRAPH_ITERATORS_H
#define RDKIT_GRAPH_ITERATORS_H

#include <cassert>
#include <cstdint>
#include <iterator>
#include <memory>
#include <type_traits>
#include <utility>

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
      noexcept(std::declval<const MolPrimaryIterator&>()   == std::declval<const MolPrimaryIterator&>())   &&
      noexcept(std::declval<const MolSecondaryIterator&>() == std::declval<const MolSecondaryIterator&>()) &&
      noexcept(std::declval<const CSRPrimaryIterator&>()   == std::declval<const CSRPrimaryIterator&>())   &&
      noexcept(std::declval<const CSRSecondaryIterator&>() == std::declval<const CSRSecondaryIterator&>());

  enum class IteratorKind : std::uint8_t {
    None = 0,
    MolPrimary,
    MolSecondary,
    CSRPrimary,
    CSRSecondary,
  };

  union IteratorStorage {
    MolPrimaryIterator   molPrimary;
    MolSecondaryIterator molSecondary;
    CSRPrimaryIterator   csrPrimary;
    CSRSecondaryIterator csrSecondary;

    IteratorStorage() {}
    ~IteratorStorage() {}
  };

  GraphElementIterator() = default;
  ~GraphElementIterator() { reset(); }

  GraphElementIterator(const GraphElementIterator &other) { copyFrom(other); }
  GraphElementIterator(GraphElementIterator &&other) noexcept { moveFrom(std::move(other)); }

  GraphElementIterator &operator=(const GraphElementIterator &other) {
    if (this != &other) {
      reset();
      copyFrom(other);
    }
    return *this;
  }

  GraphElementIterator &operator=(GraphElementIterator &&other) noexcept {
    if (this != &other) {
      reset();
      moveFrom(std::move(other));
    }
    return *this;
  }

  GraphElementIterator &operator++() noexcept(NoexceptInc) {
    assert(m_kind != IteratorKind::None && "incrementing invalid iterator");
    switch (m_kind) {
      case IteratorKind::MolPrimary:
        ++m_iter.molPrimary;
        break;
      case IteratorKind::MolSecondary:
        ++m_iter.molSecondary;
        break;
      case IteratorKind::CSRPrimary:
        ++m_iter.csrPrimary;
        break;
      case IteratorKind::CSRSecondary:
        ++m_iter.csrSecondary;
        break;
      case IteratorKind::None:
        assert(false && "incrementing invalid iterator");
        break;
    }
    return *this;
  }

  GraphElementIterator operator++(int) {
    GraphElementIterator tmp(*this);
    ++(*this);
    return tmp;
  }

  ElementPtr operator*() const {
    assert(m_kind != IteratorKind::None && "dereferencing invalid iterator");
    switch (m_kind) {
      case IteratorKind::MolPrimary:   return Traits::getElement(*static_cast<const MolGraph*>(m_graph), m_iter.molPrimary);
      case IteratorKind::MolSecondary: return Traits::getElement(*static_cast<const MolGraph*>(m_graph), m_iter.molSecondary);
      case IteratorKind::CSRPrimary:   return Traits::getElement(*static_cast<const CSRMolGraph*>(m_graph), m_iter.csrPrimary);
      case IteratorKind::CSRSecondary: return Traits::getElement(*static_cast<const CSRMolGraph*>(m_graph), m_iter.csrSecondary);
      case IteratorKind::None:
        assert(false && "dereferencing invalid iterator");
        return nullptr;
    }
    return nullptr;
  }

  ElementPtr operator->() const { return operator*(); }

  bool operator==(const GraphElementIterator &other) const noexcept(NoexceptEq) {
    if (m_kind != other.m_kind || m_graph != other.m_graph) return false;
    switch (m_kind) {
      case IteratorKind::None:
        return true;
      case IteratorKind::MolPrimary:   return m_iter.molPrimary == other.m_iter.molPrimary;
      case IteratorKind::MolSecondary: return m_iter.molSecondary == other.m_iter.molSecondary;
      case IteratorKind::CSRPrimary:   return m_iter.csrPrimary == other.m_iter.csrPrimary;
      case IteratorKind::CSRSecondary: return m_iter.csrSecondary == other.m_iter.csrSecondary;
    }
    return false;
  }

  bool operator!=(const GraphElementIterator &other) const noexcept(NoexceptEq) { return !(*this == other); }

  // Generic factory methods
  static GraphElementIterator molPrimary(const MolGraph *graph, MolPrimaryIterator iter) noexcept {
    GraphElementIterator result;
    result.m_graph = graph;
    result.m_kind  = IteratorKind::MolPrimary;
    ::new (&result.m_iter.molPrimary) MolPrimaryIterator(std::move(iter));
    return result;
  }

  static GraphElementIterator molSecondary(const MolGraph *graph, MolSecondaryIterator iter) noexcept {
    GraphElementIterator result;
    result.m_graph = graph;
    result.m_kind  = IteratorKind::MolSecondary;
    ::new (&result.m_iter.molSecondary) MolSecondaryIterator(std::move(iter));
    return result;
  }

  static GraphElementIterator csrPrimary(const CSRMolGraph *graph, CSRPrimaryIterator iter) noexcept {
    GraphElementIterator result;
    result.m_graph = graph;
    result.m_kind  = IteratorKind::CSRPrimary;
    ::new (&result.m_iter.csrPrimary) CSRPrimaryIterator(std::move(iter));
    return result;
  }

  static GraphElementIterator csrSecondary(const CSRMolGraph *graph, CSRSecondaryIterator iter) noexcept {
    GraphElementIterator result;
    result.m_graph = graph;
    result.m_kind  = IteratorKind::CSRSecondary;
    ::new (&result.m_iter.csrSecondary) CSRSecondaryIterator(std::move(iter));
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
  void reset() noexcept {
    switch (m_kind) {
      case IteratorKind::MolPrimary:
        std::destroy_at(std::addressof(m_iter.molPrimary));
        break;
      case IteratorKind::MolSecondary:
        std::destroy_at(std::addressof(m_iter.molSecondary));
        break;
      case IteratorKind::CSRPrimary:
        std::destroy_at(std::addressof(m_iter.csrPrimary));
        break;
      case IteratorKind::CSRSecondary:
        std::destroy_at(std::addressof(m_iter.csrSecondary));
        break;
      case IteratorKind::None:
        break;
    }
    m_kind  = IteratorKind::None;
    m_graph = nullptr;
  }

  void copyFrom(const GraphElementIterator &other) {
    m_kind  = other.m_kind;
    m_graph = other.m_graph;
    switch (m_kind) {
      case IteratorKind::MolPrimary:
        ::new (&m_iter.molPrimary) MolPrimaryIterator(other.m_iter.molPrimary);
        break;
      case IteratorKind::MolSecondary:
        ::new (&m_iter.molSecondary) MolSecondaryIterator(other.m_iter.molSecondary);
        break;
      case IteratorKind::CSRPrimary:
        ::new (&m_iter.csrPrimary) CSRPrimaryIterator(other.m_iter.csrPrimary);
        break;
      case IteratorKind::CSRSecondary:
        ::new (&m_iter.csrSecondary) CSRSecondaryIterator(other.m_iter.csrSecondary);
        break;
      case IteratorKind::None:
        break;
    }
  }

  void moveFrom(GraphElementIterator &&other) noexcept {
    m_kind  = other.m_kind;
    m_graph = other.m_graph;
    switch (m_kind) {
      case IteratorKind::MolPrimary:
        ::new (&m_iter.molPrimary) MolPrimaryIterator(std::move(other.m_iter.molPrimary));
        break;
      case IteratorKind::MolSecondary:
        ::new (&m_iter.molSecondary) MolSecondaryIterator(std::move(other.m_iter.molSecondary));
        break;
      case IteratorKind::CSRPrimary:
        ::new (&m_iter.csrPrimary) CSRPrimaryIterator(std::move(other.m_iter.csrPrimary));
        break;
      case IteratorKind::CSRSecondary:
        ::new (&m_iter.csrSecondary) CSRSecondaryIterator(std::move(other.m_iter.csrSecondary));
        break;
      case IteratorKind::None:
        break;
    }
    other.reset();
  }

  const void *m_graph = nullptr;
  IteratorKind m_kind = IteratorKind::None;
  IteratorStorage m_iter;
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
