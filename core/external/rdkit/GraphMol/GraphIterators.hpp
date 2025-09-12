#ifndef RDKIT_GRAPH_ITERATORS_H
#define RDKIT_GRAPH_ITERATORS_H

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include <boost/graph/graph_traits.hpp>
#include "GraphDefs.hpp"

// clang-format off
namespace RDKit {

class IAtomIterator {
public:
  virtual ~IAtomIterator() = default;
  virtual IAtomIterator *clone() const = 0;
  virtual void increment() = 0;
  virtual Atom *dereference() const = 0;
  virtual bool equal(const IAtomIterator *other) const = 0;
};

class IBondIterator {
public:
  virtual ~IBondIterator() = default;
  virtual IBondIterator *clone() const = 0;
  virtual void increment() = 0;
  virtual Bond *dereference() const = 0;
  virtual bool equal(const IBondIterator *other) const = 0;
};

template <typename GraphType, typename VertexIter>
class AtomIteratorImpl : public IAtomIterator {
private:
  const GraphType *m_graph;
  VertexIter m_iter;

public:
  AtomIteratorImpl(const GraphType *g, VertexIter it) : m_graph(g), m_iter(it) {}

  IAtomIterator *clone() const override { return new AtomIteratorImpl(m_graph, m_iter); }
  void increment() override { ++m_iter; }
  Atom *dereference() const override { return MolGraphTraits<GraphType>::getAtom(*m_graph, *m_iter); }
  bool equal(const IAtomIterator *other) const override {
    auto *typed_other = dynamic_cast<const AtomIteratorImpl<GraphType, VertexIter> *>(other);
    return typed_other && m_iter == typed_other->m_iter;
  }
};

template <typename GraphType, typename EdgeIter>
class BondIteratorImpl : public IBondIterator {
private:
  const GraphType *m_graph;
  EdgeIter m_iter;

public:
  BondIteratorImpl(const GraphType *g, EdgeIter it) : m_graph(g), m_iter(it) {}

  IBondIterator *clone() const override { return new BondIteratorImpl(m_graph, m_iter); }
  void increment() override { ++m_iter; }
  Bond *dereference() const override { return MolGraphTraits<GraphType>::getBond(*m_graph, *m_iter); }
  bool equal(const IBondIterator *other) const override {
    auto *typed_other = dynamic_cast<const BondIteratorImpl<GraphType, EdgeIter> *>(other);
    return typed_other && m_iter == typed_other->m_iter;
  }
};

namespace lt {

class AtomIterator {
private:
  std::unique_ptr<IAtomIterator> m_impl;

public:
  explicit AtomIterator(IAtomIterator *impl) : m_impl(impl) {}

  AtomIterator(const AtomIterator &other) : m_impl(other.m_impl->clone()) {}

  AtomIterator &operator=(const AtomIterator &other) {
    if (this != &other) {
      m_impl.reset(other.m_impl->clone());
    }
    return *this;
  }

  AtomIterator &operator++() {
    m_impl->increment();
    return *this;
  }

  Atom *operator*() const { return m_impl->dereference(); }
  bool operator==(const AtomIterator &other) const { return m_impl->equal(other.m_impl.get()); }
  bool operator!=(const AtomIterator &other) const { return !(*this == other); }
};

class BondIterator {
private:
  std::unique_ptr<IBondIterator> m_impl;

public:
  explicit BondIterator(IBondIterator *impl) : m_impl(impl) {}

  BondIterator(const BondIterator &other) : m_impl(other.m_impl->clone()) {}

  BondIterator &operator=(const BondIterator &other) {
    if (this != &other) {
      m_impl.reset(other.m_impl->clone());
    }
    return *this;
  }

  BondIterator &operator++() {
    m_impl->increment();
    return *this;
  }

  Bond *operator*() const { return m_impl->dereference(); }
  bool operator==(const BondIterator &other) const { return m_impl->equal(other.m_impl.get()); }
  bool operator!=(const BondIterator &other) const { return !(*this == other); }
};

// Iterator range for atoms
class AtomRange {
private:
  AtomIterator m_begin;
  AtomIterator m_end;

public:
  AtomRange(AtomIterator begin, AtomIterator end) : m_begin(std::move(begin)), m_end(std::move(end)) {}

  AtomIterator begin() const { return m_begin; }
  AtomIterator end()   const { return m_end; }
};

// Iterator range for bonds
class BondRange {
private:
  BondIterator m_begin;
  BondIterator m_end;

public:
  BondRange(BondIterator begin, BondIterator end) : m_begin(std::move(begin)), m_end(std::move(end)) {}

  BondIterator begin() const { return m_begin; }
  BondIterator end()   const { return m_end; }
};

// create iterators
template <typename GraphType>
std::pair<AtomIterator, AtomIterator> createAtomIterators(const GraphType &graph) {
  auto [begin, end] = boost::vertices(graph);
  return {
      AtomIterator(new AtomIteratorImpl<GraphType, decltype(begin)>(&graph, begin)),
      AtomIterator(new AtomIteratorImpl<GraphType, decltype(end)>(&graph, end))};
}

template <typename GraphType>
std::pair<BondIterator, BondIterator> createBondIterators(const GraphType &graph) {
  auto [begin, end] = boost::edges(graph);
  return {
      BondIterator(new BondIteratorImpl<GraphType, decltype(begin)>(&graph, begin)),
      BondIterator(new BondIteratorImpl<GraphType, decltype(end)>(&graph, end))};
}

} // namespace lt

} // namespace RDKit

#endif // RDKIT_GRAPH_ITERATORS_H
