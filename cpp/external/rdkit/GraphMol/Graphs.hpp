#ifndef RD_GRAPHS_H
#define RD_GRAPHS_H

/// Std stuff
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include <boost/graph/connected_components.hpp>
#include <cstddef>
#include <utility>

// boost stuff
#include <RDGeneral/BoostStartInclude.h>
#include <boost/dynamic_bitset.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/smart_ptr.hpp>

#include <RDGeneral/BoostEndInclude.h>

// our stuff
#include "Atom.h"
#include "Bond.h"
#include "GraphIterators.hpp"
#include <RDGeneral/RDProps.h>
#include <RDGeneral/types.h>

namespace RDKit {

// clang-format off
class GraphImpl {
public:
  virtual ~GraphImpl() = default;

  using vertex_descriptor  = boost::graph_traits<MolGraph>::vertex_descriptor;
  using edge_descriptor    = boost::graph_traits<MolGraph>::edge_descriptor;
  using out_edge_iterator  = boost::graph_traits<MolGraph>::out_edge_iterator;

  virtual unsigned int getNumVertices() const = 0;
  virtual unsigned int getNumEdges()    const = 0;

  virtual Atom *getAtom(unsigned int idx) const = 0;
  virtual Bond *getBond(unsigned int idx) const = 0;

  virtual void setBond(unsigned int idx, Bond *bond) = 0;
  virtual void setAtom(unsigned int idx, Atom *atom) = 0;

  // used by RWMol::replaceBond
  virtual void setBond(boost::graph_traits<MolGraph>::edge_descriptor edge, Bond *bond) {
    throw std::runtime_error("not implemented");
  }

  // used for operator[] and by RWMol
  virtual Bond *getBond(const edge_descriptor &edge) const { throw std::runtime_error("not implemented"); }
  virtual void removeAtom(unsigned int idx) { throw std::runtime_error("not implemented"); }

  virtual unsigned int addAtom(Atom *atom) = 0;
  virtual void addBond(unsigned int src, unsigned int dst, Bond *bond) = 0;
  virtual void removeBond(unsigned int src, unsigned int dst) = 0;

  virtual unsigned int getAtomDegree(unsigned int idx) const = 0;

  // iterators
  virtual std::pair<lt::AtomIterator, lt::AtomIterator> getAtomIterators() const = 0;
  virtual std::pair<lt::BondIterator, lt::BondIterator> getBondIterators() const = 0;
  virtual std::pair<lt::BondIterator, lt::BondIterator> getAtomBonds(unsigned int atom_idx) const = 0;
  virtual std::pair<lt::AtomIterator, lt::AtomIterator> getAtomNeighbors(unsigned int atom_idx) const = 0;

  virtual std::pair<Bond *, bool> getBondBetweenAtoms(unsigned int src, unsigned int dst) const = 0;

  virtual int getConnectedComponents(std::vector<int> &mapping) const = 0;

  virtual bool isMolGraph() const { return false; }
  virtual bool isCSRMolGraph() const { return false; }

  virtual std::unique_ptr<GraphImpl> clone() const = 0;
  virtual void clear() = 0;

  virtual std::pair<
      boost::graph_traits<MolGraph>::out_edge_iterator,
      boost::graph_traits<MolGraph>::out_edge_iterator>
  getMolGraphAtomBonds(unsigned int idx) const {
    throw std::runtime_error("Not implemented");
  }

  virtual std::pair<
      boost::graph_traits<CSRMolGraph>::out_edge_iterator,
      boost::graph_traits<CSRMolGraph>::out_edge_iterator>
  getCSRMolGraphAtomBonds(unsigned int idx) const {
    throw std::runtime_error("Not implemented");
  }

  virtual std::pair<
      boost::graph_traits<MolGraph>::adjacency_iterator,
      boost::graph_traits<MolGraph>::adjacency_iterator>
  getMolGraphAtomNeighbors(unsigned int idx) const {
    throw std::runtime_error("Not implemented");
  }

  virtual std::pair<
      boost::graph_traits<CSRMolGraph>::adjacency_iterator,
      boost::graph_traits<CSRMolGraph>::adjacency_iterator>
  getCSRMolGraphAtomNeighbors(unsigned int idx) const {
    throw std::runtime_error("Not implemented");
  }
};

class MolGraphImpl : public GraphImpl {
private:
  MolGraph m_graph;

public:
  MolGraphImpl() = default;
  explicit MolGraphImpl(const MolGraph &graph) : m_graph(graph) {}

  bool isMolGraph() const override { return true; }

  unsigned int getNumVertices() const override { return boost::num_vertices(m_graph); }
  unsigned int getNumEdges()    const override { return boost::num_edges(m_graph); }

  Atom *getAtom(unsigned int idx) const override {
    return m_graph[static_cast<MolGraphTraits<MolGraph>::vertex_descriptor>(idx)];
  }

  void setAtom(unsigned int idx, Atom *atom) override {
    m_graph[static_cast<MolGraphTraits<MolGraph>::vertex_descriptor>(idx)] = atom;
  }

  Bond *getBond(unsigned int idx) const override {
    auto [iter, end] = boost::edges(m_graph);
    for (unsigned int i = 0; i < idx; i++) {
      ++iter;
    }
    return m_graph[*iter];
  }

  void setBond(unsigned int idx, Bond *bond) override {
    auto [iter, end] = boost::edges(m_graph);
    for (unsigned int i = 0; i < idx; i++) {
      ++iter;
    }
    m_graph[*iter] = bond;
  }

  Bond *getBond(const edge_descriptor &edge) const override { return m_graph[edge]; }

  unsigned int addAtom(Atom *atom) override {
    auto v = boost::add_vertex(m_graph);
    m_graph[v] = atom;
    return v;
  }

  void addBond(unsigned int src, unsigned int dst, Bond *bond) override {
    auto result = boost::add_edge(src, dst, m_graph);
    CHECK_INVARIANT(result.second, "bond could not be added");
    m_graph[result.first] = bond;
  }

  void removeAtom(unsigned int idx) override {
    boost::clear_vertex(static_cast<MolGraphTraits<MolGraph>::vertex_descriptor>(idx), m_graph);
  }

  void removeBond(unsigned int src, unsigned int dst) override { boost::remove_edge(src, dst, m_graph); }

  unsigned int getAtomDegree(unsigned int idx) const override {
    return boost::degree(static_cast<MolGraphTraits<MolGraph>::vertex_descriptor>(idx), m_graph);
  }

  std::pair<lt::AtomIterator, lt::AtomIterator> getAtomIterators() const override {
    return lt::createAtomIterators(m_graph);
  }

  std::pair<lt::BondIterator, lt::BondIterator> getBondIterators() const override {
    return lt::createBondIterators(m_graph);
  }

  std::pair<lt::BondIterator, lt::BondIterator> getAtomBonds(unsigned int atom_idx) const override {
    auto vertex = static_cast<MolGraphTraits<MolGraph>::vertex_descriptor>(atom_idx);
    auto [begin, end] = boost::out_edges(vertex, m_graph);

    return {
        lt::BondIterator(new BondIteratorImpl<MolGraph, decltype(begin)>(&m_graph, begin)),
        lt::BondIterator(new BondIteratorImpl<MolGraph, decltype(end)>(&m_graph, end))};
  }

  std::pair<lt::AtomIterator, lt::AtomIterator> getAtomNeighbors(unsigned int atom_idx) const override {
    auto vertex = static_cast<MolGraphTraits<MolGraph>::vertex_descriptor>(atom_idx);
    auto [begin, end] = boost::adjacent_vertices(vertex, m_graph);

    return {
        lt::AtomIterator(new AtomIteratorImpl<MolGraph, decltype(begin)>(&m_graph, begin)),
        lt::AtomIterator(new AtomIteratorImpl<MolGraph, decltype(end)>(&m_graph, end))};
  }

  std::pair<
      boost::graph_traits<MolGraph>::out_edge_iterator,
      boost::graph_traits<MolGraph>::out_edge_iterator>
  getMolGraphAtomBonds(unsigned int idx) const override {
    return boost::out_edges(static_cast<boost::graph_traits<MolGraph>::vertex_descriptor>(idx), m_graph);
  }

  std::pair<
      boost::graph_traits<MolGraph>::adjacency_iterator,
      boost::graph_traits<MolGraph>::adjacency_iterator>
  getMolGraphAtomNeighbors(unsigned int idx) const override {
    return boost::adjacent_vertices(
        static_cast<boost::graph_traits<MolGraph>::vertex_descriptor>(idx),
        m_graph);
  }

  std::pair<Bond *, bool> getBondBetweenAtoms(unsigned int src, unsigned int dst) const override {
    auto [e, success] = boost::edge(src, dst, m_graph);
    if (success) return {m_graph[e], success};
    return {nullptr, false};
  }

  int getConnectedComponents(std::vector<int> &mapping) const override {
    return boost::connected_components(m_graph, &mapping[0]);
  }

  std::unique_ptr<GraphImpl> clone() const override { return std::make_unique<MolGraphImpl>(m_graph); }

  void clear() override { m_graph.clear(); }

  const MolGraph &getGraph() const { return m_graph; }
  MolGraph &getGraph() { return m_graph; }

  void build(const std::vector<Atom *> &atoms, const std::vector<std::tuple<size_t, size_t, Bond *>> &bonds) {
    m_graph.clear();

    for (Atom *atom : atoms) {
      auto v = boost::add_vertex(m_graph);
      m_graph[v] = atom;
      atom->setIdx(v);
    }

    for (const auto &[src, dst, bond] : bonds) {
      auto [e, success] = boost::add_edge(src, dst, m_graph);
      m_graph[e] = bond;
    }
  }
};

class CSRMolGraphImpl : public GraphImpl {
private:
  CSRMolGraph m_graph;

public:
  CSRMolGraphImpl() {
    std::vector<std::pair<size_t, size_t>> edge_list;
    std::vector<Bond *> edge_props;
    unsigned int n_vertices = 0;
    m_graph = CSRMolGraph(
        boost::edges_are_sorted,
        edge_list.begin(),
        edge_list.end(),
        edge_props.begin(),
        n_vertices);
  }

  explicit CSRMolGraphImpl(const CSRMolGraph &graph) : m_graph(graph) {}

  unsigned int getNumVertices() const override { return boost::num_vertices(m_graph); }
  unsigned int getNumEdges()    const override { return boost::num_edges(m_graph); }

  Atom *getAtom(unsigned int idx) const override {
    return m_graph[static_cast<MolGraphTraits<CSRMolGraph>::vertex_descriptor>(idx)];
  }

  void setAtom(unsigned int idx, Atom *atom) override {
    m_graph[static_cast<MolGraphTraits<CSRMolGraph>::vertex_descriptor>(idx)] = atom;
  }

  Bond *getBond(unsigned int idx) const override {
    auto [iter, end] = boost::edges(m_graph);
    for (unsigned int i = 0; i < idx; i++) {
      ++iter;
    }
    return m_graph[*iter];
  }

  void setBond(unsigned int idx, Bond *bond) override {
    throw std::runtime_error("Cannot modify compressed_sparse_row_graph after construction");
  }

  unsigned int addAtom(Atom *atom) override {
    throw std::runtime_error("Cannot modify compressed_sparse_row_graph after construction");
  }

  void addBond(unsigned int src, unsigned int dst, Bond *bond) override {
    throw std::runtime_error("Cannot modify compressed_sparse_row_graph after construction");
  }

  void removeBond(unsigned int src, unsigned int dst) override {
    throw std::runtime_error("Cannot modify compressed_sparse_row_graph after construction");
  }

  unsigned int getAtomDegree(unsigned int idx) const override {
    return boost::out_degree(static_cast<MolGraphTraits<CSRMolGraph>::vertex_descriptor>(idx), m_graph);
  }

  std::pair<lt::AtomIterator, lt::AtomIterator> getAtomIterators() const override {
    return lt::createAtomIterators(m_graph);
  }

  std::pair<lt::BondIterator, lt::BondIterator> getBondIterators() const override {
    return lt::createBondIterators(m_graph);
  }

  std::pair<lt::BondIterator, lt::BondIterator> getAtomBonds(unsigned int atom_idx) const override {
    auto vertex = static_cast<MolGraphTraits<CSRMolGraph>::vertex_descriptor>(atom_idx);
    auto [begin, end] = boost::out_edges(vertex, m_graph);

    return {
        lt::BondIterator(new BondIteratorImpl<CSRMolGraph, decltype(begin)>(&m_graph, begin)),
        lt::BondIterator(new BondIteratorImpl<CSRMolGraph, decltype(end)>(&m_graph, end))};
  }

  std::pair<lt::AtomIterator, lt::AtomIterator> getAtomNeighbors(unsigned int atom_idx) const override {
    auto vertex = static_cast<MolGraphTraits<CSRMolGraph>::vertex_descriptor>(atom_idx);
    auto [begin, end] = boost::adjacent_vertices(vertex, m_graph);

    return {
        lt::AtomIterator(new AtomIteratorImpl<CSRMolGraph, decltype(begin)>(&m_graph, begin)),
        lt::AtomIterator(new AtomIteratorImpl<CSRMolGraph, decltype(end)>(&m_graph, end))};
  }

  std::pair<Bond *, bool> getBondBetweenAtoms(unsigned int src, unsigned int dst) const override {
    auto [e, success] = boost::edge(
        static_cast<MolGraphTraits<CSRMolGraph>::vertex_descriptor>(src),
        static_cast<MolGraphTraits<CSRMolGraph>::vertex_descriptor>(dst),
        m_graph);
    if (success) return {m_graph[e], success};
    return {nullptr, false};
  }

  int getConnectedComponents(std::vector<int> &mapping) const override {
    return boost::connected_components(m_graph, mapping.data());
  }

  std::pair<
      boost::graph_traits<CSRMolGraph>::out_edge_iterator,
      boost::graph_traits<CSRMolGraph>::out_edge_iterator>
  getCSRMolGraphAtomBonds(unsigned int idx) const override {
    return boost::out_edges(static_cast<boost::graph_traits<CSRMolGraph>::vertex_descriptor>(idx), m_graph);
  }

  std::pair<
      boost::graph_traits<CSRMolGraph>::adjacency_iterator,
      boost::graph_traits<CSRMolGraph>::adjacency_iterator>
  getCSRMolGraphAtomNeighbors(unsigned int idx) const override {
    return boost::adjacent_vertices(
        static_cast<boost::graph_traits<CSRMolGraph>::vertex_descriptor>(idx),
        m_graph);
  }

  std::unique_ptr<GraphImpl> clone() const override { return std::make_unique<CSRMolGraphImpl>(m_graph); }

  void clear() override {
    std::vector<std::pair<size_t, size_t>> emptyEdgeList;
    std::vector<Bond *> emptyEdgeProps;
    m_graph = CSRMolGraph(
        boost::edges_are_sorted,
        emptyEdgeList.begin(),
        emptyEdgeList.end(),
        emptyEdgeProps.begin(),
        0);
  }

  const CSRMolGraph &getGraph() const { return m_graph; }
  CSRMolGraph &getGraph() { return m_graph; }

  void build(
      const std::vector<Atom *> &atoms, const std::vector<std::pair<size_t, size_t>> &edge_list,
      const std::vector<Bond *> &edge_props) {

    m_graph = CSRMolGraph(
        boost::edges_are_sorted,
        edge_list.begin(),
        edge_list.end(),
        edge_props.begin(),
        atoms.size());

    for (size_t i = 0; i < atoms.size(); ++i) {
      atoms[i]->setIdx(i);
      m_graph[i] = atoms[i];
    }
  }
};

} // namespace RDKit

#endif // RD_GRAPHS_H
