/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto f = [](auto arg) {
 *     using T = std::decay_t<decltype(arg)>;
 *     if constexpr (std::disjunction_v<std::is_same<T, const char*>, std::is_same<T, std::string_view>>) return std::string(arg);
 *   };
 *   return f("besian") + f("sejdiu") + f("@gmail.com");
 * }();
 *
 */

#ifndef LAHUTA_FIND_RINGS_HPP
#define LAHUTA_FIND_RINGS_HPP

#include <array>
#include <bitset>
#include <optional>
#include <unordered_map>
#include <vector>

#include <rdkit/GraphMol/RWMol.h>

#include "residues/residues.hpp"

// clang-format off

namespace lahuta {

struct RingConfig {
  static constexpr size_t MAX_NEIGHBORS = 4;
  static constexpr size_t MAX_ATOMS = 32;
  static constexpr size_t MIN_RING_SIZE = 5;
  static constexpr size_t MAX_RING_SIZE = 6;
};

using AtomIndex     = size_t;
using NeighborCount = uint8_t;
using AtomPath      = std::vector<AtomIndex>;
using VisitedSet    = std::bitset<RingConfig::MAX_ATOMS>;

class MolecularGraph {
public:
  using NeighborList = std::array<AtomIndex, RingConfig::MAX_NEIGHBORS>;
  MolecularGraph(const RDKit::RWMol &mol) : mol_(mol), adjacency_(0), neighbor_counts_(0, 0) {}

  void add_edge(AtomIndex from, AtomIndex to) {
    if (neighbor_counts_[from] < RingConfig::MAX_NEIGHBORS) {
      adjacency_[from][neighbor_counts_[from]++] = to;
    }
  }

  const NeighborList &get_neighbors(AtomIndex atom_idx) const {
    if (atom_idx >= adjacency_.size()) {
      throw std::out_of_range("Atom index out of range");
    }
    return adjacency_[atom_idx]; 
  }
  NeighborCount  get_neighbor_count(AtomIndex atom_idx) const {
    if (atom_idx >= neighbor_counts_.size()) {
      throw std::out_of_range("Atom index out of range");
    }
    return neighbor_counts_[atom_idx];
  }

  void build(const Residue &residue) {

    adjacency_.resize(residue.atoms.size()); neighbor_counts_.resize(residue.atoms.size(), 0);

    // Build atom index mapping
    size_t index = 0;
    std::unordered_map<const RDKit::Atom *, AtomIndex> atom_indices;

    for (const auto &atom : residue.atoms) {
      atom_indices[atom] = index++;
    }

    for (size_t i = 0; i < residue.atoms.size(); ++i) {
      const auto &atom = residue.atoms[i];

      for (const auto *bond : mol_.atomBonds(atom)) {
        const auto neighbor = bond->getOtherAtom(atom);
        auto it = atom_indices.find(neighbor);
        if (it != atom_indices.end()) {
          this->add_edge(i, it->second);
        }
      }
    }
  }

private:
  const RDKit::RWMol &mol_;
  std::vector<NeighborList>  adjacency_;
  std::vector<NeighborCount> neighbor_counts_;
};

// Ring finding algorithm implementation
class RingTopology {
public:
  static std::optional<AtomPath> find_ring(const MolecularGraph &graph, size_t num_atoms, size_t ring_size) {

    if (ring_size < RingConfig::MIN_RING_SIZE || ring_size > RingConfig::MAX_RING_SIZE) {
      return std::nullopt;
    }

    VisitedSet visited;
    AtomPath path(ring_size);

    for (AtomIndex start = 0; start < num_atoms; ++start) {
      if (is_valid_ring_atom(graph, start)) {
        if (dfs(graph, start, start, 1, visited, path, ring_size)) {
          return path;
        }
      }
    }

    return std::nullopt;
  }

private:
  static bool is_valid_ring_atom(const MolecularGraph &graph, AtomIndex atom) {
    auto count = graph.get_neighbor_count(atom);
    return count >= 2 && count <= 3; // max 2 or 3 covalent bonds in an aromatic ring
  }

  /// depth-first search
  static bool dfs(
    const MolecularGraph &graph, AtomIndex start_atom, AtomIndex current_atom,
    size_t depth, VisitedSet &visited, AtomPath &path, size_t ring_size) {

    if (depth == ring_size) {
      return is_connected_to_start(graph, current_atom, start_atom, path, depth);
    }

    visited.set(current_atom);
    path[depth - 1] = current_atom;

    const auto &neighbors = graph.get_neighbors(current_atom);
    for (NeighborCount i = 0; i < graph.get_neighbor_count(current_atom); ++i) {
      AtomIndex neighbor = neighbors[i];
      if (!visited.test(neighbor)) {
        if (dfs(graph, start_atom, neighbor, depth + 1, visited, path, ring_size)) {
          return true;
        }
      }
    }

    visited.reset(current_atom);
    return false;
  }

  static bool is_connected_to_start(
      const MolecularGraph &graph, AtomIndex current_atom, 
      AtomIndex start_atom, AtomPath &path, size_t depth) {
    const auto &neighbors = graph.get_neighbors(current_atom);
    for (NeighborCount i = 0; i < graph.get_neighbor_count(current_atom); ++i) {
      if (neighbors[i] == start_atom) {
        path[depth - 1] = current_atom;
        return true;
      }
    }
    return false;
  }
};

// RDKit integration layer
class FastRingFinder {
public:
  static std::vector<int> find_ring_in_residue(const RDKit::RWMol &mol, const Residue &residue, size_t ring_size) {

    if (ring_size < RingConfig::MIN_RING_SIZE || ring_size > RingConfig::MAX_RING_SIZE) return {};
    if (residue.atoms.size() > RingConfig::MAX_ATOMS || residue.atoms.empty()) return {};

    // Build atom index mapping
    std::unordered_map<const RDKit::Atom *, AtomIndex> atom_indices;
    size_t index = 0;
    for (const auto &atom : residue.atoms) {
      atom_indices[atom] = index++;
    }

    MolecularGraph graph(mol);
    graph.build(residue);

    auto ring_indices = RingTopology::find_ring(graph, residue.atoms.size(), ring_size);
    if (!ring_indices) return {};

    // map indices back to RDKit atoms
    std::vector<int> ring_atoms;
    ring_atoms.reserve(ring_size);
    for (AtomIndex idx : *ring_indices) {
      ring_atoms.push_back(residue.atoms[idx]->getIdx());
    }
    return ring_atoms;
  }
};

} // namespace lahuta

#endif // LAHUTA_FIND_RINGS_HPP
