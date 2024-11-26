#ifndef LAHUTA_FIND_RINGS_HPP
#define LAHUTA_FIND_RINGS_HPP

#include "GraphMol/RWMol.h"
#include "residues.hpp"
#include <array>
#include <bitset>
#include <optional>
#include <unordered_map>
#include <vector>

namespace lahuta {

struct RingConfig {
  static constexpr size_t MAX_NEIGHBORS = 4;
  static constexpr size_t MAX_ATOMS = 32;
  static constexpr size_t MIN_RING_SIZE = 5;
  static constexpr size_t MAX_RING_SIZE = 6;
};

using AtomIndex = size_t;
using NeighborCount = uint8_t;
using AtomPath = std::vector<AtomIndex>;
using VisitedSet = std::bitset<RingConfig::MAX_ATOMS>;

class MolecularGraph {
public:
  using NeighborList = std::array<AtomIndex, RingConfig::MAX_NEIGHBORS>;

  MolecularGraph(size_t num_atoms) : adjacency_(num_atoms), neighbor_counts_(num_atoms, 0) {}

  void add_edge(AtomIndex from, AtomIndex to) {
    if (neighbor_counts_[from] < RingConfig::MAX_NEIGHBORS) {
      adjacency_[from][neighbor_counts_[from]++] = to;
    }
  }

  const NeighborList &get_neighbors(AtomIndex atom) const { return adjacency_[atom]; }
  NeighborCount get_neighbor_count(AtomIndex atom) const { return neighbor_counts_[atom]; }

private:
  std::vector<NeighborList> adjacency_;
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
    return count >= 2 && count <= 3; // chemistry constraint
  }

  // clang-format off
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
// clang-format on

// RDKit integration layer
class FastRingFinder {
public:
  static std::vector<const RDKit::Atom *>
  find_ring_in_residue(const RDKit::RWMol &mol, const Residue &residue, size_t ring_size) {

    if (ring_size < RingConfig::MIN_RING_SIZE || ring_size > RingConfig::MAX_RING_SIZE) return {};
    if (residue.atoms.size() > RingConfig::MAX_ATOMS) return {};

    // Build atom index mapping
    std::unordered_map<const RDKit::Atom *, AtomIndex> atom_indices;
    size_t index = 0;
    for (const auto &atom : residue.atoms) {
      atom_indices[atom] = index++;
    }

    MolecularGraph graph(residue.atoms.size());
    build_graph(mol, residue, atom_indices, graph);

    auto ring_indices = RingTopology::find_ring(graph, residue.atoms.size(), ring_size);

    if (!ring_indices) {
      return {};
    }

    // map indices back to RDKit atoms
    std::vector<const RDKit::Atom *> ring_atoms;
    ring_atoms.reserve(ring_size);
    for (AtomIndex idx : *ring_indices) {
      ring_atoms.push_back(residue.atoms[idx]);
    }
    return ring_atoms;
  }

private:
  static void build_graph(
      const RDKit::RWMol &mol, const Residue &residue,
      const std::unordered_map<const RDKit::Atom *, AtomIndex> &atom_indices, MolecularGraph &graph) {

    for (size_t i = 0; i < residue.atoms.size(); ++i) {
      const auto &atom = residue.atoms[i];
      for (auto bondIt = mol.getAtomBonds(atom); bondIt.first != bondIt.second; ++bondIt.first) {
        const RDKit::Bond *bond = mol[*bondIt.first];
        const RDKit::Atom *neighbor = bond->getOtherAtom(atom);
        auto it = atom_indices.find(neighbor);
        if (it != atom_indices.end()) {
          graph.add_edge(i, it->second);
        }
      }
    }
  }
};

} // namespace lahuta

#endif // LAHUTA_FIND_RINGS_HPP
