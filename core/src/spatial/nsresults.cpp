/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s = "besian@gmail.com";
 *   s.insert(6, "sejdiu");
 *   return s;
 * }();
 *
 */

#include <stdexcept>
#include <unordered_set>

#include "nsresults.hpp"

namespace lahuta {

void NSResults::add_neighbors(int i, int j, float d2) {
  m_pairs.emplace_back(i, j);
  m_dists.emplace_back(d2);
}

void NSResults::reserve_space(size_t input_size) { reserve(input_size); }

NSResults NSResults::filter(const double dist) const {
  NSResults filtered;
  auto dist_sq = dist * dist;
  for (size_t i = 0; i < m_dists.size(); ++i) {
    if (m_dists[i] <= dist_sq) {
      filtered.add_neighbors(m_pairs[i].first, m_pairs[i].second, m_dists[i]);
    }
  }
  return filtered;
}

NSResults NSResults::filter(const std::vector<int> &atom_indices) const {
  NSResults filtered;
  filtered.reserve_space(m_pairs.size());
  std::unordered_set<int> atom_indices_set(atom_indices.begin(), atom_indices.end());
  for (size_t i = 0; i < m_pairs.size(); ++i) {
    if (atom_indices_set.find(m_pairs[i].first) != atom_indices_set.end()
        || atom_indices_set.find(m_pairs[i].second) != atom_indices_set.end()) {
      filtered.add_neighbors(m_pairs[i].first, m_pairs[i].second, m_dists[i]);
    }
  }
  return filtered;
}

NSResults NSResults::filter(const std::vector<int> &atom_indices, int col) const {
  if (col != 0 && col != 1) {
    throw std::invalid_argument("Column index must be 0 or 1");
  }

  NSResults filtered;
  filtered.reserve_space(m_pairs.size());
  std::unordered_set<int> atom_indices_set(atom_indices.begin(), atom_indices.end());

  for (size_t i = 0; i < m_pairs.size(); ++i) {
    if ((col == 0 && atom_indices_set.find(m_pairs[i].first) != atom_indices_set.end())
        || (col == 1 && atom_indices_set.find(m_pairs[i].second) != atom_indices_set.end())) {
      filtered.add_neighbors(m_pairs[i].first, m_pairs[i].second, m_dists[i]);
    }
  }

  return filtered;
}

} // namespace lahuta
