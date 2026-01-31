#ifndef LAHUTA_ANALYSIS_DSSP_SECONDARY_HPP
#define LAHUTA_ANALYSIS_DSSP_SECONDARY_HPP

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

#include "analysis/dssp/hbond.hpp"
#include "analysis/dssp/precompute.hpp"
#include "distances/neighbors.hpp"

namespace lahuta::analysis {

enum class BridgeType : std::uint8_t { None = 0, Parallel = 1, AntiParallel = 2 };

struct Bridge {
  BridgeType type = BridgeType::None;
  std::vector<std::uint32_t> i;
  std::vector<std::uint32_t> j;
  std::string chain_i;
  std::string chain_j;
  std::uint32_t ladder = 0;
  std::uint32_t sheet  = 0;
  std::vector<std::size_t> link_indices; // indices into bridges vector
  bool removed = false;                  // mark for deferred removal

  [[nodiscard]] bool operator<(const Bridge &other) const noexcept {
    if (i.empty() || other.i.empty()) return i.size() < other.i.size();
    if (i.front() != other.i.front()) return i.front() < other.i.front();
    if (j.empty() || other.j.empty()) return j.size() < other.j.size();
    if (j.front() != other.j.front()) return j.front() < other.j.front();
    return static_cast<int>(type) < static_cast<int>(other.type);
  }
};

[[nodiscard]] inline BridgeType test_bridge(const std::vector<DsspResidue> &residues, std::size_t i,
                                            std::size_t j) noexcept {
  const auto &r1 = residues[i];
  const auto &r2 = residues[j];
  const int a    = r1.prev;
  const int b    = static_cast<int>(i);
  const int c    = r1.next;
  const int d    = r2.prev;
  const int e    = static_cast<int>(j);
  const int f    = r2.next;

  if (a >= 0 && c >= 0 && d >= 0 && f >= 0) {
    if (detail::no_chain_break(residues, static_cast<std::size_t>(a), static_cast<std::size_t>(c)) &&
        detail::no_chain_break(residues, static_cast<std::size_t>(d), static_cast<std::size_t>(f))) {
      const auto &ra = residues[static_cast<std::size_t>(a)];
      const auto &rb = residues[static_cast<std::size_t>(b)];
      const auto &rc = residues[static_cast<std::size_t>(c)];
      const auto &rd = residues[static_cast<std::size_t>(d)];
      const auto &re = residues[static_cast<std::size_t>(e)];
      const auto &rf = residues[static_cast<std::size_t>(f)];

      if ((test_bond(rc, re) && test_bond(re, ra)) || (test_bond(rf, rb) && test_bond(rb, rd))) {
        return BridgeType::Parallel;
      }
      if ((test_bond(rc, rd) && test_bond(rf, ra)) || (test_bond(re, rb) && test_bond(rb, re))) {
        return BridgeType::AntiParallel;
      }
    }
  }
  return BridgeType::None;
}

[[nodiscard]] inline std::vector<std::pair<std::uint32_t, std::uint32_t>>
collect_ca_pairs(const std::vector<DsspResidue> &residues) {
  RDGeom::POINT3D_VECT ca_coords;
  ca_coords.reserve(residues.size());
  for (const auto &res : residues) {
    ca_coords.emplace_back(res.ca);
  }

  dist::NeighborSearchOptions opts;
  opts.cutoff               = DsspMinimalCADistance;
  opts.brute_force_fallback = true;
  opts.sort_output          = false;

  const auto results = dist::neighbors_within_radius_self(ca_coords, opts);
  std::vector<std::pair<std::uint32_t, std::uint32_t>> pairs;
  pairs.reserve(results.get_pairs().size());
  for (const auto &pair : results.get_pairs()) {
    const int i = pair.first;
    const int j = pair.second;
    if (i < 0 || j < 0) continue;
    if (i >= static_cast<int>(residues.size()) || j >= static_cast<int>(residues.size())) continue;
    if (i == j) continue;
    const auto a = static_cast<std::uint32_t>(std::min(i, j));
    const auto b = static_cast<std::uint32_t>(std::max(i, j));
    pairs.emplace_back(a, b);
  }

  //
  // Sort pairs lexicographically, whichi is needed for correct bridge extension
  // The bridge logic expects pairs in order of increasing i
  //
  std::sort(pairs.begin(), pairs.end());

  return pairs;
}

inline void calculate_beta_sheets(std::vector<DsspResidue> &residues,
                                  const std::vector<std::pair<std::uint32_t, std::uint32_t>> &pairs) {
  std::vector<Bridge> bridges;
  if (residues.empty() || pairs.empty()) return;

  struct DiagonalState {
    std::size_t bridge_index = 0;
    std::uint32_t last_i     = 0;
    std::uint32_t last_j     = 0;
    bool valid               = false;
  };

  const std::size_t residue_count = residues.size();
  std::vector<DiagonalState> parallel_diags(residue_count);
  std::vector<DiagonalState> antipar_diags(residue_count * 2);
  bridges.reserve(pairs.size());

  for (const auto &[i, j] : pairs) {
    if (i >= residue_count || j >= residue_count) continue;
    const auto type = test_bridge(residues, i, j);
    if (type == BridgeType::None) continue;

    if (type == BridgeType::Parallel) {
      const auto diag = static_cast<std::size_t>(j - i);
      auto &state     = parallel_diags[diag];

      if (state.valid && i == state.last_i + 1 && j == state.last_j + 1) {
        auto &bridge = bridges[state.bridge_index];
        bridge.i.push_back(i);
        bridge.j.push_back(j);
      } else {
        Bridge bridge;
        bridge.type    = type;
        bridge.chain_i = residues[i].asym_id;
        bridge.chain_j = residues[j].asym_id;
        bridge.i.push_back(i);
        bridge.j.push_back(j);
        bridges.push_back(std::move(bridge));
        state.bridge_index = bridges.size() - 1;
      }

      state.last_i = i;
      state.last_j = j;
      state.valid  = true;
    } else {
      const auto diag = static_cast<std::size_t>(i + j);
      auto &state     = antipar_diags[diag];

      if (state.valid && i == state.last_i + 1 && j + 1 == state.last_j) {
        auto &bridge = bridges[state.bridge_index];
        bridge.i.push_back(i);
        bridge.j.insert(bridge.j.begin(), j); // insert at front for antiparallel
      } else {
        Bridge bridge;
        bridge.type    = type;
        bridge.chain_i = residues[i].asym_id;
        bridge.chain_j = residues[j].asym_id;
        bridge.i.push_back(i);
        bridge.j.push_back(j);
        bridges.push_back(std::move(bridge));
        state.bridge_index = bridges.size() - 1;
      }

      state.last_i = i;
      state.last_j = j;
      state.valid  = true;
    }
  }

  if (bridges.empty()) return;

  std::sort(bridges.begin(), bridges.end());

  // Bulge merging
  for (std::size_t i = 0; i < bridges.size(); ++i) {
    if (bridges[i].removed) continue;

    for (std::size_t j = i + 1; j < bridges.size(); ++j) {
      if (bridges[j].removed) continue;

      const auto ibi = bridges[i].i.front();
      const auto iei = bridges[i].i.back();
      const auto jbi = bridges[i].j.front();
      const auto jei = bridges[i].j.back();

      const auto ibj = bridges[j].i.front();
      if (ibj >= iei + 6) break;

      const auto iej = bridges[j].i.back();
      const auto jbj = bridges[j].j.front();
      const auto jej = bridges[j].j.back();

      // Use unsigned arithmetic to match libdssp behavior
      if (bridges[i].type != bridges[j].type ||
          !detail::no_chain_break(residues, std::min(ibi, ibj), std::max(iei, iej)) ||
          !detail::no_chain_break(residues, std::min(jbi, jbj), std::max(jei, jej)) || ibj - iei >= 6 ||
          (iei >= ibj && ibi <= iej)) {
        continue;
      }

      bool bulge = false;
      if (bridges[i].type == BridgeType::Parallel) {
        bulge = ((jbj - jei < 6 && ibj - iei < 3) || (jbj - jei < 3));
      } else {
        bulge = ((jbi - jej < 6 && ibj - iei < 3) || (jbi - jej < 3));
      }

      if (bulge) {
        bridges[i].i.insert(bridges[i].i.end(), bridges[j].i.begin(), bridges[j].i.end());
        if (bridges[i].type == BridgeType::Parallel) {
          bridges[i].j.insert(bridges[i].j.end(), bridges[j].j.begin(), bridges[j].j.end());
        } else {
          bridges[i].j.insert(bridges[i].j.begin(), bridges[j].j.begin(), bridges[j].j.end());
        }
        bridges[j].removed = true;
      }
    }
  }

  bridges.erase(
      std::remove_if(bridges.begin(), bridges.end(), [](const Bridge &b) noexcept { return b.removed; }),
      bridges.end());
  if (bridges.empty()) return;

  struct DisjointSet {
    std::vector<std::size_t> parent;
    std::vector<std::uint8_t> rank;

    explicit DisjointSet(std::size_t n) : parent(n), rank(n, 0) {
      for (std::size_t i = 0; i < n; ++i)
        parent[i] = i;
    }

    [[nodiscard]] std::size_t find(std::size_t x) noexcept {
      while (parent[x] != x) {
        parent[x] = parent[parent[x]];
        x         = parent[x];
      }
      return x;
    }

    void unite(std::size_t a, std::size_t b) noexcept {
      a = find(a);
      b = find(b);
      if (a == b) return;
      if (rank[a] < rank[b]) std::swap(a, b);
      parent[b] = a;
      if (rank[a] == rank[b]) ++rank[a];
    }
  };

  std::vector<std::vector<std::size_t>> res_to_bridge(residue_count);
  for (std::size_t idx = 0; idx < bridges.size(); ++idx) {
    for (const auto res_idx : bridges[idx].i) {
      res_to_bridge[static_cast<std::size_t>(res_idx)].push_back(idx);
    }
    for (const auto res_idx : bridges[idx].j) {
      res_to_bridge[static_cast<std::size_t>(res_idx)].push_back(idx);
    }
  }

  DisjointSet dsu(bridges.size());
  for (const auto &bridge_list : res_to_bridge) {
    if (bridge_list.size() < 2) continue;
    const auto first = bridge_list.front();
    for (std::size_t k = 1; k < bridge_list.size(); ++k) {
      dsu.unite(first, bridge_list[k]);
    }
  }

  std::vector<std::vector<std::size_t>> components(bridges.size());
  for (std::size_t idx = 0; idx < bridges.size(); ++idx) {
    components[dsu.find(idx)].push_back(idx);
  }

  struct ComponentInfo {
    std::size_t root      = 0;
    std::uint32_t min_res = 0;
  };

  std::vector<ComponentInfo> ordered;
  ordered.reserve(bridges.size());
  for (std::size_t root = 0; root < components.size(); ++root) {
    if (components[root].empty()) continue;
    std::uint32_t min_res = std::numeric_limits<std::uint32_t>::max();
    for (const auto bridge_idx : components[root]) {
      const auto &bridge = bridges[bridge_idx];
      if (!bridge.i.empty()) min_res = std::min(min_res, bridge.i.front());
      if (!bridge.j.empty()) min_res = std::min(min_res, bridge.j.front());
    }
    ordered.push_back(ComponentInfo{root, min_res});
  }

  std::sort(ordered.begin(), ordered.end(), [](const ComponentInfo &a, const ComponentInfo &b) noexcept {
    if (a.min_res != b.min_res) return a.min_res < b.min_res;
    return a.root < b.root;
  });

  std::uint32_t sheet  = 1;
  std::uint32_t ladder = 0;
  for (const auto &entry : ordered) {
    const auto &sheet_indices = components[entry.root];
    for (const auto idx : sheet_indices) {
      bridges[idx].ladder       = ladder++;
      bridges[idx].sheet        = sheet;
      bridges[idx].link_indices = sheet_indices;
    }
    ++sheet;
  }

  // Apply assignments to residues
  for (const auto &bridge : bridges) {
    const auto ss = (bridge.i.size() > 1) ? DSSPAssignment::Strand : DSSPAssignment::BetaBridge;

    for (std::uint32_t idx = bridge.i.front(); idx <= bridge.i.back(); ++idx) {
      if (residues[idx].assignment != DSSPAssignment::Strand) {
        residues[idx].assignment = ss;
      }
      residues[idx].sheet = static_cast<int>(bridge.sheet);
    }

    for (std::uint32_t idx = bridge.j.front(); idx <= bridge.j.back(); ++idx) {
      if (residues[idx].assignment != DSSPAssignment::Strand) {
        residues[idx].assignment = ss;
      }
      residues[idx].sheet = static_cast<int>(bridge.sheet);
    }
  }
}

[[nodiscard]] inline DsspHelixPosition get_helix_flag(const DsspResidue &res, DsspHelixType type) noexcept {
  return res.helix_flags[static_cast<std::size_t>(type)];
}

inline void set_helix_flag(DsspResidue &res, DsspHelixType type, DsspHelixPosition pos) noexcept {
  res.helix_flags[static_cast<std::size_t>(type)] = pos;
}

[[nodiscard]] inline bool is_helix_start(const DsspResidue &res, DsspHelixType type) noexcept {
  const auto pos = get_helix_flag(res, type);
  return pos == DsspHelixPosition::Start || pos == DsspHelixPosition::StartAndEnd;
}

inline void calculate_helices(std::vector<DsspResidue> &residues, bool prefer_pi_helices = true) {
  if (residues.empty()) return;

  for (auto &res : residues) {
    res.is_bend = !std::isnan(res.kappa) && res.kappa > 70.0;
    res.helix_flags.fill(DsspHelixPosition::None);
  }

  for (DsspHelixType helixType : {DsspHelixType::Helix3_10, DsspHelixType::Alpha, DsspHelixType::Pi}) {
    const std::uint32_t stride = static_cast<std::uint32_t>(helixType) + 3;
    for (std::uint32_t i = 0; i + stride < residues.size(); ++i) {
      if (detail::no_chain_break(residues, i, i + stride) && test_bond(residues[i + stride], residues[i])) {
        set_helix_flag(residues[i + stride], helixType, DsspHelixPosition::End);
        for (std::uint32_t j = i + 1; j < i + stride; ++j) {
          if (get_helix_flag(residues[j], helixType) == DsspHelixPosition::None) {
            set_helix_flag(residues[j], helixType, DsspHelixPosition::Middle);
          }
        }
        if (get_helix_flag(residues[i], helixType) == DsspHelixPosition::End) {
          set_helix_flag(residues[i], helixType, DsspHelixPosition::StartAndEnd);
        } else {
          set_helix_flag(residues[i], helixType, DsspHelixPosition::Start);
        }
      }
    }
  }

  for (std::uint32_t i = 1; i + 4 < residues.size(); ++i) {
    if (is_helix_start(residues[i], DsspHelixType::Alpha) &&
        is_helix_start(residues[i - 1], DsspHelixType::Alpha)) {
      for (std::uint32_t j = i; j <= i + 3; ++j)
        residues[j].assignment = DSSPAssignment::AlphaHelix;
    }
  }

  for (std::uint32_t i = 1; i + 3 < residues.size(); ++i) {
    if (is_helix_start(residues[i], DsspHelixType::Helix3_10) &&
        is_helix_start(residues[i - 1], DsspHelixType::Helix3_10)) {
      bool empty = true;
      for (std::uint32_t j = i; empty && j <= i + 2; ++j) {
        empty = residues[j].assignment == DSSPAssignment::Coil ||
                residues[j].assignment == DSSPAssignment::Helix3_10;
      }
      if (empty) {
        for (std::uint32_t j = i; j <= i + 2; ++j)
          residues[j].assignment = DSSPAssignment::Helix3_10;
      }
    }
  }

  for (std::uint32_t i = 1; i + 5 < residues.size(); ++i) {
    if (is_helix_start(residues[i], DsspHelixType::Pi) &&
        is_helix_start(residues[i - 1], DsspHelixType::Pi)) {
      bool empty = true;
      for (std::uint32_t j = i; empty && j <= i + 4; ++j) {
        empty = residues[j].assignment == DSSPAssignment::Coil ||
                residues[j].assignment == DSSPAssignment::HelixPi ||
                (prefer_pi_helices && residues[j].assignment == DSSPAssignment::AlphaHelix);
      }
      if (empty) {
        for (std::uint32_t j = i; j <= i + 4; ++j)
          residues[j].assignment = DSSPAssignment::HelixPi;
      }
    }
  }

  for (std::uint32_t i = 1; i + 1 < residues.size(); ++i) {
    if (residues[i].assignment == DSSPAssignment::Coil) {
      bool is_turn = false;
      for (DsspHelixType helixType : {DsspHelixType::Helix3_10, DsspHelixType::Alpha, DsspHelixType::Pi}) {
        const std::uint32_t stride = 3 + static_cast<std::uint32_t>(helixType);
        for (std::uint32_t k = 1; k < stride && !is_turn; ++k) {
          if (i >= k) is_turn = is_helix_start(residues[i - k], helixType);
        }
      }
      if (is_turn) {
        residues[i].assignment = DSSPAssignment::Turn;
      } else if (residues[i].is_bend) {
        residues[i].assignment = DSSPAssignment::Bend;
      }
    }
  }
}

inline void calculate_pp_helices(std::vector<DsspResidue> &residues, int stretch_length) {
  if (residues.size() < 4) return;

  const double epsilon = 29.0;
  const double phi_min = -75.0 - epsilon;
  const double phi_max = -75.0 + epsilon;
  const double psi_min = 145.0 - epsilon;
  const double psi_max = 145.0 + epsilon;

  std::vector<double> phi(residues.size(), 360.0);
  std::vector<double> psi(residues.size(), 360.0);
  for (std::size_t i = 1; i + 1 < residues.size(); ++i) {
    phi[i] = std::isfinite(residues[i].phi) ? residues[i].phi : 360.0;
    psi[i] = std::isfinite(residues[i].psi) ? residues[i].psi : 360.0;
  }

  for (std::size_t i = 1; i + 3 < residues.size(); ++i) {
    switch (stretch_length) {
      case 2: {
        if (phi[i + 0] < phi_min || phi[i + 0] > phi_max || phi[i + 1] < phi_min || phi[i + 1] > phi_max)
          continue;
        if (psi[i + 0] < psi_min || psi[i + 0] > psi_max || psi[i + 1] < psi_min || psi[i + 1] > psi_max)
          continue;

        switch (get_helix_flag(residues[i], DsspHelixType::Ppii)) {
          case DsspHelixPosition::None:
            set_helix_flag(residues[i], DsspHelixType::Ppii, DsspHelixPosition::Start);
            break;
          case DsspHelixPosition::End:
            set_helix_flag(residues[i], DsspHelixType::Ppii, DsspHelixPosition::Middle);
            break;
          default:
            break;
        }

        set_helix_flag(residues[i + 1], DsspHelixType::Ppii, DsspHelixPosition::End);
        if (residues[i].assignment == DSSPAssignment::Coil)
          residues[i].assignment = DSSPAssignment::PolyProlineHelix;
        if (residues[i + 1].assignment == DSSPAssignment::Coil)
          residues[i + 1].assignment = DSSPAssignment::PolyProlineHelix;
      } break;
      case 3: {
        if (phi[i + 0] < phi_min || phi[i + 0] > phi_max || phi[i + 1] < phi_min || phi[i + 1] > phi_max ||
            phi[i + 2] < phi_min || phi[i + 2] > phi_max)
          continue;
        if (psi[i + 0] < psi_min || psi[i + 0] > psi_max || psi[i + 1] < psi_min || psi[i + 1] > psi_max ||
            psi[i + 2] < psi_min || psi[i + 2] > psi_max)
          continue;

        switch (get_helix_flag(residues[i], DsspHelixType::Ppii)) {
          case DsspHelixPosition::None:
            set_helix_flag(residues[i], DsspHelixType::Ppii, DsspHelixPosition::Start);
            break;
          case DsspHelixPosition::End:
            set_helix_flag(residues[i], DsspHelixType::Ppii, DsspHelixPosition::StartAndEnd);
            break;
          default:
            break;
        }

        set_helix_flag(residues[i + 1], DsspHelixType::Ppii, DsspHelixPosition::Middle);
        set_helix_flag(residues[i + 2], DsspHelixType::Ppii, DsspHelixPosition::End);

        if (residues[i + 0].assignment == DSSPAssignment::Coil)
          residues[i + 0].assignment = DSSPAssignment::PolyProlineHelix;
        if (residues[i + 1].assignment == DSSPAssignment::Coil)
          residues[i + 1].assignment = DSSPAssignment::PolyProlineHelix;
        if (residues[i + 2].assignment == DSSPAssignment::Coil)
          residues[i + 2].assignment = DSSPAssignment::PolyProlineHelix;
      } break;
      default:
        return;
    }
  }
}

inline void assign_secondary_structure(std::vector<DsspResidue> &residues,
                                       const std::vector<std::pair<std::uint32_t, std::uint32_t>> &pairs,
                                       int pp_stretch_length = 2, bool prefer_pi_helices = true) {
  if (residues.empty()) return;

  for (auto &res : residues) {
    res.assignment = DSSPAssignment::Coil;
    res.sheet      = 0;
  }

  calculate_beta_sheets(residues, pairs);
  calculate_helices(residues, prefer_pi_helices);
  calculate_pp_helices(residues, pp_stretch_length);
}

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_DSSP_SECONDARY_HPP
