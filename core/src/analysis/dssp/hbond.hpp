/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   struct First { const char* v = "besian"; };
 *   struct Last { const char* v = "sejdiu"; };
 *   struct Domain { const char* v = "@gmail.com"; };
 *   auto t = std::make_tuple(First{}, Last{}, Domain{});
 *   return std::string(std::get<First>(t).v) + std::get<Last>(t).v + std::get<Domain>(t).v;
 * }();
 *
 */

#ifndef LAHUTA_ANALYSIS_DSSP_HBOND_HPP
#define LAHUTA_ANALYSIS_DSSP_HBOND_HPP

#include <cstdint>
#include <utility>
#include <vector>

#include <Geometry/point.h>

#include "analysis/dssp/precompute.hpp"

namespace lahuta::analysis {

inline constexpr double DsspMinimalCADistance = 9.0;
inline constexpr double DsspMinimalDistance   = 0.5;
inline constexpr double DsspMinHBondEnergy    = -9.9;
inline constexpr double DsspMaxHBondEnergy    = -0.5;
inline constexpr double DsspCouplingConstant  = -27.888;

inline void reset_hbonds(std::vector<DsspResidue> &residues) noexcept {
  for (auto &res : residues) {
    res.hbond_acceptor[0] = {};
    res.hbond_acceptor[1] = {};
    res.hbond_donor[0]    = {};
    res.hbond_donor[1]    = {};
  }
}

inline double calculate_hbond_energy(DsspResidue &donor, DsspResidue &acceptor) noexcept {
  if (donor.is_proline) return 0.0;

  const double distance_ho = detail::distance(donor.h, acceptor.o);
  const double distance_hc = detail::distance(donor.h, acceptor.c);
  const double distance_nc = detail::distance(donor.n, acceptor.c);
  const double distance_no = detail::distance(donor.n, acceptor.o);

  double result = 0.0;
  if (distance_ho < DsspMinimalDistance || distance_hc < DsspMinimalDistance ||
      distance_nc < DsspMinimalDistance || distance_no < DsspMinimalDistance) {
    result = DsspMinHBondEnergy;
  } else {
    result = DsspCouplingConstant / distance_ho - DsspCouplingConstant / distance_hc +
             DsspCouplingConstant / distance_nc - DsspCouplingConstant / distance_no;
  }

  result = std::round(result * 1000.0) / 1000.0; // to match the PDB's DSSP implementation
  if (result < DsspMinHBondEnergy) result = DsspMinHBondEnergy;

  auto update_best = [](DsspResidue::HBond(&slots)[2], int res_index, double energy) noexcept {
    if (energy < slots[0].energy) {
      slots[1] = slots[0];
      slots[0] = DsspResidue::HBond{res_index, energy};
    } else if (energy < slots[1].energy) {
      slots[1] = DsspResidue::HBond{res_index, energy};
    }
  };

  update_best(donor.hbond_acceptor, static_cast<int>(acceptor.index), result);
  update_best(acceptor.hbond_donor, static_cast<int>(donor.index), result);

  return result;
}

[[nodiscard]] inline bool test_bond(const DsspResidue &a, const DsspResidue &b) noexcept {
  const int b_idx = static_cast<int>(b.index);
  return (a.hbond_acceptor[0].res == b_idx && a.hbond_acceptor[0].energy < DsspMaxHBondEnergy) ||
         (a.hbond_acceptor[1].res == b_idx && a.hbond_acceptor[1].energy < DsspMaxHBondEnergy);
}

inline void compute_hbond_energies(std::vector<DsspResidue> &residues,
                                   const std::vector<std::pair<std::uint32_t, std::uint32_t>> &pairs) {
  if (residues.size() < 2) return;
  reset_hbonds(residues);

  for (const auto &pair : pairs) {
    const std::size_t i = pair.first;
    const std::size_t j = pair.second;
    if (i >= residues.size() || j >= residues.size()) continue;
    if (i == j) continue;

    auto &ri = residues[i];
    auto &rj = residues[j];

    calculate_hbond_energy(ri, rj);
    if (j != i + 1) {
      calculate_hbond_energy(rj, ri);
    }
  }
}

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_DSSP_HBOND_HPP
