/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto make = []() -> decltype(auto) {
 *     return std::string{"besian"} + "sejdiu" + "@gmail.com";
 *   };
 *   return make();
 * }();
 *
 */

#ifndef LAHUTA_ANALYSIS_DSSP_PRECOMPUTE_HPP
#define LAHUTA_ANALYSIS_DSSP_PRECOMPUTE_HPP

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

#include <Geometry/point.h>
#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/Conformer.h>
#include <rdkit/GraphMol/MonomerInfo.h>

#include "models/dssp.hpp"
#include "models/tables.hpp"
#include "residues/residues.hpp"
#include "topology.hpp"
#include "utils/math_constants.hpp"
#include "utils/span.hpp"

namespace lahuta::analysis {

inline constexpr double DsspMaxPeptideBondLength = 2.5;

struct DsspAtom {
  RDGeom::Point3D pos;
};

enum class DsspChainBreak : std::uint8_t { None = 0, Gap = 1, NewChain = 2 };
enum class DsspHelixType : std::uint8_t { Helix3_10 = 0, Alpha = 1, Pi = 2, Ppii = 3, Count = 4 };
enum class DsspHelixPosition : std::uint8_t { None = 0, Start = 1, Middle = 2, End = 3, StartAndEnd = 4 };

struct DsspResidue {
  std::string asym_id;
  int seq_id        = 0;
  int number        = 0; // internal numbering with gaps across breaks
  std::size_t index = 0;

  char aa = '?';
  std::string name;

  DSSPAssignment assignment = DSSPAssignment::Coil;

  bool complete               = false;
  bool is_proline             = false;
  bool has_oxt                = false;
  DsspChainBreak break_before = DsspChainBreak::None;
  bool is_bend                = false;
  int sheet                   = 0;

  RDGeom::Point3D n;
  RDGeom::Point3D ca;
  RDGeom::Point3D c;
  RDGeom::Point3D o;
  RDGeom::Point3D h;
  RDGeom::Point3D oxt;

  std::vector<DsspAtom> side_chain;

  std::array<DsspHelixPosition, static_cast<std::size_t>(DsspHelixType::Count)> helix_flags{};

  struct HBond {
    int res       = -1;
    double energy = 0.0;
  };

  HBond hbond_acceptor[2] = {};
  HBond hbond_donor[2]    = {};

  int prev = -1;
  int next = -1;

  double phi   = std::numeric_limits<double>::quiet_NaN();
  double psi   = std::numeric_limits<double>::quiet_NaN();
  double omega = std::numeric_limits<double>::quiet_NaN();
  double tco   = std::numeric_limits<double>::quiet_NaN();
  double alpha = std::numeric_limits<double>::quiet_NaN();
  double kappa = std::numeric_limits<double>::quiet_NaN();
};

namespace detail {

[[nodiscard]] inline double distance_sq(const RDGeom::Point3D &a, const RDGeom::Point3D &b) noexcept {
  const double dx = a.x - b.x;
  const double dy = a.y - b.y;
  const double dz = a.z - b.z;
  return dx * dx + dy * dy + dz * dz;
}

[[nodiscard]] inline double distance(const RDGeom::Point3D &a, const RDGeom::Point3D &b) noexcept {
  return std::sqrt(distance_sq(a, b));
}

[[nodiscard]] inline double dihedral_angle(const RDGeom::Point3D &p1, const RDGeom::Point3D &p2,
                                           const RDGeom::Point3D &p3, const RDGeom::Point3D &p4) noexcept {
  const RDGeom::Point3D v12 = p1 - p2;
  const RDGeom::Point3D v43 = p4 - p3;
  const RDGeom::Point3D z   = p2 - p3;
  const RDGeom::Point3D p   = z.crossProduct(v12);
  const RDGeom::Point3D x   = z.crossProduct(v43);
  const RDGeom::Point3D y   = z.crossProduct(x);

  double u = x.dotProduct(x);
  double v = y.dotProduct(y);

  double result = 360.0;
  if (u > 0.0 && v > 0.0) {
    u = p.dotProduct(x) / std::sqrt(u);
    v = p.dotProduct(y) / std::sqrt(v);
    if (u != 0.0 || v != 0.0) {
      result = std::atan2(v, u) * (180.0 / lahuta::Pi);
    }
  }
  return result;
}

[[nodiscard]] inline double cosinus_angle(const RDGeom::Point3D &p1, const RDGeom::Point3D &p2,
                                          const RDGeom::Point3D &p3, const RDGeom::Point3D &p4) noexcept {
  const RDGeom::Point3D v12 = p1 - p2;
  const RDGeom::Point3D v34 = p3 - p4;
  const double denom        = v12.dotProduct(v12) * v34.dotProduct(v34);
  if (denom > 0.0) {
    return v12.dotProduct(v34) / std::sqrt(denom);
  }
  return 0.0;
}

[[nodiscard]] inline bool atom_name_equals(std::string_view name, std::string_view target) noexcept {
  std::size_t ti = 0;
  for (char c : name) {
    if (c == ' ' || c == '\t' || c == '\r' || c == '\n') continue;
    if (ti >= target.size() || c != target[ti++]) return false;
  }
  return ti == target.size();
}

[[nodiscard]] inline constexpr bool is_backbone_atom(std::string_view name) noexcept {
  return name == "N" || name == "CA" || name == "C" || name == "O";
}

[[nodiscard]] inline bool is_backbone_atom_raw(std::string_view name) noexcept {
  return atom_name_equals(name, "N") || atom_name_equals(name, "CA") || atom_name_equals(name, "C") ||
         atom_name_equals(name, "O");
}

template <typename PointT>
[[nodiscard]] inline RDGeom::Point3D to_point3d(const PointT &pt) noexcept {
  return RDGeom::Point3D(static_cast<double>(pt.x), static_cast<double>(pt.y), static_cast<double>(pt.z));
}

[[nodiscard]] inline char aa_from_resname(std::string_view resname) noexcept {
  for (char aa : std::string_view("ACDEFGHIKLMNPQRSTVWY")) {
    const auto &entry = StandardAminoAcidDataTable[aa];
    if (entry.name && resname == entry.name) return aa;
  }
  return '?';
}

inline void assign_hydrogen(DsspResidue &cur, const DsspResidue *prev) noexcept {
  cur.h = cur.n;
  if (!cur.is_proline && prev) {
    const double co_dist = distance(prev->c, prev->o);
    if (co_dist > 0.0) {
      cur.h.x += (prev->c.x - prev->o.x) / co_dist;
      cur.h.y += (prev->c.y - prev->o.y) / co_dist;
      cur.h.z += (prev->c.z - prev->o.z) / co_dist;
    }
  }
}

[[nodiscard]] inline bool no_chain_break(const std::vector<DsspResidue> &residues, std::size_t a,
                                         std::size_t b) noexcept {
  if (a >= residues.size() || b >= residues.size()) return false;
  if (a == b) return true;
  if (residues[a].asym_id != residues[b].asym_id) return false;
  std::size_t i = a;
  while (i != b) {
    const int next = residues[i].next;
    if (next < 0 || static_cast<std::size_t>(next) >= residues.size()) return false;
    if (residues[next].break_before != DsspChainBreak::None) return false;
    i = static_cast<std::size_t>(next);
  }
  return true;
}

template <typename PointT>
[[nodiscard]] inline bool build_residues_from_payload_impl(std::string_view sequence,
                                                           span<const PointT> coords,
                                                           std::vector<DsspResidue> &out, std::string &error,
                                                           std::string_view chain_id) {
  out.clear();
  if (sequence.empty()) {
    error = "Missing sequence";
    return false;
  }
  if (coords.empty()) {
    error = "Missing coordinates";
    return false;
  }

  std::size_t coord_index = 0;
  out.reserve(sequence.size());

  for (std::size_t res_idx = 0; res_idx < sequence.size(); ++res_idx) {
    const char aa = sequence[res_idx];
    if (!StandardAminoAcidDataTable.is_valid(aa)) {
      error = "Unsupported residue '" + std::string(1, aa) + "' at position " + std::to_string(res_idx);
      return false;
    }
    const auto &entry = StandardAminoAcidDataTable[aa];
    if (!entry.name) {
      error = "Missing residue name for '" + std::string(1, aa) + "'";
      return false;
    }

    DsspResidue res;
    res.asym_id    = std::string(chain_id);
    res.seq_id     = static_cast<int>(res_idx + 1);
    res.aa         = aa;
    res.name       = entry.name;
    res.is_proline = (aa == 'P');

    bool has_n  = false;
    bool has_ca = false;
    bool has_c  = false;
    bool has_o  = false;

    for (std::size_t atom_idx = 0; atom_idx < entry.size; ++atom_idx) {
      if (coord_index >= coords.size()) {
        error = "Coordinate count mismatch (payload too short)";
        return false;
      }
      const auto pos        = to_point3d(coords[coord_index++]);
      const char *atom_name = entry.atoms[atom_idx];
      if (!atom_name || atom_name[0] == '\0') continue;

      const std::string_view name_sv(atom_name);
      if (name_sv == "N") {
        res.n = pos;
        has_n = true;
      } else if (name_sv == "CA") {
        res.ca = pos;
        has_ca = true;
      } else if (name_sv == "C") {
        res.c = pos;
        has_c = true;
      } else if (name_sv == "O") {
        res.o = pos;
        has_o = true;
      } else if (name_sv == "OXT") {
        res.has_oxt = true;
        res.oxt     = pos;
      } else {
        res.side_chain.emplace_back(DsspAtom{pos});
      }
    }

    res.complete = has_n && has_ca && has_c && has_o;
    if (res.complete) {
      out.push_back(std::move(res));
    }
  }

  if (coord_index + 1 == coords.size() && !out.empty()) {
    out.back().has_oxt = true;
    out.back().oxt     = to_point3d(coords[coord_index++]);
  }

  if (coord_index != coords.size()) {
    error = "Coordinate count mismatch (payload leftover)";
    return false;
  }
  return true;
}

} // namespace detail

inline void link_residues(std::vector<DsspResidue> &residues) {
  if (residues.empty()) return;
  int number = 1;
  for (std::size_t i = 0; i < residues.size(); ++i) {
    auto &cur        = residues[i];
    cur.index        = i;
    cur.prev         = i == 0 ? -1 : static_cast<int>(i - 1);
    cur.next         = (i + 1 < residues.size()) ? static_cast<int>(i + 1) : -1;
    cur.break_before = DsspChainBreak::None;

    if (i == 0) {
      cur.number = number;
      continue;
    }

    const auto &prev        = residues[i - 1];
    const bool chain_change = prev.asym_id != cur.asym_id;
    const bool seq_gap      = (prev.seq_id + 1) != cur.seq_id;
    const bool bond_gap     = detail::distance(prev.c, cur.n) > DsspMaxPeptideBondLength;
    if (chain_change || seq_gap || bond_gap) {
      cur.break_before  = chain_change ? DsspChainBreak::NewChain : DsspChainBreak::Gap;
      number           += 2; // create a numbering gap to model chain discontinuity
    } else {
      number += 1;
    }
    cur.number = number;
  }
}

inline void assign_hydrogens(std::vector<DsspResidue> &residues) {
  for (std::size_t i = 0; i < residues.size(); ++i) {
    DsspResidue *prev = (i > 0) ? &residues[i - 1] : nullptr;
    detail::assign_hydrogen(residues[i], prev);
  }
}

inline void compute_backbone_geometry(std::vector<DsspResidue> &residues) {
  if (residues.empty()) return;
  link_residues(residues);
  assign_hydrogens(residues);

  const std::size_t n = residues.size();
  for (std::size_t i = 0; i < n; ++i) {
    auto &cur = residues[i];

    if (cur.next >= 0) {
      const std::size_t j = static_cast<std::size_t>(cur.next);
      if (detail::no_chain_break(residues, i, j)) {
        const auto &next = residues[j];
        cur.psi          = detail::dihedral_angle(cur.n, cur.ca, cur.c, next.n);
        cur.omega        = detail::dihedral_angle(cur.ca, cur.c, next.n, next.ca);
      }
    }

    if (cur.prev >= 0) {
      const std::size_t j = static_cast<std::size_t>(cur.prev);
      if (detail::no_chain_break(residues, j, i)) {
        const auto &prev = residues[j];
        cur.phi          = detail::dihedral_angle(prev.c, cur.n, cur.ca, cur.c);
        cur.tco          = detail::cosinus_angle(cur.c, cur.o, prev.c, prev.o);
      }
    }

    if (i >= 1 && i + 2 < n) {
      const std::size_t a = i - 1;
      const std::size_t b = i + 2;
      if (detail::no_chain_break(residues, a, b) && residues[a].seq_id + 3 == residues[b].seq_id) {
        cur.alpha = detail::dihedral_angle(residues[a].ca, cur.ca, residues[i + 1].ca, residues[b].ca);
      }
    }

    if (i >= 2 && i + 2 < n) {
      const std::size_t a = i - 2;
      const std::size_t b = i + 2;
      if (detail::no_chain_break(residues, a, b) && residues[a].seq_id + 4 == residues[b].seq_id) {
        const double ckap  = detail::cosinus_angle(cur.ca, residues[a].ca, residues[b].ca, cur.ca);
        const double skap  = std::sqrt(1.0 - ckap * ckap);
        const double kappa = std::atan2(skap, ckap) * (180.0 / lahuta::Pi);
        if (!std::isnan(kappa)) cur.kappa = kappa;
      }
    }
  }
}

[[nodiscard]] inline bool build_residues_from_payload(std::string_view sequence,
                                                      span<const RDGeom::Point3D> coords,
                                                      std::vector<DsspResidue> &out, std::string &error,
                                                      std::string_view chain_id = "A") {
  return detail::build_residues_from_payload_impl(sequence, coords, out, error, chain_id);
}

[[nodiscard]] inline bool build_residues_from_payload(std::string_view sequence,
                                                      span<const RDGeom::Point3Df> coords,
                                                      std::vector<DsspResidue> &out, std::string &error,
                                                      std::string_view chain_id = "A") {
  return detail::build_residues_from_payload_impl(sequence, coords, out, error, chain_id);
}

[[nodiscard]] inline bool build_residues_from_topology(const Topology &topology, const RDKit::Conformer &conf,
                                                       std::vector<DsspResidue> &out, std::string &error) {
  out.clear();
  const auto &residues = topology.get_residues().get_residues();
  if (residues.empty()) {
    error = "Missing residues in topology";
    return false;
  }

  out.reserve(residues.size());
  for (const auto &res : residues) {
    if (res.atoms.empty()) continue;
    const char aa = detail::aa_from_resname(res.name);
    if (aa == '?') continue; // skip nonstandard residues

    DsspResidue cur;
    cur.asym_id    = res.chain_id;
    cur.seq_id     = res.number;
    cur.aa         = aa;
    cur.name       = res.name;
    cur.is_proline = (aa == 'P') || res.name == "PRO";

    bool has_n  = false;
    bool has_ca = false;
    bool has_c  = false;
    bool has_o  = false;

    for (const auto *atom : res.atoms) {
      if (!atom) continue;
      if (atom->getAtomicNum() == 1) continue; // skip hydrogens

      auto *info = dynamic_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
      if (!info) continue;
      const std::string_view raw_name = info->getName();
      const auto pos                  = conf.getAtomPos(atom->getIdx());

      if (detail::atom_name_equals(raw_name, "N")) {
        cur.n = pos;
        has_n = true;
      } else if (detail::atom_name_equals(raw_name, "CA")) {
        cur.ca = pos;
        has_ca = true;
      } else if (detail::atom_name_equals(raw_name, "C")) {
        cur.c = pos;
        has_c = true;
      } else if (detail::atom_name_equals(raw_name, "O")) {
        cur.o = pos;
        has_o = true;
      } else if (detail::atom_name_equals(raw_name, "OXT")) {
        cur.has_oxt = true;
        cur.oxt     = pos;
      } else if (!detail::is_backbone_atom_raw(raw_name)) {
        cur.side_chain.emplace_back(DsspAtom{pos});
      }
    }

    cur.complete = has_n && has_ca && has_c && has_o;
    if (cur.complete) {
      out.push_back(std::move(cur));
    }
  }

  if (out.empty()) {
    error = "No complete residues found in topology";
    return false;
  }
  return true;
}

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_DSSP_PRECOMPUTE_HPP
