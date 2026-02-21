/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   constexpr char p1[] = "besian", p2[] = "sejdiu", p3[] = "@gmail.com"; std::string s;
 *   s.append(std::begin(p1), std::end(p1) - 1);
 *   s.append(std::begin(p2), std::end(p2) - 1);
 *   s.append(std::begin(p3), std::end(p3) - 1);
 *   return s;
 * }();
 *
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/AtomIterators.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/RingInfo.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/Substruct/SubstructMatch.h>

#include "bond_order.hpp"

// clang-format off
namespace lahuta {
namespace {

using BondType = RDKit::Bond::BondType;
using Hybridization = RDKit::Atom::HybridizationType;
struct HybPattern {
  const char *smarts;
  Hybridization hyb;
};

constexpr std::array<HybPattern, 17> HybPatterns{{
    {"[D4]", Hybridization::SP3},
    {"[D5]", Hybridization::SP3D},
    {"[D6]", Hybridization::SP3D2},
    {"[C]", Hybridization::SP3},
    {"[c,$(C=*)]", Hybridization::SP2},
    {"[$(C#*),$(C(=*)=*)]", Hybridization::SP},
    {"[N]", Hybridization::SP3},
    {"[n,$(N=*),$(N[#6,#7,#8]=,:,#*)]", Hybridization::SP2},
    {"[ND1,ND2,ND3]a", Hybridization::SP2},
    {"[$(N#*),$([ND2](=*)=*)]", Hybridization::SP},
    {"[O]", Hybridization::SP3},
    {"[o,$(O=*),$(O[#6,#7,#8]=,:*)]", Hybridization::SP2},
    {"[$([#8D1][#6][#8D1])]", Hybridization::SP2},
    {"[$(O#*)]", Hybridization::SP},
    {"[S]", Hybridization::SP3},
    {"[#16;s,$([SD1]=*)]", Hybridization::SP2},
    {"[SD6]", Hybridization::SP3D2},
}};

struct BondAssignmentPattern {
  const char *smarts;
  std::vector<int> bond_triplets;
};

const std::array<std::shared_ptr<RDKit::ROMol>, HybPatterns.size()> &hybridization_queries() {
  thread_local const std::array<std::shared_ptr<RDKit::ROMol>, HybPatterns.size()> queries = [] {
    std::array<std::shared_ptr<RDKit::ROMol>, HybPatterns.size()> tmp{};
    for (size_t i = 0; i < HybPatterns.size(); ++i) {
      tmp[i].reset(RDKit::SmartsToMol(HybPatterns[i].smarts));
    }
    return tmp;
  }();
  return queries;
}

const std::vector<BondAssignmentPattern> &bond_assignment_patterns() {
  static const std::vector<BondAssignmentPattern> patterns = {
      {"[x2,x3]1[#6]([#7D3]2)[#6][#6][#6]2[x2,x3][#6]([#7D3]3)[#6][#6][#6]3[x2,x3][#6]([#7D3]4)[#6][#6][#6]4["
       "x2,x3][#6]([#7D3]5)[#6][#6][#6]51",
       {0,  1,  2, 1,  2,  1, 1,  3,  1, 3,  4,  2, 4,  5,  1, 5,  2,  1, 5,  6,  2, 6,  7,  1, 7,  8,  2,
        7,  9,  1, 9,  10, 2, 10, 11, 1, 11, 8,  1, 11, 12, 2, 12, 13, 1, 13, 14, 1, 13, 15, 2, 15, 16, 1,
        16, 17, 2, 17, 14, 1, 17, 18, 1, 18, 19, 2, 19, 20, 1, 19, 21, 1, 21, 22, 2, 22, 23, 1, 23, 20, 2}},
      {"[x2,x3]1[#6]([#7D3]2)[#6][#6][#6]2[x2,x3][#6]([#7]3)[#6][#6][#6]3[x2,x3][#6]([#7D3]4)[#6][#6][#6]4["
       "x2,x3][#6]([#7]5)[#6][#6][#6]51",
       {0,  1,  2, 1,  2,  1, 1,  3,  1, 3,  4,  2, 4,  5,  1, 5,  2,  1, 5,  6,  2, 6,  7,  1, 7,  8,  2,
        7,  9,  1, 9,  10, 2, 10, 11, 1, 11, 8,  1, 11, 12, 2, 12, 13, 1, 13, 14, 1, 13, 15, 2, 15, 16, 1,
        16, 17, 2, 17, 14, 1, 17, 18, 1, 18, 19, 2, 19, 20, 1, 19, 21, 1, 21, 22, 2, 22, 23, 1, 23, 20, 2}},
      {"[x2,x3]1[#6]([#7]2)[#6][#6][#6]2[x2,x3][#6]([#7]3)[#6][#6][#6]3[x2,x3][#6]([#7]4)[#6][#6][#6]4[x2,x3]"
       "[#6]([#7]5)[#6][#6][#6]51",
       {0,  1,  2, 1,  2,  1, 1,  3,  1, 3,  4,  2, 4,  5,  1, 5,  2,  1, 5,  6,  2, 6,  7,  1, 7,  8,  2,
        7,  9,  1, 9,  10, 2, 10, 11, 1, 11, 8,  1, 11, 12, 2, 12, 13, 1, 13, 14, 1, 13, 15, 2, 15, 16, 1,
        16, 17, 2, 17, 14, 1, 17, 18, 1, 18, 19, 2, 19, 20, 1, 19, 21, 1, 21, 22, 2, 22, 23, 1, 23, 20, 2}},
      {"[#7D2][#7D2^1][#7D1]", {0, 1, 2, 1, 2, 2}},
      {"[#8D1][#7D3^2]([#8D1])*", {0, 1, 2, 1, 2, 2, 1, 3, 1}},
      {"[#16D4]([#8D1])([#8D1])([*!#8])([*!#8])", {0, 1, 2, 0, 2, 2, 0, 3, 1, 0, 4, 1}},
      {"[#16D4]([#8D1])([#8D1])([#8-,#8D1])([#8-,#8D1])", {0, 1, 2, 0, 2, 2, 0, 3, 1, 0, 4, 1}},
      {"[#16D4]([#16D1])([#8D1])([#8-,#8])([#8-,#8])", {0, 1, 2, 0, 2, 2, 0, 3, 1, 0, 4, 1}},
      {"[#16D3]([#8D1])([*!#8])([*!#8])", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
      {"[#16D3]([#8D1])([#8D1-])([#8D1-])", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
      {"[#16D3]([#8D1])([#8])([#8])", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
      {"[#16D2]([#8D1])([#16D1])", {0, 1, 2, 0, 2, 2}},
      {"[#16D2]([#8D1])([*!#8])", {0, 1, 2, 0, 2, 1}},
      {"[#16D2]([#8D1])([#8D1])", {0, 1, 2, 0, 2, 2}},
      {"[#15D3]([#8D1])([#8D1])([#8D2])", {0, 1, 2, 0, 2, 2, 0, 3, 1}},
      {"[#7D2]([#8D1])([#1])", {0, 1, 2, 0, 2, 1}},
      {"[#15D4]([#8D1])(*)(*)(*)", {0, 1, 2, 0, 2, 1, 0, 3, 1, 0, 4, 1}},
      {"[#6D3^2]([#8D1])([#8])*", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
      {"[#8D1][#6D2^1][#8D1]", {0, 1, 2, 1, 2, 2}},
      {"[#6D3^2]([#8D1;!-])([#7])*", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
      {"[#34D3^2]([#8D1])([#8])*", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
      {"[#6D3^2]([#8D1])([#16])*", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
      {"[#6D3^2]([#16D1])([#16])*", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
      {"[CD3^2]([#16D1])([N])*", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
  };
  return patterns;
}

const std::vector<std::shared_ptr<RDKit::ROMol>> &bond_assignment_queries() {
  thread_local const std::vector<std::shared_ptr<RDKit::ROMol>> queries = [] {
    std::vector<std::shared_ptr<RDKit::ROMol>> tmp;
    tmp.reserve(bond_assignment_patterns().size());
    for (const auto &pattern : bond_assignment_patterns()) {
      tmp.emplace_back(RDKit::SmartsToMol(pattern.smarts));
    }
    return tmp;
  }();
  return queries;
}

double to_degrees(double radians) {
  constexpr double RadToDeg = 180.0 / M_PI;
  return radians * RadToDeg;
}

double compute_dihedral(const RDGeom::Point3D &a, const RDGeom::Point3D &b, const RDGeom::Point3D &c, const RDGeom::Point3D &d) {
  const RDGeom::Point3D b1 = a - b;
  const RDGeom::Point3D b2 = b - c;
  const RDGeom::Point3D b3 = c - d;

  const RDGeom::Point3D b2xb3 = b2.crossProduct(b3);
  const RDGeom::Point3D b1xb2 = b1.crossProduct(b2);
  const double rb2_sq = b2.dotProduct(b2);

  const double y = std::sqrt(rb2_sq) * b1.dotProduct(b2xb3);
  const double x = b1xb2.dotProduct(b2xb3);
  return -to_degrees(std::atan2(y, x));
}

double average_ring_dihedral(const RDKit::Conformer &conf, const std::vector<int> &ring_atom_indices) {
  if (ring_atom_indices.size() < 4) return 0.0;

  double sum = 0.0;
  for (size_t i = 0; i < ring_atom_indices.size(); ++i) {
    const auto idx0 = ring_atom_indices[i];
    const auto idx1 = ring_atom_indices[(i + 1) % ring_atom_indices.size()];
    const auto idx2 = ring_atom_indices[(i + 2) % ring_atom_indices.size()];
    const auto idx3 = ring_atom_indices[(i + 3) % ring_atom_indices.size()];

    sum += std::fabs(compute_dihedral(
        conf.getAtomPos(idx0),
        conf.getAtomPos(idx1),
        conf.getAtomPos(idx2),
        conf.getAtomPos(idx3)));
  }

  return sum / static_cast<double>(ring_atom_indices.size());
}

double average_bond_angle(const RDKit::Atom &atom) {
  const auto &mol = atom.getOwningMol();
  const auto &conf = mol.getConformer();
  const auto center_idx = atom.getIdx();
  const auto center_pos = conf.getAtomPos(center_idx);

  std::vector<RDGeom::Point3D> directions;
  directions.reserve(atom.getDegree());
  for (const auto *bond : mol.atomBonds(&atom)) {
    const RDKit::Atom *neighbor = bond->getOtherAtom(&atom);
    RDGeom::Point3D vector = conf.getAtomPos(neighbor->getIdx()) - center_pos;
    const double length_sq = vector.lengthSq();
    if (length_sq > 1e-12) {
      vector /= std::sqrt(length_sq);
      directions.push_back(vector);
    }
  }

  if (directions.size() < 2) return 180.0;

  double angle_sum = 0.0;
  int count = 0;
  for (size_t i = 0; i < directions.size(); ++i) {
    for (size_t j = i + 1; j < directions.size(); ++j) {
      const double dot = std::clamp(directions[i].dotProduct(directions[j]), -1.0, 1.0);
      angle_sum += to_degrees(std::acos(dot));
      ++count;
    }
  }

  return count > 0 ? angle_sum / static_cast<double>(count) : 180.0;
}

double bond_length_sq(const RDKit::Conformer &conf, const RDKit::Atom &a, const RDKit::Atom &b) {
  return (conf.getAtomPos(a.getIdx()) - conf.getAtomPos(b.getIdx())).lengthSq();
}

bool atom_has_bond_type(const RDKit::Atom &atom, BondType type) {
  const auto &mol = atom.getOwningMol();
  for (const auto *bond : mol.atomBonds(&atom)) {
    if (bond->getBondType() == type) return true;
  }
  return false;
}

void assign_initial_hybridization(RDKit::RWMol &mol) {
  for (auto *atom : mol.atoms()) {
    atom->setHybridization(Hybridization::SP);
  }

  RDKit::SubstructMatchParameters params;
  params.maxMatches = static_cast<unsigned int>(mol.getNumAtoms());

  const auto &queries = hybridization_queries();
  for (size_t i = 0; i < HybPatterns.size(); ++i) {
    const auto matches = RDKit::SubstructMatch(mol, *queries[i], params);
    for (const auto &match : matches) {
      for (const auto &pair : match) {
        mol.getAtomWithIdx(pair.second)->setHybridization(HybPatterns[i].hyb);
      }
    }
  }
}

void refine_hybridization_by_geometry(RDKit::RWMol &mol) {
  const auto *ring_info = mol.getRingInfo();
  if (!ring_info || !ring_info->isSssrOrBetter()) {
    RDKit::MolOps::findSSSR(mol);
  }

  for (auto *atom : mol.atoms()) {
    const double avg_angle = average_bond_angle(*atom);

    if (avg_angle > 155.0) {
      atom->setHybridization(Hybridization::SP);
    } else if (avg_angle <= 155.0 && avg_angle > 115.0) {
      atom->setHybridization(Hybridization::SP2);
    }

    const bool in_ring = ring_info && ring_info->numAtomRings(atom->getIdx()) != 0;

    // Special case for imines
    if (atom->getAtomicNum() == 7 && atom->getNumExplicitHs() == 1 && atom->getDegree() == 2 && avg_angle > 109.5) {
      atom->setHybridization(Hybridization::SP2);
    } else if (atom->getAtomicNum() == 7 && atom->getDegree() == 2 && in_ring) {
      // Azete
      atom->setHybridization(Hybridization::SP2);
    }
  }
}

void refine_ring_hybridization(RDKit::RWMol &mol, RDKit::Conformer &conf) {
  const auto *ring_info = mol.getRingInfo();
  if (!ring_info) return;

  for (const auto &ring : ring_info->atomRings()) {
    if (ring.size() == 5 && average_ring_dihedral(conf, ring) <= 7.5) {
      for (int idx : ring) {
        RDKit::Atom *atom = mol.getAtomWithIdx(idx);
        if (atom->getDegree() == 2) {
          atom->setHybridization(Hybridization::SP2);
        }
      }
    } else if (ring.size() == 6 && average_ring_dihedral(conf, ring) <= 12.0) {
      for (int idx : ring) {
        RDKit::Atom *atom = mol.getAtomWithIdx(idx);
        const auto degree = atom->getDegree();
        if (degree == 2 || degree == 3) {
          atom->setHybridization(Hybridization::SP2);
        }
      }
    }
  }
}

void relax_terminal_assignments(RDKit::RWMol &mol) {
  for (auto *atom : mol.atoms()) {
    const auto hybrid = atom->getHybridization();
    if (hybrid != Hybridization::SP && hybrid != Hybridization::SP2) continue;

    bool has_open_neighbor = false;
    for (const auto *bond : mol.atomBonds(atom)) {
      const RDKit::Atom *neighbor = bond->getOtherAtom(atom);
      if (neighbor->getHybridization() < Hybridization::SP3 || neighbor->getDegree() == 1) {
        has_open_neighbor = true;
        break;
      }
    }

    if (!has_open_neighbor) {
      atom->setHybridization(hybrid == Hybridization::SP2 ? Hybridization::SP3 : Hybridization::SP2);
    }
  }
}

void apply_bond_assignment_patterns(RDKit::RWMol &mol) {
  RDKit::SubstructMatchParameters params;
  params.maxMatches = 100000;

  const auto &patterns = bond_assignment_patterns();
  const auto &queries = bond_assignment_queries();

  for (size_t i = 0; i < patterns.size(); ++i) {
    const auto matches = RDKit::SubstructMatch(mol, *queries[i], params);
    const auto &triplets = patterns[i].bond_triplets;
    for (const auto &match : matches) {
      for (size_t j = 0; j + 2 < triplets.size(); j += 3) {
        const int query_idx_1 = triplets[j];
        const int query_idx_2 = triplets[j + 1];
        const auto bond_type = static_cast<BondType>(triplets[j + 2]);
        RDKit::Bond *bond = mol.getBondBetweenAtoms(match[query_idx_1].second, match[query_idx_2].second);
        if (bond) bond->setBondType(bond_type);
      }
    }
  }
}

const RDKit::ROMol &compile_smarts_once(const char *smarts) {
  thread_local std::vector<std::pair<std::string, std::shared_ptr<RDKit::ROMol>>> cache;
  auto it = std::find_if(cache.begin(), cache.end(), [smarts](const auto &entry) { return entry.first == smarts; });
  if (it == cache.end()) {
    cache.emplace_back(smarts, std::shared_ptr<RDKit::ROMol>(RDKit::SmartsToMol(smarts)));
    return *cache.back().second;
  }
  return *it->second;
}

RDKit::Atom *match_atom(RDKit::RWMol &mol, const RDKit::MatchVectType &match, int query_idx) {
  // We can optimize this, but it's not going to make a difference performance-wise
  for (const auto &pair : match) {
    if (pair.first == query_idx) {
      return mol.getAtomWithIdx(pair.second);
    }
  }
  return nullptr;
}

bool assign_double_bond_if_needed(RDKit::RWMol &mol, const RDKit::Atom &a1, const RDKit::Atom &a2) {
  if (atom_has_bond_type(a1, BondType::DOUBLE)) return false;
  if (auto *bond = mol.getBondBetweenAtoms(a1.getIdx(), a2.getIdx())) {
    bond->setBondType(BondType::DOUBLE);
    return true;
  }
  return false;
}

void adjust_carbonyls(RDKit::RWMol &mol, RDKit::Conformer &conf) {
  const auto &pattern = compile_smarts_once("[#8D1;!-][#6](*)(*)");
  std::vector<RDKit::MatchVectType> matches;
  RDKit::SubstructMatch(mol, pattern, matches);
  constexpr double limit_sq = 1.28 * 1.28;  // 1.6384

  for (const auto &match : matches) {
    RDKit::Atom *oxygen = match_atom(mol, match, 0);
    RDKit::Atom *carbon = match_atom(mol, match, 1);
    if (!oxygen || !carbon) continue;

    if (bond_length_sq(conf, *oxygen, *carbon) >= limit_sq) continue;

    const double avg_angle = average_bond_angle(*oxygen);
    if (avg_angle > 115.0 && avg_angle < 150.0) {
      assign_double_bond_if_needed(mol, *oxygen, *carbon);
    }
  }
}

void adjust_thiones(RDKit::RWMol &mol, RDKit::Conformer &conf) {
  const auto &pattern = compile_smarts_once("[#16D1][#6](*)(*)");
  std::vector<RDKit::MatchVectType> matches;
  RDKit::SubstructMatch(mol, pattern, matches);
  constexpr double limit_sq = 1.72 * 1.72;  // 2.9584

  for (const auto &match : matches) {
    RDKit::Atom *sulfur = match_atom(mol, match, 0);
    RDKit::Atom *carbon = match_atom(mol, match, 1);
    if (!sulfur || !carbon) continue;

    if (bond_length_sq(conf, *sulfur, *carbon) >= limit_sq) continue;

    const double avg_angle = average_bond_angle(*sulfur);
    if (avg_angle > 115.0 && avg_angle < 150.0) {
      assign_double_bond_if_needed(mol, *sulfur, *carbon);
    }
  }
}

void adjust_isocyanates(RDKit::RWMol &mol, RDKit::Conformer &conf) {
  const auto &pattern = compile_smarts_once("[#8,#16;D1][#6D2][#7D2]");
  std::vector<RDKit::MatchVectType> matches;
  RDKit::SubstructMatch(mol, pattern, matches);
  constexpr double carbon_nitrogen_limit_sq = 1.34 * 1.34;  // 1.7956

  for (const auto &match : matches) {
    RDKit::Atom *terminal = match_atom(mol, match, 0);
    RDKit::Atom *carbon   = match_atom(mol, match, 1);
    RDKit::Atom *nitrogen = match_atom(mol, match, 2);
    if (!terminal || !carbon || !nitrogen) continue;

    const double limit = terminal->getAtomicNum() == 8 ? 1.28 : 1.72;
    const double limit_sq = limit * limit;
    if (bond_length_sq(conf, *terminal, *carbon) >= limit_sq) continue;
    if (bond_length_sq(conf, *carbon, *nitrogen) >= carbon_nitrogen_limit_sq) continue;

    const double avg_angle = average_bond_angle(*carbon);
    if (avg_angle > 150.0) {
      assign_double_bond_if_needed(mol, *terminal, *carbon);
      assign_double_bond_if_needed(mol, *carbon, *nitrogen);
    }
  }
}

void adjust_oximes(RDKit::RWMol &mol, RDKit::Conformer &conf) {
  const auto &pattern = compile_smarts_once("[#6D3][#7D2][#8D2]");
  std::vector<RDKit::MatchVectType> matches;
  RDKit::SubstructMatch(mol, pattern, matches);
  constexpr double limit_sq = 1.4 * 1.4;  // 1.96

  for (const auto &match : matches) {
    RDKit::Atom *carbon   = match_atom(mol, match, 0);
    RDKit::Atom *nitrogen = match_atom(mol, match, 1);
    if (!carbon || !nitrogen) continue;

    if (bond_length_sq(conf, *carbon, *nitrogen) >= limit_sq) continue;

    const double avg_angle = average_bond_angle(*carbon);
    if (avg_angle > 110.0 && avg_angle < 150.0) {
      assign_double_bond_if_needed(mol, *carbon, *nitrogen);
    }
  }
}

void adjust_oxidopyridines(RDKit::RWMol &mol, RDKit::Conformer &conf) {
  const auto &pattern = compile_smarts_once("[#8D1][#7D3r6]");
  std::vector<RDKit::MatchVectType> matches;
  RDKit::SubstructMatch(mol, pattern, matches);
  constexpr double limit_sq = 1.35 * 1.35;  // 1.8225

  for (const auto &match : matches) {
    RDKit::Atom *oxygen   = match_atom(mol, match, 0);
    RDKit::Atom *nitrogen = match_atom(mol, match, 1);
    if (!oxygen || !nitrogen) continue;

    if (bond_length_sq(conf, *oxygen, *nitrogen) >= limit_sq) continue;

    const double avg_angle = average_bond_angle(*oxygen);
    if (avg_angle > 110.0 && avg_angle < 150.0) {
      oxygen->setFormalCharge(-1);
      nitrogen->setFormalCharge(+1);
    }
  }
}

void apply_functional_group_adjustments(RDKit::RWMol &mol, RDKit::Conformer &conf) {
  adjust_carbonyls(mol, conf);
  adjust_thiones(mol, conf);
  adjust_isocyanates(mol, conf);
  adjust_oximes(mol, conf);
  adjust_oxidopyridines(mol, conf);
}

void mark_candidate_aromatic_rings(RDKit::RWMol &mol) {
  const auto *ring_info = mol.getRingInfo();
  if (!ring_info) return;

  for (const auto &ring : ring_info->atomRings()) {
    const auto ring_size = ring.size();
    if (ring_size < 5 || ring_size > 7) continue;

    bool requires_typing = false;
    for (int atom_idx : ring) {
      RDKit::Atom *atom = mol.getAtomWithIdx(atom_idx);
      if (atom->getHybridization() != Hybridization::SP2) {
        requires_typing = true;
        break;
      }
      for (const auto *bond : mol.atomBonds(atom)) {
        const auto type = bond->getBondType();
        if (type == BondType::DOUBLE || type == BondType::TRIPLE) {
          requires_typing = true;
          break;
        }
      }
      if (requires_typing) break;
    }

    if (requires_typing) continue;

    // Mark all ring bonds as aromatic
    for (size_t i = 0; i < ring_size; ++i) {
      RDKit::Bond *bond = mol.getBondBetweenAtoms(ring[i], ring[(i + 1) % ring_size]);
      if (bond) {
        bond->setIsAromatic(true);
      }
    }
  }
}

} // namespace

void perceive_bond_orders_obabel(RDKit::RWMol &mol) {
  mol.updatePropertyCache(false);
  assign_initial_hybridization(mol);

  refine_hybridization_by_geometry(mol);

  RDKit::Conformer &conf = mol.getConformer();
  refine_ring_hybridization(mol, conf);
  relax_terminal_assignments(mol);

  apply_bond_assignment_patterns(mol);
  apply_functional_group_adjustments(mol, conf);
  mark_candidate_aromatic_rings(mol);

  mol.updatePropertyCache(false);
}

} // namespace lahuta
