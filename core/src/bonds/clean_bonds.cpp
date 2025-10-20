#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <vector>

#include <rdkit/GraphMol/AtomIterators.h>
#include <rdkit/GraphMol/BondIterators.h>
#include <rdkit/GraphMol/PeriodicTable.h>

#include "clean_bonds.hpp"

// clang-format off
namespace lahuta {
namespace {

constexpr double Pi = 3.14159265358979323846;
constexpr double RadiansToDegrees = 180.0 / Pi;
constexpr double AngleCutoffDeg   = 40.0;
constexpr double NormEpsilon      = 1e-12;

const std::array<int, 119> &biological_sigma_limits() {
  static const std::array<int, 119> limits = [] {
    std::array<int, 119> data{};
    data.fill(-1);

    data[1] = 1;  // H
    data[5] = 4;  // B
    data[6] = 4;  // C
    data[7] = 4;  // N (permits quaternary ammonium)
    data[8] = 2;  // O
    data[9] = 1;  // F
    data[15] = 6; // P (allows hypervalent intermediates)
    data[16] = 6; // S (covers sulfonium/sulfate)
    data[17] = 1; // Cl
    data[34] = 2; // Se
    data[35] = 1; // Br
    data[53] = 1; // I

    return data;
  }();
  return limits;
}

const RDKit::PeriodicTable &periodic_table() {
  static const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
  return *tbl;
}

unsigned allowed_sigma_bonds(const RDKit::Atom &atom) {
  const unsigned atomic_num = atom.getAtomicNum();
  if (atomic_num == 0) return atom.getDegree();

  const auto &limits = biological_sigma_limits();
  if (atomic_num < limits.size()) {
    const int predefined = limits[atomic_num];
    if (predefined >= 0) {
      return static_cast<unsigned>(predefined);
    }
  }

  unsigned fallback = 0;
  for (int valence : periodic_table().getValenceList(atomic_num)) {
    if (valence >= 0) {
      fallback = std::max(fallback, static_cast<unsigned>(valence));
    }
  }

  if (fallback == 0) return atom.getDegree();
  return fallback;
}

double squared_distance(const RDKit::Conformer &conf, unsigned a_idx, unsigned b_idx) {
  const auto &a = conf.getAtomPos(a_idx);
  const auto &b = conf.getAtomPos(b_idx);
  const double dx = a.x - b.x;
  const double dy = a.y - b.y;
  const double dz = a.z - b.z;
  return dx * dx + dy * dy + dz * dz;
}

double angle_degrees(const RDKit::Conformer &conf, unsigned center_idx, unsigned first_idx, unsigned second_idx) {
  const auto &pc = conf.getAtomPos(center_idx);
  const auto &p1 = conf.getAtomPos(first_idx);
  const auto &p2 = conf.getAtomPos(second_idx);

  const double v1x = p1.x - pc.x;
  const double v1y = p1.y - pc.y;
  const double v1z = p1.z - pc.z;

  const double v2x = p2.x - pc.x;
  const double v2y = p2.y - pc.y;
  const double v2z = p2.z - pc.z;

  const double len1_sq = v1x * v1x + v1y * v1y + v1z * v1z;
  const double len2_sq = v2x * v2x + v2y * v2y + v2z * v2z;

  if (len1_sq < NormEpsilon || len2_sq < NormEpsilon) {
    return 180.0;
  }

  const double inv_norm = 1.0 / (std::sqrt(len1_sq) * std::sqrt(len2_sq));
  const double dot = v1x * v2x + v1y * v2y + v1z * v2z;
  const double cosine = std::clamp(dot * inv_norm, -1.0, 1.0);

  return std::acos(cosine) * RadiansToDegrees;
}

double removal_priority(const RDKit::Atom *center, const RDKit::Bond &bond, const RDKit::Conformer &conf) {
  const auto *nbr = bond.getOtherAtom(center);
  const bool neighbor_is_h = nbr && nbr->getAtomicNum() == 1;
  const bool center_is_h = center->getAtomicNum() == 1;

  double priority = squared_distance(conf, bond.getBeginAtomIdx(), bond.getEndAtomIdx());

  if (neighbor_is_h) priority += 25.0;
  if (center_is_h && neighbor_is_h) priority += 10.0;

  switch (bond.getBondType()) {
    case RDKit::Bond::BondType::AROMATIC:
    case RDKit::Bond::BondType::DOUBLE:
    case RDKit::Bond::BondType::TRIPLE:
    case RDKit::Bond::BondType::QUADRUPLE:
      priority -= 6.0;
      break;
    default:
      break;
  }

  return priority;
}

RDKit::Bond *choose_bond_to_remove(const RDKit::Atom *center, const std::vector<RDKit::Bond *> &candidates, const RDKit::Conformer &conf) {
  if (candidates.empty()) return nullptr;

  if (center->getAtomicNum() == 1) {
    auto it = std::find_if(candidates.begin(), candidates.end(), [center](const RDKit::Bond *bond) {
      const RDKit::Atom *other = bond->getOtherAtom(center);
      return other && other->getAtomicNum() == 1;
    });
    if (it != candidates.end()) {
      return *it;
    }
  }

  return *std::max_element(
      candidates.begin(),
      candidates.end(),
      [center, &conf](const RDKit::Bond *lhs, const RDKit::Bond *rhs) {
        return removal_priority(center, *lhs, conf) < removal_priority(center, *rhs, conf);
      });
}

bool prune_excess_valence(const RDKit::Atom *atom, RDKit::RWMol &mol, RDKit::Conformer &conf) {
  std::vector<RDKit::Bond *> bonds;
  bonds.reserve(atom->getDegree());
  for (const auto &bond : mol.atomBonds(atom)) {
    bonds.push_back(bond);
  }
  RDKit::Bond *selected = choose_bond_to_remove(atom, bonds, conf);
  if (!selected) return false;

  mol.removeBond(selected->getBeginAtomIdx(), selected->getEndAtomIdx());
  return true;
}

struct AngleCandidate {
  double angle_deg    = 180.0;
  RDKit::Bond *first  = nullptr;
  RDKit::Bond *second = nullptr;
};

std::optional<AngleCandidate> tightest_angle(const RDKit::Atom *atom, const RDKit::RWMol &mol, const RDKit::Conformer &conf) {
  const auto degree = atom->getDegree();
  if (degree < 2) return std::nullopt;

  std::vector<RDKit::Bond *> bonds;
  bonds.reserve(degree);
  for (const auto &bond : mol.atomBonds(atom)) {
    bonds.push_back(bond);
  }

  AngleCandidate best{};
  best.angle_deg = 180.0;

  for (size_t i = 0; i < bonds.size(); ++i) {
    const RDKit::Atom *first = bonds[i]->getOtherAtom(atom);
    if (!first) continue;
    for (size_t j = i + 1; j < bonds.size(); ++j) {
      const RDKit::Atom *second = bonds[j]->getOtherAtom(atom);
      if (!second) continue;

      double angle = angle_degrees(conf, atom->getIdx(), first->getIdx(), second->getIdx());
      if (angle < best.angle_deg) {
        best.angle_deg = angle;
        best.first     = bonds[i];
        best.second    = bonds[j];
      }
    }
  }

  if (!best.first || !best.second) return std::nullopt;
  return best;
}

bool prune_crowded_geometry(const RDKit::Atom *atom, RDKit::RWMol &mol, RDKit::Conformer &conf) {
  auto candidate = tightest_angle(atom, mol, conf);
  if (!candidate || candidate->angle_deg >= AngleCutoffDeg) return false;

  std::vector<RDKit::Bond *> options{candidate->first, candidate->second};
  RDKit::Bond *selected = choose_bond_to_remove(atom, options, conf);
  if (!selected) return false;

  mol.removeBond(selected->getBeginAtomIdx(), selected->getEndAtomIdx());
  return true;
}

} // namespace

void clean_bonds(RDKit::RWMol &mol, RDKit::Conformer &conf) {
  if (mol.getNumBonds() == 0) return;

  bool updated;
  do {
    updated = false;

    for (const auto &atom : mol.atoms()) {
      if (!atom) continue;

      while (true) {
        const unsigned sigma_limit = allowed_sigma_bonds(*atom);
        bool removed = false;

        if (atom->getDegree() > sigma_limit) {
          removed = prune_excess_valence(atom, mol, conf);
        } else {
          removed = prune_crowded_geometry(atom, mol, conf);
        }

        if (!removed) break;

        updated = true;
        if (atom->getDegree() == 0) break;

      }
    }
  } while (updated);
}

} // namespace lahuta
