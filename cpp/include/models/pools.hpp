#ifndef LAHUTA_POOLS_HPP
#define LAHUTA_POOLS_HPP

#include "GraphMol/Atom.h"
#include "GraphMol/Bond.h"
#include "GraphMol/MonomerInfo.h"
#include "models/obj_pool.hpp"

namespace lahuta {

struct AtomPoolTraits {
  static RDKit::Atom *create() { return new RDKit::Atom(); }
  static RDKit::Atom *create(int atomic_num) { return new RDKit::Atom(atomic_num); }

  static void reset(RDKit::Atom *atom) { atom->resetState(); }
};

class AtomPool {
public:
  explicit AtomPool(std::size_t initial_capacity = 2000) : pool_(initial_capacity) {}

  auto *createAtom(int atomic_num) { return pool_.create(atomic_num); }
  auto *createAtom() { return pool_.create(); }

  auto prepareAtoms(std::size_t count) {
    auto atoms = pool_.prepare(count);
    return atoms;
  }

  void reset() { pool_.reset(); }
  void clear() { pool_.clear(); }

private:
  ObjectPool<RDKit::Atom, AtomPoolTraits> pool_;
};

struct BondPoolTraits {
  static RDKit::Bond *create() { return new RDKit::Bond(); }

  static void reset(RDKit::Bond *bond) { bond->resetState(); }
};

// NOTE: calls to BondPool::createBond() seem to be more expensive than their counterparts.
//       Could we be making unintended copies somewhere?
class BondPool {
public:
  explicit BondPool(std::size_t initial_capacity = 2000) : pool_(initial_capacity) {}

  auto *createBond(int begin_idx, int end_idx, RDKit::Bond::BondType bond_type) {

    // NOTE: we have some redundancy on the assignment here, but its' irrelevant
    RDKit::Bond *b = pool_.create();
    b->setBeginAtomIdx(begin_idx);
    b->setEndAtomIdx(end_idx);
    b->setBondType(bond_type);
    return b;
  }

  auto prepareBonds(const std::vector<std::tuple<int, int, RDKit::Bond::BondType>> &bond_data) {
    auto bonds = pool_.prepare(bond_data.size());
    for (std::size_t i = 0; i < bond_data.size(); ++i) {
      const auto &[begin_idx, end_idx, bond_type] = bond_data[i];
      RDKit::Bond *b = bonds[i];
      // BondPoolTraits::reset() is already called inside prepare(...),

      b->setBeginAtomIdx(begin_idx);
      b->setEndAtomIdx(end_idx);
      b->setBondType(bond_type);
    }
    return bonds;
  }

  void reset() { pool_.reset(); }
  void clear() { pool_.clear(); }

private:
  ObjectPool<RDKit::Bond, BondPoolTraits> pool_;
};

struct InfoPoolTraits {
  static auto *create() { return new RDKit::pAtomPDBResidueInfo(); }

  static void reset(RDKit::pAtomPDBResidueInfo *info) {
    info->setName("");
    info->setSerialNumber(-1);
    info->setResidueName("");
    info->setResidueNumber(-1);
  }
};

class InfoPool {
public:
  explicit InfoPool(std::size_t initial_capacity = 1000) : pool_(initial_capacity) {}

  auto *createAtomInfo(const char *atom_name, int serial, const char *res_name, int res_number) {
    RDKit::pAtomPDBResidueInfo *info = pool_.create();
    info->setName(atom_name);
    info->setSerialNumber(serial);
    info->setResidueName(res_name);
    info->setResidueNumber(res_number);
    return info;
  }

  auto prepareInfos(std::size_t count) { return pool_.prepare(count); }

  void reset() { pool_.reset(); }
  void clear() { pool_.clear(); }

private:
  ObjectPool<RDKit::pAtomPDBResidueInfo, InfoPoolTraits> pool_;
};

} // namespace lahuta

#endif // LAHUTA_POOLS_HPP
