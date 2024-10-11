#ifndef LAHUTA_CONTACTS_HPP
#define LAHUTA_CONTACTS_HPP

#include "GraphMol/RWMol.h"
#include "atom_types.hpp"
#include "hydrogen_bonds.hpp"
#include "valence.hpp"

namespace lahuta {

struct AtomTypeParams {
  AtomType atom_type;
  std::string name;
  std::vector<AtomType> components;
};

class AtomTypeStrategy {
public:
  virtual ~AtomTypeStrategy() = default;
  virtual AtomType identify(const RDKit::RWMol &mol,
                            const RDKit::Atom &atom) const = 0;
};

class HBondAcceptorStrategy : public AtomTypeStrategy {
public:
  AtomType identify(const RDKit::RWMol &mol,
                    const RDKit::Atom &atom) const override {
    return add_hydrogen_acceptor(mol, atom);
  }
};

class HBondDonorStrategy : public AtomTypeStrategy {
public:
  AtomType identify(const RDKit::RWMol &mol,
                    const RDKit::Atom &atom) const override {
    return add_hydrogen_donor(mol, atom);
  }
};

class CompositeAtomTypeStrategy {
private:
  std::vector<std::unique_ptr<AtomTypeStrategy>> strategies;

public:
  // Add strategy to the composite
  template <typename T> void addStrategy() {
    static_assert(std::is_base_of<AtomTypeStrategy, T>::value,
                  "T must derive from AtomTypeStrategy");
    strategies.push_back(std::make_unique<T>());
  }

  // Apply strategies and combine results
  AtomType identify(const RDKit::RWMol &mol, const RDKit::Atom &atom) const {
    AtomType result = AtomType::NONE;
    for (const auto &strategy : strategies) {
      result |= strategy->identify(mol, atom);
    }
    return result;
  }
};

class AtomTypeStrategyFactory {
public:
  static CompositeAtomTypeStrategy create() {
    CompositeAtomTypeStrategy composite;
    // only two for now
    composite.addStrategy<HBondAcceptorStrategy>();
    composite.addStrategy<HBondDonorStrategy>();
    return composite;
  }
};

class AtomTypeAnalyzer {
public:
  static std::vector<AtomType> analyzeAtomTypes(const RDKit::RWMol &mol) {
    std::vector<AtomType> atom_types = {mol.getNumAtoms(), AtomType::NONE};

    auto strategy = AtomTypeStrategyFactory::create();

    for (const auto &atom : mol.atoms()) {
      compute_valence(mol, *atom, true, true);
      atom_types[atom->getIdx()] = strategy.identify(mol, *atom);
    }
    return atom_types;
  }
};

} // namespace lahuta

#endif // LAHUTA_CONTACTS_HPP
