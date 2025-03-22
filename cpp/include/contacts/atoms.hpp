#ifndef LAHUTA_CONTACTS_HPP
#define LAHUTA_CONTACTS_HPP

#include "GraphMol/RWMol.h"
#include "atom_types.hpp"
#include "contacts/halogen_bonds.hpp"
#include "contacts/hydrophobic.hpp"
#include "contacts/metals.hpp"
#include "hydrogen_bonds.hpp"

namespace lahuta {

class AtomTypeBase {
public:
  virtual ~AtomTypeBase() = default;
  virtual AtomType identify(const RDKit::RWMol &mol, const RDKit::Atom &atom) const = 0;
};

class HBondAcceptorAtom : public AtomTypeBase {
public:
  AtomType identify(const RDKit::RWMol &mol, const RDKit::Atom &atom) const override {
    return add_hydrogen_acceptor(mol, atom);
  }
};

class HBondDonorAtom : public AtomTypeBase {
public:
  AtomType identify(const RDKit::RWMol &mol, const RDKit::Atom &atom) const override {
    return add_hydrogen_donor(mol, atom);
  }
};

class WeakHBondDonorAtom : public AtomTypeBase {
public:
  AtomType identify(const RDKit::RWMol &mol, const RDKit::Atom &atom) const override {
    return add_weak_hydrogen_donor(mol, atom);
  }
};

class HydrophobicAtom : public AtomTypeBase {
public:
  AtomType identify(const RDKit::RWMol &mol, const RDKit::Atom &atom) const override {
    return add_hydrophobic_atom(mol, atom);
  }
};

class HalogenDonorAtom : public AtomTypeBase {
public:
  AtomType identify(const RDKit::RWMol &mol, const RDKit::Atom &atom) const override {
    return add_halogen_donor(mol, atom);
  }
};

class HalogenAcceptorAtom : public AtomTypeBase {
public:
  AtomType identify(const RDKit::RWMol &mol, const RDKit::Atom &atom) const override {
    return add_halogen_acceptor(mol, atom);
  }
};

class MetalAtom : public AtomTypeBase {
public:
  AtomType identify(const RDKit::RWMol &mol, const RDKit::Atom &atom) const override {
    return add_metal(mol, atom);
  }
};

class MetalBindingAtom : public AtomTypeBase {
public:
  AtomType identify(const RDKit::RWMol &mol, const RDKit::Atom &atom) const override {
    return add_metal_binding(mol, atom);
  }
};

class AtomTypeStrategy {
private:
  std::vector<std::unique_ptr<AtomTypeBase>> strategies;

public:
  template <typename T> void add_strategy() {
    static_assert(std::is_base_of<AtomTypeBase, T>::value, "T must derive from AtomTypeStrategy");
    strategies.push_back(std::make_unique<T>());
  }

  AtomType identify(const RDKit::RWMol &mol, const RDKit::Atom &atom) const {
    AtomType result = AtomType::NONE;
    for (const auto &strategy : strategies) {
      result |= strategy->identify(mol, atom);
    }
    return result;
  }
};

class AtomTypeFactory {
public:
  static AtomTypeStrategy create() {
    AtomTypeStrategy at;

    at.add_strategy<HBondDonorAtom>();
    at.add_strategy<HBondAcceptorAtom>();
    at.add_strategy<WeakHBondDonorAtom>();
    at.add_strategy<HydrophobicAtom>();
    at.add_strategy<HalogenDonorAtom>();
    at.add_strategy<HalogenAcceptorAtom>();
    at.add_strategy<MetalAtom>();
    at.add_strategy<MetalBindingAtom>();

    return at;
  }
};

class AtomTypeAnalysis {
public:
  static AtomEntityCollection analyze(const RDKit::RWMol &mol) {
    AtomEntityCollection atom_types;
    atom_types.reserve(mol.getNumAtoms());

    auto strategy = AtomTypeFactory::create();

    // FIX: We ignore here atom typing added by OpenBabel typing system.
    for (const auto &atom : mol.atoms()) {
      AtomType at = strategy.identify(mol, *atom);
      atom_types.add_data(mol, atom, at);
    }
    return atom_types;
  }
};

} // namespace lahuta

#endif // LAHUTA_CONTACTS_HPP
