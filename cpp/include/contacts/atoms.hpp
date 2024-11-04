#ifndef LAHUTA_CONTACTS_HPP
#define LAHUTA_CONTACTS_HPP

#include "GraphMol/RWMol.h"
#include "atom_types.hpp"
#include "contacts/halogen_bonds.hpp"
#include "contacts/hydrophobic.hpp"
#include "contacts/metals.hpp"
#include "hydrogen_bonds.hpp"
#include "valence_model.hpp"

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
    AtomTypeStrategy composite;

    composite.add_strategy<HBondDonorAtom>();
    composite.add_strategy<HBondAcceptorAtom>();
    composite.add_strategy<WeakHBondDonorAtom>();
    composite.add_strategy<HydrophobicAtom>();
    composite.add_strategy<HalogenDonorAtom>();
    composite.add_strategy<HalogenAcceptorAtom>();
    composite.add_strategy<MetalAtom>();
    composite.add_strategy<MetalBindingAtom>();

    return composite;
  }
};

class AtomTypeAnalysis {
public:
  static std::vector<AtomType> analyze(const RDKit::RWMol &mol) {
    std::vector<AtomType> atom_types = {mol.getNumAtoms(), AtomType::NONE};

    auto strategy = AtomTypeFactory::create();

    bool assign_charge = true, assign_h = true;
    ValenceModel valence_model{assign_charge, assign_h};
    valence_model.apply(mol);

    for (const auto &atom : mol.atoms()) {
      // FIX: double check assign_h default value and how to handle it
      atom_types[atom->getIdx()] = strategy.identify(mol, *atom);
    }
    return atom_types;
  }
};

} // namespace lahuta

#endif // LAHUTA_CONTACTS_HPP
