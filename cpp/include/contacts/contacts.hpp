#ifndef LAHUTA_CONTACTS_HPP
#define LAHUTA_CONTACTS_HPP

#include "Geometry/point.h"
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
  // Add strategy to the composite
  template <typename T> void add_strategy() {
    static_assert(std::is_base_of<AtomTypeBase, T>::value, "T must derive from AtomTypeStrategy");
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

class AtomTypeFactory {
public:
  static AtomTypeStrategy create() {
    AtomTypeStrategy composite;

    // default strategies
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
      // FIX: potentially expose options to the user
      atom_types[atom->getIdx()] = strategy.identify(mol, *atom);
    }
    return atom_types;
  }
};

struct AtomData {
  AtomType type;
  const RDKit::Atom *atom;
  const RDGeom::Point3D *pos;
  size_t idx;

  size_t get_idx() const { return idx; }
  int get_id() const { return idx; }

  /*AtomData() = default;*/
  AtomData(AtomType type, const RDKit::Atom *atom, const RDGeom::Point3D *pos, size_t idx)
      : type(type), atom(atom), pos(pos), idx(idx) {}
  /*AtomData(const AtomData &other) = default;*/
  /*AtomData(AtomData &&other) = default;*/
  /*AtomData &operator=(const AtomData &other) = default;*/
  /*AtomData &operator=(AtomData &&other) = default;*/

};

struct AtomDataVec {
private:
  std::vector<AtomData> data;
public:
  AtomDataVec() = default;
  /*AtomDataVec(const std::vector<AtomData> &data_) : data(data_) {}*/

  const std::vector<AtomData> &get_data() const { return data; }

  /*void add_data(const AtomData &atom_data) { data.push_back(atom_data); }*/
  /*void add_data(AtomType type, const RDKit::Atom *atom, const RDGeom::Point3D *pos, size_t idx) {*/
  /*  data.push_back(AtomData(type, atom, pos, idx));*/
  /*}*/
  void add_data(const RDKit::RWMol &mol, const RDKit::Atom *atom, AtomType type) {
    data.push_back(AtomData(type, atom, &mol.getConformer().getAtomPos(atom->getIdx()), atom->getIdx()));
  }

  AtomData &operator[](size_t index) { return data[index]; }
  const AtomData &operator[](size_t index) const { return data[index]; }
  int size() const { return data.size(); }

  RDGeom::POINT3D_VECT positions() const {
    RDGeom::POINT3D_VECT pos;
    pos.reserve(data.size());
    for (const auto& atom_data : data) {
      pos.push_back(*atom_data.pos);
    }
    return pos;
  }

  std::vector<size_t> atom_ids() const {
    std::vector<size_t> ids;
    ids.reserve(data.size());
    for (const auto& atom_data : data) {
      ids.push_back(atom_data.idx);
    }
    return ids;
  }

  const std::vector<const RDKit::Atom*> atoms() const {
    std::vector<const RDKit::Atom*> atoms_vec;
    atoms_vec.reserve(data.size());
    for (const auto& atom_data : data) {
      atoms_vec.push_back(atom_data.atom);
    }
    return atoms_vec;
  }
};

using FeatureTypeCheckFunc = std::function<bool(const AtomType &, const AtomType &)>;
const AtomDataVec get_atom_data(const Luni* luni, AtomType type, FeatureTypeCheckFunc check_func = AtomTypeFlags::has_any);
std::vector<const RDKit::Atom *> get_atom_types(const Luni *luni, AtomType type);

} // namespace lahuta

#endif // LAHUTA_CONTACTS_HPP
