#include "luni_props.hpp"
#include "lahuta.hpp"
#include "properties/registry.hpp"
#include "properties/types.hpp"

// clang-format off

namespace lahuta {

template <typename T>
std::vector<T> atom_attrs(const Luni &luni, std::function<T(const RDKit::Atom *)> func) {
  std::vector<T> attrs;
  auto &mol = luni.get_molecule();
  attrs.reserve(mol.getNumAtoms());
  for (const auto *atom : mol.atoms()) {
    attrs.push_back(func(atom));
  }
  return attrs;
}

void LuniProperties::register_all() {

  PropertyRegistry::register_property<Luni, PropertyKey::Indices>(
      [](const Luni &L) { return atom_attrs<int>(L, [](const RDKit::Atom *a) { return a->getIdx(); }); });

  PropertyRegistry::register_property<Luni, PropertyKey::Names>([](const Luni &L) {
    return atom_attrs<std::string>(L, [&L](const RDKit::Atom *a) -> const std::string & {
      return L.get_info(a->getIdx())->getName();
    });
  });

  PropertyRegistry::register_property<Luni, PropertyKey::Elements>([](const Luni &L) {
    const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
    return atom_attrs<std::string>(L, [&tbl](const RDKit::Atom *a) {
      return tbl->getElementSymbol(a->getAtomicNum());
    });
  });

  PropertyRegistry::register_property<Luni, PropertyKey::Positions>(
      [](const Luni &L) -> std::vector<RDGeom::Point3D> {
        /*return std::move(l.positions());*/
        return L.positions(); // FIX: Return a copy
      });
}

} // namespace lahuta
