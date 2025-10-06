#include "subgraph.hpp"

namespace lahuta::bonds::subgraph {

RDKit::RWMol build_rdkit_submol(const RDKit::RWMol &source, span<const int> indices, bool include_bonds) {
  RDKit::RWMol sub;
  const RDKit::Conformer &conf = source.getConformer();
  auto *new_conf = new RDKit::Conformer();

  std::vector<int> map(source.getNumAtoms(), -1);
  map.reserve(source.getNumAtoms());

  for (size_t i = 0; i < indices.size(); ++i) {
    int src_idx = indices[i];
    const RDKit::Atom *src_atom = source.getAtomWithIdx(src_idx);
    auto *new_atom = new RDKit::Atom(src_atom->getAtomicNum());
    int dst_idx = sub.addAtom(new_atom, false, true);
    map[src_idx] = dst_idx;
    new_conf->setAtomPos(dst_idx, conf.getAtomPos(src_idx));
  }

  if (include_bonds) {
    for (const auto &bond : source.bonds()) {
      int a = bond->getBeginAtomIdx();
      int b = bond->getEndAtomIdx();
      if (a < 0 || b < 0 || a >= static_cast<int>(map.size()) || b >= static_cast<int>(map.size())) continue;
      int na = map[a];
      int nb = map[b];
      if (na >= 0 && nb >= 0) {
        sub.addBond(na, nb, bond->getBondType());
      }
    }
  }

  sub.addConformer(new_conf, true);
  // sub.updatePropertyCache(false);
  return sub;
}

} // namespace lahuta::bonds::subgraph
