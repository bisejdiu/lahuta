#include "lahuta.hpp"

namespace lahuta {

const std::vector<std::string> Luni::symbols() const {
  return atom_attrs<std::string>(
      [](const RDKit::Atom *atom) { return atom->getSymbol(); });
}

const std::vector<std::string> Luni::names() const {
  return atom_attrs<std::string>(
      [](const RDKit::Atom *atom) -> const std::string & {
        auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(
            atom->getMonomerInfo());
        return info->getName();
      });
}

const std::vector<int> Luni::indices() const {
  return atom_attrs<int>(
      [](const RDKit::Atom *atom) { return atom->getIdx(); });
}

const std::vector<int> Luni::atomic_numbers() const {
  return atom_attrs<int>(
      [](const RDKit::Atom *atom) { return atom->getAtomicNum(); });
}

const std::vector<std::string> Luni::elements() const {
  const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
  return atom_attrs<std::string>([&tbl](const RDKit::Atom *atom) {
    return tbl->getElementSymbol(atom->getAtomicNum());
  });
}

const std::vector<std::string> Luni::resnames() const {
  return atom_attrs<std::string>(
      [](const RDKit::Atom *atom) -> const std::string & {
        auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(
            atom->getMonomerInfo());
        return info->getResidueName();
      });
}

const std::vector<int> Luni::resids() const {
  return atom_attrs<int>([](const RDKit::Atom *atom) -> int {
    auto *info =
        static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
    return info->getResidueNumber();
  });
}

const std::vector<int> Luni::resindices() const {
  return atom_attrs<int>([](const RDKit::Atom *atom) -> int {
    auto *info =
        static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
    return info->getSegmentNumber();
  });
}

const std::vector<std::string> Luni::chainlabels() const {
  return atom_attrs<std::string>(
      [](const RDKit::Atom *atom) -> const std::string & {
        auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(
            atom->getMonomerInfo());
        return info->getChainId();
      });
}

template <typename T>
std::vector<T>
Luni::atom_attrs(std::function<T(const RDKit::Atom *)> func) const {
  std::vector<T> attrs;
  attrs.reserve(mol->getNumAtoms());
  for (const auto atom : mol->atoms()) {
    attrs.push_back(func(atom));
  }
  return attrs;
}

template <typename T>
std::vector<std::reference_wrapper<const T>>
Luni::atom_attrs_ref(std::function<const T &(const RDKit::Atom *)> func) const {
  std::vector<std::reference_wrapper<const T>> attributes;
  attributes.reserve(mol->getNumAtoms());
  for (const auto atom : mol->atoms()) {
    attributes.push_back(std::cref(func(atom)));
  }
  return attributes;
}

std::vector<RDKit::MatchVectType>
Luni::match_smarts_string(std::string sm, std::string atype,
                          bool log_values) const {

  if (!mol->getRingInfo()->isInitialized()) {
    RDKit::MolOps::symmetrizeSSSR(*mol);
  }
  std::vector<RDKit::MatchVectType> match_list;
  auto sm_mol = RDKit::SmartsToMol(sm);
  mol->updatePropertyCache(false);
  RDKit::SubstructMatch(*mol, *sm_mol, match_list);
  return match_list;
};

NSResults Luni::find_neighbors_opt(double cutoff) {

  if (cutoff == _cutoff) {
    return bonded_nps;
  } else if (cutoff < _cutoff) {
    return bonded_nps.filter(cutoff);
  }

  grid.update_cutoff(cutoff);
  auto ns = grid.self_search();
  return ns;
}

NSResults Luni::remove_adjascent_residueid_pairs(NSResults &results,
                                                 int res_diff) {
  Pairs filtered;
  Distances dists;
  for (size_t i = 0; i < results.get_pairs().size(); ++i) {
    auto *fatom = get_molecule().getAtomWithIdx(results.get_pairs()[i].first);
    auto *satom = get_molecule().getAtomWithIdx(results.get_pairs()[i].second);

    auto *finfo =
        static_cast<const RDKit::AtomPDBResidueInfo *>(fatom->getMonomerInfo());
    auto *sinfo =
        static_cast<const RDKit::AtomPDBResidueInfo *>(satom->getMonomerInfo());

    if (fatom->getAtomicNum() == 1 || satom->getAtomicNum() == 1)
      continue; // skip H atoms (hydrogens)

    auto f_resid = finfo->getResidueNumber();
    auto s_resid = sinfo->getResidueNumber();

    if (std::abs(f_resid - s_resid) > res_diff) {
      filtered.push_back(results.get_pairs()[i]);
      dists.push_back(results.get_distances()[i]);
    }
  }
  return NSResults(filtered, dists);
}

} // namespace lahuta
