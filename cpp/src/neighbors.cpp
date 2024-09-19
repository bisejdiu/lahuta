#include "GraphMol/RWMol.h"
#include "gemmi/cif.hpp"
#include "lahuta.hpp"
#include "neighbors.hpp"

namespace lahuta {


ContextProvider<AtomAtomPair>::ContextProvider(const Luni& mainCtx)
    : molecule(&mainCtx.get_molecule()) {}

ContextProvider<AtomRingPair>::ContextProvider(const Luni& mainCtx)
        : molecule(&mainCtx.get_molecule()), rings(&mainCtx.get_rings()) {}

std::string AtomAtomPair::names(const ContextProvider<AtomAtomPair> &ctx) const {
    return ctx.molecule->getAtomWithIdx(i)->getSymbol() + " " +
           ctx.molecule->getAtomWithIdx(j)->getSymbol();
}


std::string AtomRingPair::names(const ContextProvider<AtomRingPair> &ctx) const {
    std::vector<int> atom_ids;
    auto ring = ctx.rings->rings[j];
    for (const auto &atom : ring.atom_ids) {
      atom_ids.push_back(atom);
    }
    std::string atom_ids_str;
    for (const auto &id : atom_ids) {
      atom_ids_str += std::to_string(id) + " ";
    }

    return "->" + ctx.molecule->getAtomWithIdx(i)->getSymbol() + " " + atom_ids_str;
  }


template class Neighbors<AtomAtomPair>;
template class Neighbors<AtomRingPair>;

// template <typename T>
// Neighbors<T>::Neighbors(const Luni &luni, Pairs pairs, Distances dists, ContactType c, bool is_sorted):
//     m_luni(&luni) {
//     if (pairs.size() != dists.size()) {
//       throw std::runtime_error("Pairs and distances must have the same size");
//     }
//     using RefType = typename T::RefType;
//     // const RefType *ctx = nullptr;
//     // switch (c) {
//     //   case ContactType::AtomAtom: {
//     //     const RDKit::RWMol* mol = &luni.get_molecule();
//     //     ctx = mol; 
//     //     break;
//     //   }
//     //   case ContactType::AtomRing:
//     //   const RDKit::RWMol* mol = &luni.get_molecule();
//     //   AtomRingPairType ap = {mol, &luni.get_rings()};
//     //     ctx = ap;
//     //     break;
//     // }
//
//     for (size_t i = 0; i < pairs.size(); ++i) {
//       data.push_back({pairs[i].first, pairs[i].second, dists[i]});
//     }
//     if (!is_sorted) {
//       std::sort(this->data.begin(), this->data.end());
//     }
//     
//     
//   }


template <typename T>
Neighbors<T> Neighbors<T>::type_filter(AtomType type, int partner) {
  if (partner != 0 && partner != 1) {
  throw std::runtime_error("Invalid partner: " + std::to_string(partner) +
                           ". Must be 0 or 1.");
  }
  if (m_luni == nullptr) {
    throw std::runtime_error("Luni is not set");
  }

  Pairs filtered;
  Distances dists;
  for (size_t i = 0; i < data.size(); ++i) {
    if (partner == 0) {
      auto atype = m_luni->topology.atom_types[data[i].i];
      if (AtomTypeFlags::has(atype, type)) {
        filtered.push_back(data[i].get_pair());
        dists.push_back(data[i].d);
      }
    } else if (partner == 1) {
      auto atype = m_luni->topology.atom_types[data[i].j];
      if (AtomTypeFlags::has(atype, type)) {
        filtered.push_back(data[i].get_pair());
        dists.push_back(data[i].d);
      }
    } else {
      std::cerr << "Invalid partner: " << partner << std::endl;
    }
  }
  return Neighbors(*this->get_luni(), filtered, dists, false);
}

template <typename T>
Neighbors<T>  Neighbors<T>::remove_adjascent_residueid_pairs(int res_diff) {
  Pairs filtered;
  Distances dists;
  for (size_t i = 0; i < data.size(); ++i) {
    auto *fatom = m_luni->get_molecule().getAtomWithIdx(data[i].i);
    auto *finfo =
        static_cast<const RDKit::AtomPDBResidueInfo *>(fatom->getMonomerInfo());
    auto *satom = m_luni->get_molecule().getAtomWithIdx(data[i].j);
    auto *sinfo =
        static_cast<const RDKit::AtomPDBResidueInfo *>(satom->getMonomerInfo());

    if (fatom->getAtomicNum() == 1 || satom->getAtomicNum() == 1)
      continue; // skip H atoms (hydrogens)

    auto f_resid = finfo->getResidueNumber();
    auto s_resid = sinfo->getResidueNumber();

    if (std::abs(f_resid - s_resid) > res_diff) {
      filtered.push_back(data[i].get_pair());
      dists.push_back(data[i].d);
    }
  }
  return Neighbors(*this->get_luni(), filtered, dists, false);
}

} // namespace lahuta
