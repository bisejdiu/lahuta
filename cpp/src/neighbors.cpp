#include "neighbors.hpp"
#include "GraphMol/RWMol.h"
#include "lahuta.hpp"

namespace lahuta {

const RDKit::RWMol &ContextProvider<AtomAtomPair>::molecule() const {
  return luni->get_molecule();
}

const RDKit::RWMol &ContextProvider<AtomRingPair>::molecule() const {
  return luni->get_molecule();
}

const RingEntityCollection &ContextProvider<AtomRingPair>::rings() const {
  return luni->get_rings();
}

ContextProvider<AtomAtomPair>::ContextProvider(const Luni &ctx) : luni(&ctx){};
ContextProvider<AtomRingPair>::ContextProvider(const Luni &ctx) : luni(&ctx){};

std::string
AtomAtomPair::names(const ContextProvider<AtomAtomPair> &ctx) const {
  // FIX: Luni should support a member function: Luni::atom_name(int idx)
  auto atom1 = ctx.molecule().getAtomWithIdx(i);
  auto atom2 = ctx.molecule().getAtomWithIdx(j);
  auto res1 =
      static_cast<const RDKit::AtomPDBResidueInfo *>(atom1->getMonomerInfo());
  auto res2 =
      static_cast<const RDKit::AtomPDBResidueInfo *>(atom2->getMonomerInfo());
  return res1->getName() + " " + res2->getName();
}

std::string
AtomRingPair::names(const ContextProvider<AtomRingPair> &ctx) const {
  std::vector<int> atom_ids;
  auto ring = ctx.rings().get_data()[i];
  for (const auto &atom : ring.atoms) {
    atom_ids.push_back(atom->getIdx());
  }
  std::string atom_ids_str;
  for (const auto &id : atom_ids) {
    atom_ids_str += std::to_string(id) + " ";
  }

  auto res1 = static_cast<const RDKit::AtomPDBResidueInfo *>(ring.atoms.front()->getMonomerInfo());

  return res1->getResidueName() + " " + atom_ids_str;
}

// Neighbors class
template <typename T>
std::vector<T> Neighbors<T>::intersection(std::vector<T> data,
                                          std::vector<T> other) {
  std::vector<T> result;
  std::set_intersection(data.begin(), data.end(), other.begin(), other.end(),
                        std::back_inserter(result));
  return result;
}

template <typename T>
std::vector<T> Neighbors<T>::difference(std::vector<T> data,
                                        std::vector<T> other) {
  std::vector<T> result;
  std::set_difference(data.begin(), data.end(), other.begin(), other.end(),
                      std::back_inserter(result));
  return result;
}

template <typename T>
std::vector<T> Neighbors<T>::union_(std::vector<T> data, std::vector<T> other) {
  std::vector<T> result;
  std::set_union(data.begin(), data.end(), other.begin(), other.end(),
                 std::back_inserter(result));
  return result;
}

template <typename T>
std::vector<T> Neighbors<T>::symmetric_difference(std::vector<T> data,
                                                  std::vector<T> other) {
  std::vector<T> result;
  std::set_symmetric_difference(data.begin(), data.end(), other.begin(),
                                other.end(), std::back_inserter(result));
  return result;
}

template <typename T>
Neighbors<T> Neighbors<T>::intersection(const Neighbors<T> &other) const {
  std::vector<T> result = intersection(_data, other._data);
  return Neighbors<T>(*m_luni, result);
}

template <typename T>
Neighbors<T> Neighbors<T>::difference(const Neighbors<T> &other) const {
  std::vector<T> result = difference(_data, other._data);
  return Neighbors<T>(*m_luni, result);
}

template <typename T>
Neighbors<T>
Neighbors<T>::symmetric_difference(const Neighbors<T> &other) const {
  std::vector<T> result = symmetric_difference(_data, other._data);
  return Neighbors<T>(*m_luni, result);
}

template <typename T>
Neighbors<T> Neighbors<T>::union_(const Neighbors<T> &other) const {
  std::vector<T> result = union_(_data, other._data);
  return Neighbors<T>(*m_luni, result);
}

template <typename T>
Neighbors<T>
Neighbors<T>::filter(std::function<bool(const T &)> predicate) const {
  std::vector<T> result;
  std::copy_if(_data.begin(), _data.end(), std::back_inserter(result),
               predicate);
  return Neighbors<T>(*m_luni, result);
}

template <typename T> Pairs Neighbors<T>::get_pairs() const {
  std::vector<std::pair<int, int>> result;
  for (const T &p : _data) {
    result.push_back({p.i, p.j});
  }
  return result;
}

template <typename T> Distances Neighbors<T>::get_distances() const {
  std::vector<float> result;
  for (const T &p : _data) {
    result.push_back(p.d);
  }
  return result;
}

template <typename T> std::vector<std::string> Neighbors<T>::names() const {
  std::vector<std::string> result;
  for (const T &p : _data) {
    result.push_back(p.names(ctx));
  }
  return result;
}

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
  for (size_t i = 0; i < _data.size(); ++i) {
    if (partner == 0) {
      auto atype = m_luni->topology.atom_types[_data[i].i];
      if (AtomTypeFlags::has(atype, type)) {
        filtered.push_back(_data[i].get_pair());
        dists.push_back(_data[i].d);
      }
    } else if (partner == 1) {
      auto atype = m_luni->topology.atom_types[_data[i].j];
      if (AtomTypeFlags::has(atype, type)) {
        filtered.push_back(_data[i].get_pair());
        dists.push_back(_data[i].d);
      }
    } else {
      std::cerr << "Invalid partner: " << partner << std::endl;
    }
  }
  // FIX: does this cause problems on the python side?
  return Neighbors(*this->get_luni(), std::move(filtered), std::move(dists),
                   true);
}


// void add_neighbor(int i, int j, float d, bool sort) {
//   // data.push_back({i, j, d});
//   // if (sort) {
//   //   std::sort(data.begin(), data.end());
//   // }
//   T new_val(i, j, d);
//   auto it = std::lower_bound(data.begin(), data.end(), new_val);
//   data.insert(it, new_val);
// }

// explicit instantiation
template class Neighbors<AtomAtomPair>;
template class Neighbors<AtomRingPair>;

} // namespace lahuta
