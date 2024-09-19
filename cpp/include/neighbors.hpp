#ifndef LAHUTA_NEIGHBORS_HPP
#define LAHUTA_NEIGHBORS_HPP

#include <algorithm>
#include <vector>

#include "atom_types.hpp"
#include "nsgrid.hpp"
#include "rings.hpp"
#include <rdkit/GraphMol/RWMol.h>

namespace lahuta {

class Luni;

template<typename T>
class ContextProvider;

struct AtomRingPairType {
  const RDKit::RWMol *i;
  const RingDataVec *j;
};

enum class ContactType { AtomAtom, AtomRing };


template <typename T>
class BasePair {
public:
    int i;
    int j;
    float d;

    using RefType = const T;

    BasePair(int i, int j, float d) : d(d) {
      std::tie(this->i, this->j) = std::minmax(i, j);
    }

    bool operator==(const BasePair &other) const {
        return (i == other.i) && (j == other.j);
    }

    bool operator!=(const BasePair &other) const {
        return i != other.i || j != other.j;
    }

    bool operator<(const BasePair &other) const {
        return i != other.i ? i < other.i : j < other.j;
    }

    bool operator>(const BasePair &other) const {
        return i != other.i ? i > other.i : j > other.j;
    }

    bool operator<=(const BasePair &other) const {
        return i != other.i ? i <= other.i : j <= other.j;
    }

    bool operator>=(const BasePair &other) const {
        return i != other.i ? i >= other.i : j >= other.j;
    }

    std::pair<int, int> get_pair() const { return {i, j}; }
    float get_distance() const { return d; }
};

class AtomAtomPair : public BasePair<RDKit::RWMol> {
public:
    using BasePair::BasePair;

    const RDKit::Atom* get_i(RefType *mol) const {
        return mol->getAtomWithIdx(i);
    }

    const RDKit::Atom* get_j(RefType *mol) const {
        return mol->getAtomWithIdx(j);
    }

    std::string names(const ContextProvider<AtomAtomPair> &ctx) const;
};

class AtomRingPair : public BasePair<AtomRingPairType> {
public:
  using BasePair::BasePair;

  const RDKit::Atom* get_i(RefType *ref) const {
    return ref->i->getAtomWithIdx(i);
  }

  const RingData* get_j(RefType *ref) const {
    return &(ref->j->rings[j]);
    // return &(ref->j->rings.at(j)); // FIXME: use [] instead of at
  }

  // FIXME: not yet implemented. 
  std::string names(const ContextProvider<AtomRingPair> &ctx) const;
  // std::string names(RefType *ref) const {
  //   std::vector<int> atom_ids;
  //   auto ring = ref->j->rings.at(j);
  //   for (const auto &atom : ring.atom_ids) {
  //     atom_ids.push_back(atom);
  //   }
  //   std::string atom_ids_str;
  //   for (const auto &id : atom_ids) {
  //     atom_ids_str += std::to_string(id) + " ";
  //   }
  //
  //   return "->" + get_i(ref)->getSymbol() + " " + atom_ids_str;
  // }
};

///////////////////////////////////////////////////////
///////////////// Neighbors class /////////////////////
///////////////////////////////////////////////////////

template<>
class ContextProvider<AtomAtomPair> {
public:
    const RDKit::RWMol *molecule;
    ContextProvider(const Luni& mainCtx);
};

template<>
class ContextProvider<AtomRingPair> {
public:
    const RDKit::RWMol *molecule;
    const RingDataVec *rings;

    ContextProvider(const Luni& mainCtx);
};

template <typename T>
class Neighbors {

public:
  using RefType = typename T::RefType;

  std::vector<T> data;
  const ContextProvider<T> contextProvider;
  
public:
  Neighbors(const Luni &luni, std::vector<T> data, bool is_sorted = false): 
    m_luni(&luni), data(std::move(data)), contextProvider(luni) {
    if (!is_sorted) {
      std::sort(this->data.begin(), this->data.end());
    }
  }

  Neighbors(const Luni &luni, Pairs pairs, Distances dists, bool is_sorted = false): 
    m_luni(&luni), contextProvider(luni) {
    if (pairs.size() != dists.size()) {
      throw std::runtime_error("Pairs and distances must have the same size");
    }
    for (size_t i = 0; i < pairs.size(); ++i) {
      data.push_back({pairs[i].first, pairs[i].second, dists[i]});
    }
    if (!is_sorted) {
      std::sort(this->data.begin(), this->data.end());
    }
  }

  // FIXME: here I am making a copy
  static inline std::vector<T> intersection(std::vector<T> data, std::vector<T> other) {
    std::vector<T> result;
    // std::sort(data.begin(), data.end());
    // std::sort(other.begin(), other.end());

    std::set_intersection(data.begin(), data.end(),
      other.begin(), other.end(),
      std::back_inserter(result) 
    );
    return result;
  }

  static inline std::vector<T> difference(std::vector<T> data, std::vector<T> other) {
    std::vector<T> result;
    std::set_difference(data.begin(), data.end(),
      other.begin(), other.end(),
      std::back_inserter(result)
    );
    return result;
  }

  static inline std::vector<T> union_(std::vector<T> data, std::vector<T> other) {
    std::vector<T> result;
    std::set_union(data.begin(), data.end(),
      other.begin(), other.end(),
      std::back_inserter(result)
    );
    return result;
  }

  static inline std::vector<T> symmetric_difference(std::vector<T> data, std::vector<T> other) {
    std::vector<T> result;
    std::set_symmetric_difference(data.begin(), data.end(),
      other.begin(), other.end(),
      std::back_inserter(result)
    );
    return result;
  }

  Neighbors<T> intersection(const Neighbors<T> &other) const {
    std::vector<T> result = intersection(data, other.data);
    return Neighbors<T>(*m_luni, result);
  }

  Neighbors<T> difference(const Neighbors<T> &other) const {
    std::vector<T> result = difference(data, other.data);
    return Neighbors<T>(*m_luni, result);
  }

  Neighbors<T> symmetric_difference(const Neighbors<T> &other) const {
    std::vector<T> result = symmetric_difference(data, other.data);
    return Neighbors<T>(*m_luni, result);
  }

  Neighbors<T> union_(const Neighbors<T> &other) const {
    std::vector<T> result = union_(data, other.data);
    return Neighbors<T>(*m_luni, result);
  }

  Neighbors<T> filter(std::function<bool(const T &)> predicate) const {
    std::vector<T> result;
    std::copy_if(data.begin(), data.end(), std::back_inserter(result), predicate);
    return Neighbors<T>(*m_luni, result);
  }

  // FIXME: Should FastNS be resonsible for defining these methods? 
  Neighbors<T> type_filter(AtomType type, int partner);
  Neighbors<T> remove_adjascent_residueid_pairs(int res_diff);

  // void add_neighbor(int i, int j, float d, bool sort = true) {
  //   // data.push_back({i, j, d});
  //   // if (sort) {
  //   //   std::sort(data.begin(), data.end());
  //   // }
  //   T new_val(i, j, d);
  //   auto it = std::lower_bound(data.begin(), data.end(), new_val);
  //   data.insert(it, new_val);
  //
  // }


  std::vector<std::string> names() const {
    std::vector<std::string> result;
    for (const T &p : data) {
      result.push_back(p.names(contextProvider));
    }
    return result;
  }

  Pairs get_pairs() const { 
    std::vector<std::pair<int, int>> result;
    for (const T &p : data) {
      result.push_back({p.i, p.j});
    }
    return result;
  }

  Distances get_distances() const {
    std::vector<float> result;
    for (const T &p : data) {
      result.push_back(p.d);
    }
    return result;
  }

  size_t size() const { return data.size(); }
  auto get_luni() const { return m_luni; }

  friend class Luni;

public:
  const Luni* const m_luni = nullptr;
  // Luni &m_luni;
  // std::unique_ptr<Luni> m_luni = nullptr;

};

} // namespace lahuta

#endif // LAUTA_NEIGHBORS_HPP
