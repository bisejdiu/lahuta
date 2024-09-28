#ifndef LAHUTA_NEIGHBORS_HPP
#define LAHUTA_NEIGHBORS_HPP

#include <algorithm>
#include <vector>

#include <rdkit/GraphMol/RWMol.h>

#include "atom_types.hpp"
#include "nsgrid.hpp"
#include "rings.hpp"

namespace lahuta {

class Luni;

template<typename T>
class ContextProvider;

struct AtomRingPairType {
  const RDKit::RWMol *atom;
  const RingDataVec *ring;
};

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

    // FIXME: RefType, along with get_i and get_j are not being used. 
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
    return ref->atom->getAtomWithIdx(i);
  }

  const RingData* get_j(RefType *ref) const {
    return &(ref->ring->rings[j]);
  }

  std::string names(const ContextProvider<AtomRingPair> &ctx) const;
};

///////////////////////////////////////////////////////
///////////////// Neighbors class /////////////////////
///////////////////////////////////////////////////////

template<>
class ContextProvider<AtomAtomPair> {
public:
    explicit ContextProvider(const Luni& ctx);
    const RDKit::RWMol& molecule() const;

    const Luni* get_luni() const { return luni; }
private:
  const Luni* luni;
};

template<>
class ContextProvider<AtomRingPair> {
public:
    ContextProvider(const Luni& mainCtx);

    const RDKit::RWMol& molecule() const;
    const RingDataVec& rings() const;

    const Luni* get_luni() const { return luni; }

private:
  const Luni *luni;
};

template <typename T>
class Neighbors {

public:
  Neighbors(const Luni &luni, std::vector<T> data, bool is_sorted = false): 
    _data(std::move(data)), ctx(luni) {
    if (!is_sorted) {
      std::sort(this->_data.begin(), this->_data.end());
    }
  }

  /*Neighbors(const Luni &luni, Pairs pairs, Distances dists, bool is_sorted = false): */
  /*  contextProvider(luni) {*/
  /*  if (pairs.size() != dists.size()) {*/
  /*    throw std::runtime_error("Pairs and distances must have the same size");*/
  /*  }*/
  /*  for (size_t i = 0; i < pairs.size(); ++i) {*/
  /*    _data.push_back({pairs[i].first, pairs[i].second, dists[i]});*/
  /*  }*/
  /*  if (!is_sorted) {*/
  /*    std::sort(this->_data.begin(), this->_data.end());*/
  /*  }*/
  /*}*/
  // copy assignment
  Neighbors &operator=(const Neighbors &other) {
    if (this != &other) {
      _data = other._data;
    }
    return *this;
  }
  
  Neighbors(Neighbors &&other) noexcept: 
    ctx(other.ctx), m_luni(other.m_luni), _data(std::move(other._data)) {}

  Neighbors(const Luni &luni, const Pairs &&pairs, const Distances &&dists, bool is_sorted = false): 
    ctx(luni) {
    if (pairs.size() != dists.size()) {
      throw std::runtime_error("Pairs and distances must have the same size");
    }
    for (size_t i = 0; i < pairs.size(); ++i) {
      _data.push_back({pairs[i].first, pairs[i].second, dists[i]});
    }
    if (!is_sorted) {
      std::sort(this->_data.begin(), this->_data.end());
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
    std::vector<T> result = intersection(_data, other._data);
    return Neighbors<T>(*m_luni, result);
  }

  Neighbors<T> difference(const Neighbors<T> &other) const {
    std::vector<T> result = difference(_data, other._data);
    return Neighbors<T>(*m_luni, result);
  }

  Neighbors<T> symmetric_difference(const Neighbors<T> &other) const {
    std::vector<T> result = symmetric_difference(_data, other._data);
    return Neighbors<T>(*m_luni, result);
  }

  Neighbors<T> union_(const Neighbors<T> &other) const {
    std::vector<T> result = union_(_data, other._data);
    return Neighbors<T>(*m_luni, result);
  }

  Neighbors<T> filter(std::function<bool(const T &)> predicate) const {
    std::vector<T> result;
    std::copy_if(_data.begin(), _data.end(), std::back_inserter(result), predicate);
    return Neighbors<T>(*m_luni, result);
  }

  // NOTE: type_filter needs to be a supported member
  Neighbors<T> type_filter(AtomType type, int partner);

  // FIXME: Should FastNS be resonsible for defining these methods? 
  /*Neighbors<T> remove_adjascent_residueid_pairs(int res_diff);*/

  // void add_neighbor(int i, int j, float d, bool sort = true) {
  //   // data.push_back({i, j, d});
  //   // if (sort) {
  //   //   std::sort(data.begin(), data.end());
  //   // }
  //   T new_val(i, j, d);
  //   auto it = std::lower_bound(data.begin(), data.end(), new_val);
  //   data.insert(it, new_val);
  // }


  std::vector<std::string> names() const {
    std::vector<std::string> result;
    for (const T &p : _data) {
      result.push_back(p.names(ctx));
    }
    return result;
  }

  Pairs get_pairs() const { 
    std::vector<std::pair<int, int>> result;
    for (const T &p : _data) {
      result.push_back({p.i, p.j});
    }
    return result;
  }

  Distances get_distances() const {
    std::vector<float> result;
    for (const T &p : _data) {
      result.push_back(p.d);
    }
    return result;
  }

  size_t size() const { return _data.size(); }
  auto get_luni() const { return m_luni; }
  std::vector<T> data() const { return _data; }

  friend class Luni;

private:
  using RefType = typename T::RefType;
  const ContextProvider<T> ctx;
  const Luni* const m_luni = ctx.get_luni();
  std::vector<T> _data;
};

} // namespace lahuta

#endif // LAUTA_NEIGHBORS_HPP
