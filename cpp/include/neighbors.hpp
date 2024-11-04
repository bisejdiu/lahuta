#ifndef LAHUTA_NEIGHBORS_HPP
#define LAHUTA_NEIGHBORS_HPP

#include <algorithm>
#include <vector>

#include <rdkit/GraphMol/RWMol.h>

#include "atom_types.hpp"
#include "nsgrid.hpp"
#include "entities.hpp"

namespace lahuta {

class Luni;

template <typename T> class ContextProvider;

struct AtomRingPairType {
  const RDKit::RWMol *mol;
  const RingEntityCollection *ring;
};

template <typename T> class BasePair {
public:
  int i;
  int j;
  float d;

  using RefType = const T;

  BasePair(int i, int j, float d, bool sort = true) : d(d) {
    if (sort) {
      std::tie(this->i, this->j) = std::minmax(i, j);
    } else {
      this->i = i;
      this->j = j;
    }
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
  const RDKit::Atom *get_i(RefType *mol) const {
    return mol->getAtomWithIdx(i);
  }

  const RDKit::Atom *get_j(RefType *mol) const {
    return mol->getAtomWithIdx(j);
  }

  std::string names(const ContextProvider<AtomAtomPair> &ctx) const;
};

class AtomRingPair : public BasePair<AtomRingPairType> {
public:
  // FIXME:
  /*using BasePair::BasePair;*/
  AtomRingPair(int i, int j, float d, bool sort = false)
      : BasePair(i, j, d, false) {}
  /*AtomRingPair(int i, int j, float d) : BasePair(i, j, d, false) {}*/

  /*const RDKit::Atom *get_i(RefType *ref) const {*/
  /*  return ref->mol->getAtomWithIdx(i);*/
  /*}*/

  /*const RingData *get_j(RefType *ref) const { return &(ref->ring->rings[j]);
   * }*/

  std::string names(const ContextProvider<AtomRingPair> &ctx) const;
};

///////////////////////////////////////////////////////
///////////////// Neighbors class /////////////////////
///////////////////////////////////////////////////////

template <> class ContextProvider<AtomAtomPair> {
public:
  explicit ContextProvider(const Luni &ctx);
  const RDKit::RWMol &molecule() const;

  const Luni *get_luni() const { return luni; }

private:
  const Luni *luni;
};

template <> class ContextProvider<AtomRingPair> {
public:
  ContextProvider(const Luni &mainCtx);

  const RDKit::RWMol &molecule() const;
  const RingEntityCollection &rings() const;

  const Luni *get_luni() const { return luni; }

private:
  const Luni *luni;
};

template <typename T> class Neighbors {

public:
  Neighbors() = delete;
  Neighbors(const Luni &luni) : ctx(luni), m_luni(&luni) {}
  Neighbors(const Neighbors &other)
      : ctx(other.ctx), m_luni(other.m_luni), _data(other._data) {}

  Neighbors(Neighbors &&other) noexcept
      : ctx(other.ctx), m_luni(other.m_luni), _data(std::move(other._data)) {}

  Neighbors &operator=(const Neighbors &other) {
    if (this != &other) {
      _data = other._data;
    }
    return *this;
  }

  Neighbors &operator=(Neighbors &&other) noexcept {
    if (this != &other) {
      _data = std::move(other._data);
    }
    return *this;
  }

  Neighbors(const Luni &luni, std::vector<T> data, bool is_sorted = false)
      : _data(std::move(data)), ctx(luni) {
    if (!is_sorted) {
      std::sort(this->_data.begin(), this->_data.end());
    }
  }

  Neighbors(const Luni &luni, const Pairs &&pairs, const Distances &&dists,
            bool is_sorted = false)
      : ctx(luni) {
    if (pairs.size() != dists.size()) {
      throw std::runtime_error("Pairs and distances must have the same size");
    }
    _data.reserve(pairs.size());
    for (size_t i = 0; i < pairs.size(); ++i) {
      _data.push_back({std::move(pairs[i].first), std::move(pairs[i].second),
                       std::move(dists[i]), !is_sorted});
    }
    if (!is_sorted) {
      std::sort(this->_data.begin(), this->_data.end());
    }
  }

  Neighbors(const Luni &luni, NSResults &&results, bool is_sorted = false)
      : ctx(luni) {
    if (results.get_pairs().size() != results.get_distances().size()) {
      throw std::runtime_error("Pairs and distances must have the same size");
    }
    _data.reserve(results.get_pairs().size());
    for (size_t i = 0; i < results.get_pairs().size(); ++i) {
      _data.push_back({std::move(results.get_pairs()[i].first),
                       std::move(results.get_pairs()[i].second),
                       std::move(results.get_distances()[i])});
    }
    if (!is_sorted) {
      std::sort(this->_data.begin(), this->_data.end());
    }
  }

  // FIX: Does the copy-based logic make a significant performance difference?
  static std::vector<T> intersection(std::vector<T> data, std::vector<T> other);
  static std::vector<T> difference(std::vector<T> data, std::vector<T> other);
  static std::vector<T> union_(std::vector<T> data, std::vector<T> other);
  static std::vector<T> symmetric_difference(std::vector<T> data,
                                             std::vector<T> other);

  Neighbors<T> intersection(const Neighbors<T> &other) const;
  Neighbors<T> difference(const Neighbors<T> &other) const;
  Neighbors<T> symmetric_difference(const Neighbors<T> &other) const;
  Neighbors<T> union_(const Neighbors<T> &other) const;

  Neighbors<T> operator+(const Neighbors<T> &other) const {
    return union_(other);
  }
  Neighbors<T> operator-(const Neighbors<T> &other) const {
    return difference(other);
  }
  Neighbors<T> operator^(const Neighbors<T> &other) const {
    return symmetric_difference(other);
  }
  Neighbors<T> operator&(const Neighbors<T> &other) const {
    return intersection(other);
  }
  /*Neighbors<T> operator|(const Neighbors<T> &other) const {*/
  /*  return union_(other);*/
  /*}*/

  Neighbors<T> filter(std::function<bool(const T &)> predicate) const;
  Neighbors<T> type_filter(AtomType type, int partner);

  // void add_neighbor(int i, int j, float d, bool sort = true);

  std::vector<std::string> names() const;
  Pairs get_pairs() const;
  Distances get_distances() const;

  size_t size() const { return _data.size(); }
  auto get_luni() const { return m_luni; }
  std::vector<T> get_data() const { return _data; }

  friend class Luni;

private:
  using RefType = typename T::RefType;
  const ContextProvider<T> ctx;
  const Luni *const m_luni = ctx.get_luni();
  std::vector<T> _data;
};

} // namespace lahuta

#endif // LAUTA_NEIGHBORS_HPP
