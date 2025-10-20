#ifndef LAHUTA_NSRESULTS_HPP
#define LAHUTA_NSRESULTS_HPP

#include <cassert>
#include <vector>

#include <rdkit/Geometry/point.h>

namespace lahuta {

using Pairs = std::vector<std::pair<int, int>>;
using Distances = std::vector<float>;

struct NSResults {
  struct Iterator {
    using PairIterator = Pairs::const_iterator;
    using DistanceIterator = Distances::const_iterator;

    Iterator(PairIterator pair_it, DistanceIterator dist_it) : pair_it_(pair_it), dist_it_(dist_it) {}

    std::pair<std::pair<int, int>, float> operator*() const { return {*pair_it_, *dist_it_}; }

    Iterator &operator++() {
      ++pair_it_;
      ++dist_it_;
      return *this;
    }

    bool operator!=(const Iterator &other) const { return pair_it_ != other.pair_it_; }
    bool operator==(const Iterator &other) const { return pair_it_ == other.pair_it_; }

  private:
    PairIterator pair_it_;
    DistanceIterator dist_it_;
  };

  NSResults() = default;
  NSResults(const NSResults &other) = default;
  NSResults(NSResults &&other) = default;
  NSResults &operator=(const NSResults &other) = default;
  NSResults &operator=(NSResults &&other) = default;

  NSResults(Pairs &&pairs, std::vector<float> &&dists)
      : m_pairs(std::move(pairs)), m_dists(std::move(dists)) {}
  NSResults(Pairs &pairs, std::vector<float> &dists) : m_pairs(pairs), m_dists(dists) {}

  explicit NSResults(std::initializer_list<std::pair<int, int>> pairs, std::initializer_list<float> dists)
      : m_pairs(pairs.begin(), pairs.end()), m_dists(dists.begin(), dists.end()) {

    if (pairs.size() != dists.size()) {
      throw std::invalid_argument("Number of pairs must match number of distances");
    }
  }
  ~NSResults() {
#ifndef NDEBUG
    assert(m_pairs.size() == m_dists.size());
#endif
  }

  void add_neighbors(int i, int j, float d2);
  void reserve(size_t n) {
    m_pairs.reserve(n);
    m_dists.reserve(n);
  }
  void reserve_space(size_t input_size);

  //
  // These filter methods are old leftovers. They are used by Python tests and examples,
  // but not in the main C++ code. I don't find these methods particularly useful, but
  // I don't want to deal with these now. - Besian, October 2025
  //
  [[nodiscard]] NSResults filter(const double distance) const;
  [[nodiscard]] NSResults filter(const std::vector<int> &atom_indices) const;
  [[nodiscard]] NSResults filter(const std::vector<int> &atom_indices, int col) const;

  void clear() {
    m_pairs.clear();
    m_dists.clear();
  }
  size_t size() const { return m_pairs.size(); }

  const Pairs &get_pairs() const { return m_pairs; }
  const Distances &get_distances() const { return m_dists; }

  Iterator begin() const { return Iterator(m_pairs.begin(), m_dists.begin()); }
  Iterator end() const { return Iterator(m_pairs.end(), m_dists.end()); }

private:
  Pairs m_pairs;
  Distances m_dists;
};

} // namespace lahuta

#endif // LAHUTA_NSRESULTS_HPP
