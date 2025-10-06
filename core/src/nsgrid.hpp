#ifndef LAHUTA_NSGRID_HPP
#define LAHUTA_NSGRID_HPP

#include <array>
#include <limits>
#include <vector>
#include <rdkit/Geometry/point.h>

namespace lahuta {

// Inspired and adapted from MDAnalysis NSResults and FastNS Cython implementation.
// See: MDAnalysis/lib/nsgrid.pyx (6a75df0)

constexpr int DIMENSIONS = 3;
constexpr double MIN_VAL = std::numeric_limits<double>::lowest();
constexpr double MAX_VAL = std::numeric_limits<double>::max();

using Pairs     = std::vector<std::pair<int, int>>;
using Distances = std::vector<float>;

struct NSResults {
  struct Iterator {
    using PairIterator     = Pairs::const_iterator;
    using DistanceIterator = Distances::const_iterator;

    Iterator(PairIterator pair_it, DistanceIterator dist_it) : pair_it_(pair_it), dist_it_(dist_it) {}

    std::pair<std::pair<int, int>, float> operator*() const {
      return {*pair_it_, *dist_it_};
    }

    Iterator &operator++() {
      ++pair_it_;
      ++dist_it_;
      return *this;
    }

    bool operator!=(const Iterator &other) const {return pair_it_ != other.pair_it_;}
    bool operator==(const Iterator &other) const {return pair_it_ == other.pair_it_;}

  private:
    PairIterator     pair_it_;
    DistanceIterator dist_it_;
  };

  NSResults() = default;
  NSResults(const NSResults &other) = default;
  NSResults(NSResults &&other) = default;
  NSResults &operator=(const NSResults &other) = default;
  NSResults &operator=(NSResults &&other) = default;

  NSResults(Pairs &&pairs, std::vector<float> &&dists): m_pairs(std::move(pairs)), m_dists(std::move(dists)) {}
  NSResults(Pairs  &pairs, std::vector<float>  &dists): m_pairs(pairs),            m_dists(dists) {}

  explicit NSResults(std::initializer_list<std::pair<int, int>> pairs, std::initializer_list<float> dists)
      : m_pairs(pairs.begin(), pairs.end()), m_dists(dists.begin(), dists.end()) {

    if (pairs.size() != dists.size()) {
      throw std::invalid_argument("Number of pairs must match number of distances");
    }
  }
  ~NSResults() {}

  void add_neighbors(int i, int j, float d2);
  void reserve_space(size_t input_size); // FIX: not accurate. should use a smarter estimate

  [[nodiscard]] NSResults filter(const double distance) const;
  [[nodiscard]] NSResults filter(const std::vector<int> &atom_indices) const;
  [[nodiscard]] NSResults filter(const std::vector<int> &atom_indices, int col) const;

  void clear() {
    m_pairs.clear();
    m_dists.clear();
  }
  size_t size() const { return m_pairs.size(); }

  const Pairs     &get_pairs()     const { return m_pairs; }
  const Distances &get_distances() const { return m_dists; }

  Iterator begin() const { return Iterator(m_pairs.begin(), m_dists.begin());}
  Iterator end()   const { return Iterator(m_pairs.end(),   m_dists.end());}

private:
  Pairs     m_pairs;
  Distances m_dists;
};


class FastNS {
public:
  FastNS() = default; // used by Luni, for performance reasons, to store the grid
  FastNS(const RDGeom::POINT3D_VECT &coords);
  FastNS(const std::vector<std::vector<double>> &coords);
  FastNS(const double *coords_ptr, std::size_t n_points);
  ~FastNS() = default;

  bool build(double cutoff, bool brute_force_fallback = true);

  // find all neighbors among the input coordinates within the cutoff distance
  NSResults self_search() const;

  // find all neighbors between provided search coords and the input coordinates within the cutoff distance
  NSResults search(const RDGeom::POINT3D_VECT &search_coords) const;
  NSResults search(const std::vector<std::vector<double>> &search_coords) const;

  double get_cutoff() const { return cutoff; }

  // accessors useful for debugging and testing
  std::array<int,   DIMENSIONS> get_ncells()   const { return ncells; }
  std::array<float, DIMENSIONS> get_cellsize() const { return cellsize; }
  std::array<float, DIMENSIONS> get_box()      const { return box; }
  bool is_grid_ready() const { return grid_ready; }

  template <typename T>
  static inline T dist_sq(const T* __restrict a, const T* __restrict b) {
      T dx = a[0] - b[0];
      T dy = a[1] - b[1];
      T dz = a[2] - b[2];
      return dx * dx + dy * dy + dz * dz;
  }

private:
  bool configure_grid(const std::array<float, DIMENSIONS> &extents, float padding, float min_cell);
  void reset_state();
  void pack_grid();
  void build_grid();
  void log_occupancy_stats() const;

  NSResults brute_force_self(double cutoff_sq) const;

  inline int  coord_to_cell_id (const float *__restrict coord) const;
  inline void coord_to_cell_xyz(const float *__restrict coord, std::array<int, 3> &xyz) const;
  inline int  cell_xyz_to_cell_id(int cx, int cy, int cz) const;

  // Note: only internal bounds (_lmin/_lmax) are kept for packing and grid setup

  std::array<float, DIMENSIONS> box          = {0.0f, 0.0f, 0.0f};
  std::array<float, DIMENSIONS> cellsize     = {0.0f, 0.0f, 0.0f};
  std::array<int,   DIMENSIONS> ncells       = {0, 0, 0};
  std::array<int,   DIMENSIONS> cell_offsets = {0, 0, 0};

  std::vector<double> _lmin;
  std::vector<double> _lmax;
  std::vector<int>    head_id;
  std::vector<int>    next_id;

  std::vector<float>  coords_bbox;
  std::size_t n_points = 0;

  double cutoff = 0.0;
  bool grid_ready = false;
  bool brute_force_mode = false;
  NSResults brute_force_results;

};

class Luni;
namespace ns_utils {

NSResults remove_adjascent_residueid_pairs(const Luni &luni, NSResults &results, int res_diff);

} // namespace ns_utils

} // namespace lahuta

#endif // LAHUTA_NSGRID_HPP
