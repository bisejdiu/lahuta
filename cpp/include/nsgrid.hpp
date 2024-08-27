#ifndef LAHUTA_NSGRID_HPP
#define LAHUTA_NSGRID_HPP

#include "atom_types.hpp"
#include <array>
#include <rdkit/Geometry/point.h>
#include <vector>

namespace lahuta {

const int kDIMENSIONS = 3;
using Pairs = std::vector<std::pair<int, int>>;
using Distances = std::vector<float>;

void transform_coordinates(std::vector<RDGeom::Point3D> &coords,
                           std::array<float, 3> &pseudobox);
std::vector<float> flatten_coordinates(std::vector<RDGeom::Point3D> &coords);

class Luni;

struct NSResults {

  NSResults() = default;
  NSResults(const NSResults &other) = default;
  NSResults(NSResults &&other) = default;
  NSResults &operator=(const NSResults &other) = default;
  NSResults &operator=(NSResults &&other) = default;

  NSResults(Pairs &&pairs, std::vector<float> &&dists)
      : m_pairs(std::move(pairs)), m_dists(std::move(dists)) {}

  NSResults(Pairs &pairs, std::vector<float> &dists)
      : m_pairs(pairs), m_dists(dists) {}

  NSResults(Luni &luni, Pairs &&pairs, std::vector<float> &&dists)
      : m_pairs(std::move(pairs)), m_dists(std::move(dists)) {
    this->m_luni = &luni;
  }

  NSResults(Luni &luni, Pairs &pairs, std::vector<float> &dists)
      : m_pairs(pairs), m_dists(dists) {
    this->m_luni = &luni;
  }

  // NSResults results = {{1, 2}, {3, 4}, {5, 6}}, {0.1f, 0.2f, 0.3f}};
  explicit NSResults(std::initializer_list<std::pair<int, int>> pairs,
                     std::initializer_list<float> dists)
      : m_pairs(pairs.begin(), pairs.end()),
        m_dists(dists.begin(), dists.end()) {
    if (pairs.size() != dists.size()) {
      throw std::invalid_argument(
          "Number of pairs must match number of distances");
    }
  }

  [[nodiscard]] Luni *get_luni() const { return m_luni; }

  void add_neighbors(int i, int j, float d2);

  void reserve_space(size_t input_size);

  size_t size() const { return m_pairs.size(); }

  [[nodiscard]] NSResults filter(float distance) const;

  void clear() {
    m_pairs.clear();
    m_dists.clear();
  }

  NSResults type_filter(AtomType type, int partner);

  const Pairs &get_pairs() const { return m_pairs; }
  const Distances &get_distances() const { return m_dists; }

  auto begin() { return m_pairs.begin(); }
  auto end() { return m_pairs.end(); }

  friend class Luni;

private:
  Pairs m_pairs;
  Distances m_dists;
  Luni *m_luni = nullptr;
};

class FastNS {
public:
  FastNS() = default;
  FastNS(const RDGeom::POINT3D_VECT &coords, float cutoff);

  NSResults self_search() const;

  void update_cutoff(float new_cutoff);

private:
  float cutoff;
  std::array<float, kDIMENSIONS> box;
  std::array<int, kDIMENSIONS> ncells;
  std::array<float, kDIMENSIONS> cellsize;
  std::array<int, kDIMENSIONS> cell_offsets;

  std::vector<float> coords_bbox;
  std::vector<int> head_id;
  std::vector<int> next_id;

  void prepare_box();
  void pack_grid();
  void build_grid();

  inline int _coord_to_cell_id(const float *__restrict coord) const;

  inline void _coord_to_cell_xyz(const float *__restrict coord,
                                 std::array<int, 3> &xyz) const;

  inline int _cell_xyz_to_cell_id(int cx, int cy, int cz) const;

  inline float dist_sq(const float *__restrict a,
                       const float *__restrict b) const;
  inline bool is_within_cutoff(const float *__restrict a,
                               const float *__restrict b, float cutoff2) const;
};

} // namespace lahuta

#endif // LAHUTA_NSGRID_HPP
