#ifndef LAHUTA_NSGRID_HPP
#define LAHUTA_NSGRID_HPP

#include <array>
#include <vector>
#include <rdkit/Geometry/point.h>
#include "atom_types.hpp"

namespace lahuta {

const int kDIMENSIONS = 3;
using FlatCoords = std::vector<float>;
using _NeighborPairs = std::vector<std::pair<int, int>>;

void transform_coordinates(std::vector<RDGeom::Point3D> &coords,
                           std::array<float, 3> &pseudobox);
FlatCoords flatten_coordinates(std::vector<RDGeom::Point3D> &coords);

class Luni;

struct NSResults {
  _NeighborPairs m_pairs;
  std::vector<float> m_distances;

public:
  NSResults() = default;
  NSResults(const NSResults &other) = default;
  NSResults(NSResults &&other) = default;
  NSResults &operator=(const NSResults &other) = default;
  NSResults &operator=(NSResults &&other) = default;

  NSResults(_NeighborPairs &&pairs, std::vector<float> &&dists)
      : m_pairs(std::move(pairs)), m_distances(std::move(dists)) {}

  NSResults(_NeighborPairs &pairs, std::vector<float> &dists)
      : m_pairs(pairs), m_distances(dists) {}

  NSResults(Luni &luni, _NeighborPairs &&pairs, std::vector<float> &&dists)
      : m_pairs(std::move(pairs)), m_distances(std::move(dists)) {
    this->m_luni = &luni;
  }

  explicit NSResults(Luni &luni, _NeighborPairs &pairs,
                     std::vector<float> &dists)
      : m_pairs(pairs), m_distances(dists) {
    this->m_luni = &luni;
  }

  // NSResults results = {{1, 2}, {3, 4}, {5, 6}}, {0.1f, 0.2f, 0.3f}};
  NSResults(std::initializer_list<std::pair<int, int>> pairs,
            std::initializer_list<float> dists)
      : m_pairs(pairs.begin(), pairs.end()),
        m_distances(dists.begin(), dists.end()) {
    if (pairs.size() != dists.size()) {
      throw std::invalid_argument(
          "Number of pairs must match number of distances");
    }
  }

  [[nodiscard]] Luni *get_luni() const { return m_luni; }

  void add_neighbors(int i, int j, float d2);

  void reserve_space(size_t input_size);

  const _NeighborPairs &get_neighbors() const { return m_pairs; }
  const std::vector<float> &get_distances() const { return m_distances; }

  size_t size() const { return m_pairs.size(); }

  [[nodiscard]] NSResults filter(float distance) const;

  void clear() {
    m_pairs.clear();
    m_distances.clear();
  }

  NSResults filter_by_atom_type(AtomType type, int partner);

  friend class Luni;

private:
  Luni *m_luni = nullptr;
  // std::vector<AtomType> &atom_types = m_luni->atom_types;
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

  FlatCoords coords_bbox;
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
