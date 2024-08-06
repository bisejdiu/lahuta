#pragma once
#include <Geometry/point.h>
#include <array>
// #include <cmath>
#include <vector>

const int kDIMENSIONS = 3;
using FlatCoords = std::vector<float>;
using NeighborPairs = std::vector<std::pair<int, int>>;

void transform_coordinates(std::vector<RDGeom::Point3D> &coords,
                     std::array<float, 3> &pseudobox);
FlatCoords flatten_coordinates(std::vector<RDGeom::Point3D> &coords);

struct NSResults {
  std::vector<std::pair<int, int>> neighbor_pairs;
  std::vector<float> distances;

  void add_neighbors(int i, int j, float d2);

  void reserve_space(size_t input_size);

  const NeighborPairs &get_neighbors() const { return neighbor_pairs; }

  size_t size() const { return neighbor_pairs.size(); }

  NSResults filter(float distance) const;
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
