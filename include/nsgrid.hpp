#pragma once
#include <Geometry/point.h>
#include <array>
// #include <cmath>
#include <vector>

const int kDIMENSIONS = 3;
using FlatCoords = std::vector<float>;
using NeighborPairs = std::vector<std::pair<int, int>>;

void transformCoords(std::vector<RDGeom::Point3D> &coords,
                     std::array<float, 3> &pseudobox);
FlatCoords flattenCoords(std::vector<RDGeom::Point3D> &coords);

struct NSResults {
  std::vector<std::pair<int, int>> neighbor_pairs;
  std::vector<float> distances;

  void addNeighbors(int i, int j, float d2);

  void reserveSpace(size_t input_size);

  const NeighborPairs &getNeighbors() const { return neighbor_pairs; }

  size_t getNeighborPairsSize() const { return neighbor_pairs.size(); }
};

class FastNS {
public:
  FastNS(const RDGeom::POINT3D_VECT &coords, const float cutoff);

  NSResults selfSearch() const;

  void updateCutoff(float new_cutoff);

private:
  float cutoff;
  std::array<float, kDIMENSIONS> box;
  std::array<int, kDIMENSIONS> ncells;
  std::array<float, kDIMENSIONS> cellsize;
  std::array<int, kDIMENSIONS> cell_offsets;

  FlatCoords coords_bbox;
  std::vector<int> head_id;
  std::vector<int> next_id;

  void prepareBox();
  void packGrid();
  void buildGrid();

  inline int _coord2CellId(const float *__restrict coord) const;

  inline void _coord2CellXYZ(const float *__restrict coord,
                             std::array<int, 3> &xyz) const;

  inline int _cellXYZ2CellId(int cx, int cy, int cz) const;

  inline float calcDistSq(const float *__restrict a,
                          const float *__restrict b) const;
  inline bool IsWithinCutoff(const float *__restrict a,
                             const float *__restrict b, float cutoff2) const;
};
