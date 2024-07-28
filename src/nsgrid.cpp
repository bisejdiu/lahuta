#include <Geometry/point.h>
#include <algorithm>
#include <array>
// #include <cmath>
#include <vector>

#include "nsgrid.hpp"

static const int END = -1;
static const int MAX_GRID_DIM = 1290;

// FIX: To be used with bonding code
static const double EXTDIST = 0.45;
static const double MAXDIST = 5.45;     // 2 * MAXRAD + EXTDIST
static const double MINDIST2 = 0.16;    // MINDIST * MINDIST
static const double MAXDIST2 = 29.7025; // MAXDIST * MAXDIST

constexpr std::array<std::array<int, kDIMENSIONS>, 13> neighborCells = {
    {{1, 0, 0},
     {1, 1, 0},
     {0, 1, 0},
     {-1, 1, 0},
     {1, 0, -1},
     {1, 1, -1},
     {0, 1, -1},
     {-1, 1, -1},
     {1, 0, 1},
     {1, 1, 1},
     {0, 1, 1},
     {-1, 1, 1},
     {0, 0, 1}}};

FastNS::FastNS(const std::vector<RDGeom::Point3D> &coords, const float cutoff)
    : cutoff(cutoff) {

  std::cout << "Creating FastNS with cutoff: " << cutoff << std::endl;
  std::array<float, kDIMENSIONS> pbox = {0.0f, 0.0f, 0.0f};
  std::vector<RDGeom::Point3D> _coords(coords);
  transformCoords(_coords, pbox);

  if (pbox[0] <= 0.0f || pbox[1] <= 0.0f || pbox[2] <= 0.0f) {
    throw std::runtime_error(
        "Failed creating valid box from input coordinates");
  }

  if (pbox[0] < cutoff || pbox[1] < cutoff || pbox[2] < cutoff) {
    throw std::runtime_error(
        "Cutoff is larger than the smallest box dimension");
  }

  std::vector<float> flat_coords = flattenCoords(_coords);

  coords_bbox = std::move(flat_coords);
  box = std::move(pbox);

  buildGrid();
}

void FastNS::buildGrid() {
  prepareBox();
  packGrid();
}

NSResults FastNS::selfSearch() const {
  NSResults results;
  results.reserveSpace(coords_bbox.size() / 3);

  float cutoff2 = cutoff * cutoff;
  // FIX: use neighborCells
  std::array<std::array<int, 3>, 13> route = {{{1, 0, 0},
                                               {1, 1, 0},
                                               {0, 1, 0},
                                               {-1, 1, 0},
                                               {1, 0, -1},
                                               {1, 1, -1},
                                               {0, 1, -1},
                                               {-1, 1, -1},
                                               {1, 0, 1},
                                               {1, 1, 1},
                                               {0, 1, 1},
                                               {-1, 1, 1},
                                               {0, 0, 1}}};
  for (int cx = 0; cx < ncells[0]; ++cx) {
    for (int cy = 0; cy < ncells[1]; ++cy) {
      for (int cz = 0; cz < ncells[2]; ++cz) {
        int ci = _cellXYZ2CellId(cx, cy, cz);

        int i = head_id[ci];
        while (i != END) {
          int j = next_id[i];
          const float *coord_i = &coords_bbox[3 * i];
          while (j != END) {
            float d2 = calcDistSq(coord_i, &coords_bbox[3 * j]);
            if (d2 <= cutoff2) {
              // FIX: add another check for minimum distance
              results.addNeighbors(i, j, d2);
            }
            j = next_id[j];
          }

          for (const auto &offset : route) {
            int ox = cx + offset[0];
            int oy = cy + offset[1];
            int oz = cz + offset[2];

            int cj = _cellXYZ2CellId(ox, oy, oz);
            if (cj == END)
              continue;

            j = head_id[cj];
            while (j != END) {
              float d2 = calcDistSq(coord_i, &coords_bbox[3 * j]);
              if (d2 <= cutoff2) {
                results.addNeighbors(i, j, d2);
              }
              j = next_id[j];
            }
          }

          i = next_id[i];
        }
      }
    }
  }
  return results;
}

void FastNS::updateCutoff(float new_cutoff) {
  cutoff = new_cutoff;
  buildGrid();
};

void FastNS::prepareBox() {
  double min_cellsize = cutoff + 0.001;

  for (int i = 0; i < 3; ++i) {
    ncells[i] = std::min(static_cast<int>(std::floor(box[i] / min_cellsize)),
                         MAX_GRID_DIM);
  }

  cellsize[0] = box[0] / ncells[0];
  cellsize[1] = box[1] / ncells[1];
  cellsize[2] = box[2] / ncells[2];

  cell_offsets[0] = 0;
  cell_offsets[1] = ncells[0];
  cell_offsets[2] = ncells[0] * ncells[1];
};

void FastNS::packGrid() {
  head_id.assign(cell_offsets[2] * ncells[2], END);
  next_id.assign(coords_bbox.size() / 3, END);

  for (size_t i = 0; i < coords_bbox.size() / 3; ++i) {
    int j = _coord2CellId(&coords_bbox[3 * i]);
    next_id[i] = head_id[j];
    head_id[j] = i;
  }
};

inline int FastNS::_coord2CellId(const float *__restrict coord) const {
  std::array<int, kDIMENSIONS> xyz;
  _coord2CellXYZ(coord, xyz);
  return xyz[0] + xyz[1] * cell_offsets[1] + xyz[2] * cell_offsets[2];
}

inline void FastNS::_coord2CellXYZ(const float *__restrict coord,
                                   std::array<int, kDIMENSIONS> &xyz) const {
  xyz[2] = static_cast<int>(coord[2] / cellsize[2]); // % ncells[2];
  xyz[1] = static_cast<int>(coord[1] / cellsize[1]); // % ncells[1];
  xyz[0] = static_cast<int>(coord[0] / cellsize[0]); // % ncells[0];

  xyz[2] %= ncells[2];
  xyz[1] %= ncells[1];
  xyz[0] %= ncells[0];
}

inline int FastNS::_cellXYZ2CellId(int cx, int cy, int cz) const {
  if (cx < 0 || cx == ncells[0] || cy < 0 || cy == ncells[1] || cz < 0 ||
      cz == ncells[2]) {
    return END;
  }
  return cx + cy * cell_offsets[1] + cz * cell_offsets[2];
}

inline float FastNS::calcDistSq(const float *__restrict a,
                                const float *__restrict b) const {
  float dx = a[0] - b[0];
  float dy = a[1] - b[1];
  float dz = a[2] - b[2];
  return dx * dx + dy * dy + dz * dz;
}

inline bool FastNS::IsWithinCutoff(const float *__restrict a,
                                   const float *__restrict b,
                                   float cutoff2) const {
  float dx = a[0] - b[0];
  float dist2 = dx * dx;
  if (dist2 > cutoff2)
    return false;

  float dy = a[1] - b[1];
  dist2 += dy * dy;
  if (dist2 > cutoff2)
    return false;

  float dz = a[2] - b[2];
  dist2 += dz * dz;
  return dist2 <= cutoff2;
}

void NSResults::addNeighbors(int i, int j, float d2) {
  neighbor_pairs.emplace_back(i, j);
  distances.push_back(d2);
}

void NSResults::reserveSpace(size_t input_size) {
  neighbor_pairs.reserve(input_size);
  distances.reserve(input_size);
}

NSResults NSResults::filterByDistance(float dist) const {
  std::cout << "Filtering by distance: " << dist << std::endl;
  NSResults filtered; // we'll not reserve space 
  for (size_t i = 0; i < distances.size(); ++i) {
    if (distances[i] >= dist) {
      filtered.addNeighbors(neighbor_pairs[i].first, neighbor_pairs[i].second,
                            distances[i]);
    }
  }
  return filtered;
}

void transformCoords(std::vector<RDGeom::Point3D> &coords,
                     std::array<float, kDIMENSIONS> &pseudobox) {

  std::vector<double> lmax(3, std::numeric_limits<double>::lowest());
  std::vector<double> lmin(3, std::numeric_limits<double>::max());

  for (const auto &coord : coords) {
    for (int i = 0; i < 3; ++i) {
      lmax[i] = std::max(lmax[i], coord[i]);
      lmin[i] = std::min(lmin[i], coord[i]);
    }
  }

  // pseudobox
  for (int i = 0; i < 3; ++i) {
    pseudobox[i] = 1.1f * (lmax[i] - lmin[i]);
  }

  // shift coordinates
  for (auto &coord : coords) {
    for (int i = 0; i < 3; ++i) {
      coord[i] -= lmin[i];
    }
  }
}

std::vector<float> flattenCoords(std::vector<RDGeom::Point3D> &coords) {
  std::vector<float> flat_coords;
  flat_coords.reserve(coords.size() * 3);
  for (const auto &coord : coords) {
    flat_coords.insert(flat_coords.end(), {static_cast<float>(coord.x),
                                           static_cast<float>(coord.y),
                                           static_cast<float>(coord.z)});
  }
  return flat_coords;
}
