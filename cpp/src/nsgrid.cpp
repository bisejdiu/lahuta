#include <algorithm>
#include <unordered_set>

#include <rdkit/Geometry/point.h>

#include "nsgrid.hpp"
#include "spdlog/spdlog.h"

namespace lahuta {

static const int END = -1;
static const int MAX_GRID_DIM = 1290;
static const int SMALL_SYSTEM_THRESHOLD = 10'000; // upper cutoff for brute-force search

constexpr std::array<std::array<int, DIMENSIONS>, 13> NeighborCells = {{
  {1, 0,  0}, {1, 1,  0}, {0, 1,  0}, {-1, 1,  0},
  {1, 0, -1}, {1, 1, -1}, {0, 1, -1}, {-1, 1, -1},
  {1, 0,  1}, {1, 1,  1}, {0, 1,  1}, {-1, 1,  1},
  {0, 0,  1}
}};

FastNS::FastNS(const RDGeom::POINT3D_VECT &coords, float scale_factor) : _coords(coords), _scale_factor(scale_factor) {

  // initialize min/max bounds
  _lmin.resize(DIMENSIONS, MAX_VAL);
  _lmax.resize(DIMENSIONS, MIN_VAL);

  // calculate bounds
  for (const auto &coord : _coords) {
    for (int i = 0; i < DIMENSIONS; ++i) {
      _lmax[i] = std::max(_lmax[i], coord[i]);
      _lmin[i] = std::min(_lmin[i], coord[i]);
    }
  }
}

FastNS::FastNS(const std::vector<std::vector<double>> &coords, float scale_factor) : _scale_factor(scale_factor) {

  _coords.reserve(coords.size());
  for (const auto &coord : coords) {
    _coords.emplace_back(coord[0], coord[1], coord[2]);
  }

  // initialize min/max bounds
  _lmin.resize(DIMENSIONS, MAX_VAL);
  _lmax.resize(DIMENSIONS, MIN_VAL);

  // calculate bounds
  for (const auto &coord : _coords) {
    for (int i = 0; i < DIMENSIONS; ++i) {
      _lmax[i] = std::max(_lmax[i], coord[i]);
      _lmin[i] = std::min(_lmin[i], coord[i]);
    }
  }
}



bool FastNS::build(double cutoff) {
  this->cutoff = cutoff;

  // compute box dimensions based on the current scale factor
  for (int i = 0; i < DIMENSIONS; ++i) {
    box[i] = _scale_factor * static_cast<float>(_lmax[i] - _lmin[i]);
  }

  // validate box dimensions against cutoff
  if (box[0] < cutoff || box[1] < cutoff || box[2] < cutoff) {

    if (_coords.size() < SMALL_SYSTEM_THRESHOLD) { // fall back to brute-force search
      return adaptive_build(cutoff);
    }
    spdlog::critical("Failed to build grid: likely because the cutoff is larger than the box dimensions.");
    return false;
  }

  // shift coordinates
  for (auto &coord : _coords) {
    for (int i = 0; i < DIMENSIONS; ++i) {
      coord[i] -= _lmin[i];
    }
  }

  coords_bbox = flatten_coordinates(_coords);
  lmin = _lmin;
  lmax = _lmax;
  build_grid();

  return true;
}


// An (unnecessarily?) complicated way to fall back to what is essentially a brute-force search.
bool FastNS::adaptive_build(double cutoff, int max_retries) {

  this->cutoff = cutoff;
  float updated_scale_factor = _scale_factor;

  // For box[i] to meet the cutoff, we require:
  //   _scale_factor * ( _lmax[i] - _lmin[i] ) >= cutoff
  //   _scale_factor >= cutoff / (_lmax[i] - _lmin[i])

  const float epsilon = 1e-2f; // a smaller epsilon would be (technically) better, but less robust.
  float min_required_scale = 0.0f;
  for (int i = 0; i < DIMENSIONS; ++i) {
    float delta = static_cast<float>(_lmax[i] - _lmin[i]); // we may need to check if delta > 0
    float required_scale = static_cast<float>(cutoff) / delta;
    min_required_scale = std::max(min_required_scale, required_scale);
  }
  _scale_factor = min_required_scale + epsilon;

  for (int i = 0; i < DIMENSIONS; ++i) {
    box[i] = _scale_factor * static_cast<float>(_lmax[i] - _lmin[i]);
  }

  // shift coordinates
  for (auto &coord : _coords) {
    for (int i = 0; i < DIMENSIONS; ++i) {
      coord[i] -= _lmin[i];
    }
  }

  coords_bbox = flatten_coordinates(_coords);
  lmin = _lmin, lmax = _lmax;
  build_grid();

  return true;
}


void FastNS::build_grid() {
  prepare_box();
  pack_grid();
}

NSResults FastNS::self_search() const {
  NSResults results;
  results.reserve_space(coords_bbox.size() / 3);

  float cutoff2 = cutoff * cutoff;
  for (int cx = 0; cx < ncells[0]; ++cx) {
    for (int cy = 0; cy < ncells[1]; ++cy) {
      for (int cz = 0; cz < ncells[2]; ++cz) {
        int ci = cell_xyz_to_cell_id(cx, cy, cz);

        int i = head_id[ci];
        while (i != END) {
          int j = next_id[i];
          const float *coord_i = &coords_bbox[3 * i];
          while (j != END) {
            float d2 = dist_sq(coord_i, &coords_bbox[3 * j]);
            if (d2 <= cutoff2) {
              results.add_neighbors(i, j, d2);
            }
            j = next_id[j];
          }

          for (const auto &offset : NeighborCells) {
            int ox = cx + offset[0];
            int oy = cy + offset[1];
            int oz = cz + offset[2];

            int cj = cell_xyz_to_cell_id(ox, oy, oz);
            if (cj == END) continue;

            j = head_id[cj];
            while (j != END) {
              float d2 = dist_sq(coord_i, &coords_bbox[3 * j]);
              if (d2 <= cutoff2) {
                results.add_neighbors(i, j, d2);
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

NSResults FastNS::search(const RDGeom::POINT3D_VECT &search_coords) const {
  NSResults results;
  results.reserve_space(search_coords.size());

  float cutoff_sq = cutoff * cutoff;

  if (search_coords.empty()) {
    return results;
  }

  RDGeom::POINT3D_VECT scoords(search_coords);
  for (auto it = scoords.begin(); it != scoords.end(); ++it) {
    for (int i = 0; i < 3; ++i) {
      (*it)[i] -= lmin[i];
    }
  }

  for (size_t i = 0; i < search_coords.size(); ++i) {
    std::array<float, 3> tmpcoord = {
        static_cast<float>(scoords[i].x),
        static_cast<float>(scoords[i].y),
        static_cast<float>(scoords[i].z)};

    std::array<int, 3> cellcoord;
    coord_to_cell_xyz(tmpcoord.data(), cellcoord);
    for (int xi = 0; xi < 3; ++xi) {
      for (int yi = 0; yi < 3; ++yi) {
        for (int zi = 0; zi < 3; ++zi) {
          int cx = cellcoord[0] - 1 + xi;
          int cy = cellcoord[1] - 1 + yi;
          int cz = cellcoord[2] - 1 + zi;

          int cellid = cell_xyz_to_cell_id(cx, cy, cz);
          if (cellid == END) {
            continue;
          }

          int j = head_id[cellid];
          while (j != END) {
            float d2 = dist_sq(tmpcoord.data(), &coords_bbox[3 * j]);
            if (d2 <= cutoff_sq) {
              results.add_neighbors(i, j, d2);
            }
            j = next_id[j];
          }
        }
      }
    }
  }

  return results;
}

NSResults FastNS::search(const std::vector<std::vector<double>> &search_coords) const {
  RDGeom::POINT3D_VECT scoords;
  scoords.reserve(search_coords.size());
  for (const auto &coord : search_coords) {
    scoords.emplace_back(coord[0], coord[1], coord[2]);
  }

  return search(scoords);
}


std::vector<float> FastNS::flatten_coordinates(const RDGeom::POINT3D_VECT &coords) {
    std::vector<float> flat_coords;
    flat_coords.reserve(coords.size() * 3);
    for (const auto &coord : coords) {
      flat_coords.insert(
          flat_coords.end(),
          {static_cast<float>(coord.x), static_cast<float>(coord.y), static_cast<float>(coord.z)});
    }
    return flat_coords;
}

bool FastNS::update(double cutoff) {
  if (cutoff == this->cutoff) {
    return true;
  }

  if (box[0] < cutoff || box[1] < cutoff || box[2] < cutoff) {
    return false;
  }

  ncells       = {0, 0, 0};
  cellsize     = {0.0f, 0.0f, 0.0f};
  cell_offsets = {0, 0, 0};

  head_id.clear();
  next_id.clear();

  this->cutoff = cutoff;

  build_grid();

  return true;
};

void FastNS::prepare_box() {
  double min_cellsize = cutoff + 0.001; // TODO: test the stability if we use no offset

  for (int i = 0; i < 3; ++i) {
    ncells[i] = std::min(static_cast<int>(std::floor(box[i] / min_cellsize)), MAX_GRID_DIM);
  }

  cellsize[0] = box[0] / ncells[0];
  cellsize[1] = box[1] / ncells[1];
  cellsize[2] = box[2] / ncells[2];

  cell_offsets[0] = 0;
  cell_offsets[1] = ncells[0];
  cell_offsets[2] = ncells[0] * ncells[1];
};

void FastNS::pack_grid() {
  head_id.assign(cell_offsets[2] * ncells[2], END);
  next_id.assign(coords_bbox.size() / 3, END);

  for (size_t i = 0; i < coords_bbox.size() / 3; ++i) {
    int j = coord_to_cell_id(&coords_bbox[3 * i]);
    next_id[i] = head_id[j];
    head_id[j] = i;
  }
};

inline int FastNS::coord_to_cell_id(const float *__restrict coord) const {
  std::array<int, DIMENSIONS> xyz;
  coord_to_cell_xyz(coord, xyz);
  return xyz[0] + xyz[1] * cell_offsets[1] + xyz[2] * cell_offsets[2];
}

inline void
FastNS::coord_to_cell_xyz(const float *__restrict coord, std::array<int, DIMENSIONS> &xyz) const {
  xyz[2] = static_cast<int>(coord[2] / cellsize[2]) % ncells[2];
  xyz[1] = static_cast<int>(coord[1] / cellsize[1]) % ncells[1];
  xyz[0] = static_cast<int>(coord[0] / cellsize[0]) % ncells[0];
}

int FastNS::cell_xyz_to_cell_id(int cx, int cy, int cz) const {
  if (cx < 0 || cx == ncells[0] || cy < 0 || cy == ncells[1] || cz < 0 || cz == ncells[2]) {
    return END;
  }
  return cx + cy * cell_offsets[1] + cz * cell_offsets[2];
}

void NSResults::add_neighbors(int i, int j, float d2) {
  m_pairs.emplace_back(i, j);
  m_dists.emplace_back(d2);
}

void NSResults::reserve_space(size_t input_size) {
  m_pairs.reserve(input_size);
  m_dists.reserve(input_size);
}

NSResults NSResults::filter(const double dist) const {
  NSResults filtered; // we'll not reserve space
  auto dist_sq = dist * dist;
  for (size_t i = 0; i < m_dists.size(); ++i) {
    if (m_dists[i] <= dist_sq) {
      filtered.add_neighbors(m_pairs[i].first, m_pairs[i].second, m_dists[i]);
    }
  }
  return filtered;
}

NSResults NSResults::filter(const std::vector<int> &atom_indices) const {
  NSResults filtered;
  filtered.reserve_space(m_pairs.size());
  std::unordered_set<int> atom_indices_set(atom_indices.begin(), atom_indices.end());
  for (size_t i = 0; i < m_pairs.size(); ++i) {
    if (atom_indices_set.find(m_pairs[i].first) != atom_indices_set.end()
        || atom_indices_set.find(m_pairs[i].second) != atom_indices_set.end()) {
      filtered.add_neighbors(m_pairs[i].first, m_pairs[i].second, m_dists[i]);
    }
  }
  return filtered;
}

NSResults NSResults::filter(const std::vector<int> &atom_indices, int col) const {
  if (col != 0 && col != 1) {
    throw std::invalid_argument("Column index must be 0 or 1");
  }

  NSResults filtered;
  filtered.reserve_space(m_pairs.size());
  std::unordered_set<int> atom_indices_set(atom_indices.begin(), atom_indices.end());

  for (size_t i = 0; i < m_pairs.size(); ++i) {
    if ((col == 0 && atom_indices_set.find(m_pairs[i].first) != atom_indices_set.end()) ||
        (col == 1 && atom_indices_set.find(m_pairs[i].second) != atom_indices_set.end())) {
      filtered.add_neighbors(m_pairs[i].first, m_pairs[i].second, m_dists[i]);
    }
  }

  return filtered;
}


} // namespace lahuta
