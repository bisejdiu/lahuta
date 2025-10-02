#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <unordered_set>

#include <rdkit/Geometry/point.h>

#include "definitions.hpp"
#include "lahuta.hpp"
#include "logging.hpp"
#include "nsgrid.hpp"

// clang-format off
namespace lahuta {

namespace {

constexpr int   END = -1;
constexpr int   MAX_GRID_DIM = 1290;
constexpr std::size_t BRUTE_FORCE_THRESHOLD = 5'000; // upper cutoff for internal brute-force fallback
constexpr float CELL_EPSILON = 1e-3f;
constexpr std::int64_t MAX_TOTAL_CELLS = static_cast<std::int64_t>(64) * 1024 * 1024; // 64M cells

constexpr std::array<std::array<int, DIMENSIONS>, 13> NeighborCells = {{
  {1, 0,  0}, {1, 1,  0}, {0, 1,  0}, {-1, 1,  0},
  {1, 0, -1}, {1, 1, -1}, {0, 1, -1}, {-1, 1, -1},
  {1, 0,  1}, {1, 1,  1}, {0, 1,  1}, {-1, 1,  1},
  {0, 0,  1}
}};

}

FastNS::FastNS(const RDGeom::POINT3D_VECT &coords) {

  // min/max bounds
  _lmin.resize(DIMENSIONS, MAX_VAL);
  _lmax.resize(DIMENSIONS, MIN_VAL);

  for (const auto &coord : coords) {
    for (int i = 0; i < DIMENSIONS; ++i) {
      _lmax[i] = std::max(_lmax[i], coord[i]);
      _lmin[i] = std::min(_lmin[i], coord[i]);
    }
  }

  // pack into coords_bbox
  coords_bbox.reserve(coords.size() * 3);
  for (const auto &c : coords) {
    coords_bbox.push_back(static_cast<float>(c.x - _lmin[0]));
    coords_bbox.push_back(static_cast<float>(c.y - _lmin[1]));
    coords_bbox.push_back(static_cast<float>(c.z - _lmin[2]));
  }
  n_points = coords.size();
  assert(coords_bbox.size() == 3 * n_points);
}

FastNS::FastNS(const std::vector<std::vector<double>> &coords) {

  // min/max bounds
  _lmin.resize(DIMENSIONS, MAX_VAL);
  _lmax.resize(DIMENSIONS, MIN_VAL);

  for (const auto &coord : coords) {
    for (int i = 0; i < DIMENSIONS; ++i) {
      _lmax[i] = std::max(_lmax[i], coord[i]);
      _lmin[i] = std::min(_lmin[i], coord[i]);
    }
  }

  // pack into coords_bbox
  coords_bbox.reserve(coords.size() * 3);
  for (const auto &c : coords) {
    coords_bbox.push_back(static_cast<float>(c[0] - _lmin[0]));
    coords_bbox.push_back(static_cast<float>(c[1] - _lmin[1]));
    coords_bbox.push_back(static_cast<float>(c[2] - _lmin[2]));
  }
  n_points = coords.size();
  assert(coords_bbox.size() == 3 * n_points);
}

FastNS::FastNS(const double *coords_ptr, std::size_t npts) {
  if (!coords_ptr || npts == 0) return;

  _lmin.assign(DIMENSIONS, MAX_VAL);
  _lmax.assign(DIMENSIONS, MIN_VAL);

  // pass 1: bounds
  const double *p = coords_ptr;
  for (std::size_t i = 0; i < npts; ++i, p += DIMENSIONS) {
    const double x = p[0];
    const double y = p[1];
    const double z = p[2];
    if (x < _lmin[0]) _lmin[0] = x; if (x > _lmax[0]) _lmax[0] = x;
    if (y < _lmin[1]) _lmin[1] = y; if (y > _lmax[1]) _lmax[1] = y;
    if (z < _lmin[2]) _lmin[2] = z; if (z > _lmax[2]) _lmax[2] = z;
  }

  // pass 2: pack into coords_bbox
  coords_bbox.reserve(npts * 3);
  p = coords_ptr;
  for (std::size_t i = 0; i < npts; ++i, p += DIMENSIONS) {
    coords_bbox.push_back(static_cast<float>(p[0] - _lmin[0]));
    coords_bbox.push_back(static_cast<float>(p[1] - _lmin[1]));
    coords_bbox.push_back(static_cast<float>(p[2] - _lmin[2]));
  }
  n_points = npts;
  assert(coords_bbox.size() == 3 * n_points);
}

void FastNS::reset_state() {
  grid_ready = false;
  brute_force_mode = false;
  brute_force_results.clear();
  head_id.clear();
  next_id.clear();
  ncells = {0, 0, 0};
  cell_offsets = {0, 0, 0};
  cellsize = {0.0f, 0.0f, 0.0f};
}

bool FastNS::build(double cutoff, bool brute_force_fallback) {
  this->cutoff = cutoff;
  reset_state();

  // invalid cutoff must fail regardless of fallback
  if (!std::isfinite(cutoff)) {
    Logger::get_logger()->warn("FastNS.build: invalid cutoff ({}). Must be finite.", cutoff);
    return false;
  }

  // coordinates must already be packed by a constructor
  if (coords_bbox.empty()) return false;

  std::array<float, DIMENSIONS> extents{};
  for (int i = 0; i < DIMENSIONS; ++i) {
    extents[i] = static_cast<float>(_lmax[i] - _lmin[i]);
  }

  const float min_cell = static_cast<float>(cutoff) + CELL_EPSILON;
  constexpr std::array<float, 4> padding_attempts = {1.05f, 1.25f, 1.5f, 2.0f};

  bool configured = false;
  for (float padding : padding_attempts) {
    if (configure_grid(extents, padding, min_cell)) {
      configured = true;
      break;
    }
  }

  if (!configured) {
    if (brute_force_fallback && n_points <= BRUTE_FORCE_THRESHOLD) {
      brute_force_mode = true;
      const double cutoff_sq = cutoff * cutoff;
      brute_force_results = brute_force_self(cutoff_sq);
      Logger::get_logger()->debug(
          "FastNS.build: n={} cutoff={:.6f} using brute-force fallback ({} pairs)",
          n_points, cutoff, brute_force_results.size());
      return true;
    }

    Logger::get_logger()->warn(
        "FastNS.build: failed to configure grid for n={} cutoff={:.6f} (fallback disabled or too large)",
        n_points, cutoff);
    return false;
  }

  Logger::get_logger()->debug(
      "FastNS.build: n={} cutoff={:.6f} box=({:.6f}, {:.6f}, {:.6f}) ncells=({}, {}, {})",
      n_points, cutoff, box[0], box[1], box[2], ncells[0], ncells[1], ncells[2]);

  grid_ready = true;
  return true; // self_search() will lazily build the grid when needed
}

bool FastNS::configure_grid(const std::array<float, DIMENSIONS> &extents, float padding, float min_cell) {
  if (!(padding > 0.0f) || !(min_cell > 0.0f)) return false;

  std::array<float, DIMENSIONS> candidate_box{};
  std::array<int,   DIMENSIONS> candidate_ncells{};

  const float min_cellsize = std::max(min_cell, CELL_EPSILON);

  double total_cells = 1.0;
  for (int i = 0; i < DIMENSIONS; ++i) {
    float extent = extents[i];
    if (!std::isfinite(extent) || extent < 0.0f) extent = 0.0f;

    float candidate = extent * padding;
    if (!std::isfinite(candidate) || candidate < min_cellsize) candidate = min_cellsize;

    candidate_box[i] = candidate;

    double cells_d = std::floor(static_cast<double>(candidate) / static_cast<double>(min_cellsize));
    if (!std::isfinite(cells_d) || cells_d < 1.0) cells_d = 1.0;
    if (cells_d > static_cast<double>(MAX_GRID_DIM)) cells_d = static_cast<double>(MAX_GRID_DIM);

    candidate_ncells[i] = static_cast<int>(cells_d);
    total_cells *= cells_d;
    if (!std::isfinite(total_cells) || total_cells <= 0.0) return false;
  }

  if (total_cells > static_cast<double>(MAX_TOTAL_CELLS)) {
    const double scale = std::cbrt(total_cells / static_cast<double>(MAX_TOTAL_CELLS));
    for (int i = 0; i < DIMENSIONS; ++i) {
      double reduced = static_cast<double>(candidate_ncells[i]) / scale;
      if (reduced < 1.0) reduced = 1.0;
      candidate_ncells[i] = static_cast<int>(std::floor(reduced));
      if (candidate_ncells[i] < 1) candidate_ncells[i] = 1;
    }

    total_cells = 1.0;
    for (int i = 0; i < DIMENSIONS; ++i) total_cells *= static_cast<double>(candidate_ncells[i]);

    // If still above the cap, greedily reduce the largest dimension until acceptable
    int guard = 0;
    while (total_cells > static_cast<double>(MAX_TOTAL_CELLS) && guard < 1024) {
      int max_dim = 0;
      for (int i = 1; i < DIMENSIONS; ++i) {
        if (candidate_ncells[i] > candidate_ncells[max_dim]) max_dim = i;
      }
      if (candidate_ncells[max_dim] == 1) break;
      --candidate_ncells[max_dim];
      total_cells = 1.0;
      for (int i = 0; i < DIMENSIONS; ++i) {
        total_cells *= static_cast<double>(candidate_ncells[i]);
      }
      ++guard;
    }

    if (total_cells > static_cast<double>(MAX_TOTAL_CELLS)) return false;
  }

  for (int i = 0; i < DIMENSIONS; ++i) {
    if (candidate_ncells[i] <= 0) return false;

    float size = candidate_box[i] / static_cast<float>(candidate_ncells[i]);
    if (!std::isfinite(size) || size <= 0.0f) return false;
    cellsize[i] = size;
  }

  box = candidate_box;
  ncells = candidate_ncells;
  cell_offsets[0] = 0;
  cell_offsets[1] = ncells[0];
  cell_offsets[2] = ncells[0] * ncells[1];

  return true;
}

void FastNS::build_grid() {
  if (!grid_ready) {
    throw std::runtime_error("FastNS grid configuration missing. Call build() before self_search().");
  }
  if (brute_force_mode) return;

  pack_grid();
  log_occupancy_stats();
}

NSResults FastNS::self_search() const {
  if (brute_force_mode) {
    return brute_force_results;
  }
  if (!grid_ready) {
    throw std::runtime_error("FastNS grid not built. Call build() before self_search().");
  }

  NSResults results;
  if (head_id.empty() || next_id.empty()) {
    const_cast<FastNS*>(this)->build_grid();
  }
  results.reserve_space(n_points);

  const double cutoff_sq = cutoff * cutoff;
  for (int cx = 0; cx < ncells[0]; ++cx) {
    for (int cy = 0; cy < ncells[1]; ++cy) {
      for (int cz = 0; cz < ncells[2]; ++cz) {
        int ci = cell_xyz_to_cell_id(cx, cy, cz);

        int i = head_id[ci];
        while (i != END) {
          int j = next_id[i];
          const float *coord_i = &coords_bbox[3 * i];
          while (j != END) {
            const float *coord_j = &coords_bbox[3 * j];
            const double dx = static_cast<double>(coord_i[0]) - coord_j[0];
            const double dy = static_cast<double>(coord_i[1]) - coord_j[1];
            const double dz = static_cast<double>(coord_i[2]) - coord_j[2];
            const double d2 = dx * dx + dy * dy + dz * dz;
            if (d2 <= cutoff_sq) {
              results.add_neighbors(i, j, static_cast<float>(d2));
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
              const float *coord_j = &coords_bbox[3 * j];
              const double dx = static_cast<double>(coord_i[0]) - coord_j[0];
              const double dy = static_cast<double>(coord_i[1]) - coord_j[1];
              const double dz = static_cast<double>(coord_i[2]) - coord_j[2];
              const double d2 = dx * dx + dy * dy + dz * dz;
              if (d2 <= cutoff_sq) {
                results.add_neighbors(i, j, static_cast<float>(d2));
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

  if (search_coords.empty()) {
    return results;
  }

  const double cutoff_sq = cutoff * cutoff;

  if (brute_force_mode) {
    for (std::size_t i = 0; i < search_coords.size(); ++i) {
      const double sx = static_cast<double>(search_coords[i].x) - _lmin[0];
      const double sy = static_cast<double>(search_coords[i].y) - _lmin[1];
      const double sz = static_cast<double>(search_coords[i].z) - _lmin[2];

      for (std::size_t j = 0; j < n_points; ++j) {
        const float *coord_j = &coords_bbox[3 * j];
        const double dx = sx - static_cast<double>(coord_j[0]);
        const double dy = sy - static_cast<double>(coord_j[1]);
        const double dz = sz - static_cast<double>(coord_j[2]);
        const double d2 = dx * dx + dy * dy + dz * dz;
        if (d2 <= cutoff_sq) {
          results.add_neighbors(static_cast<int>(i), static_cast<int>(j), static_cast<float>(d2));
        }
      }
    }
    return results;
  }

  if (!grid_ready) {
    throw std::runtime_error("FastNS grid not built. Call build() before search().");
  }

  if (head_id.empty() || next_id.empty()) {
    const_cast<FastNS*>(this)->build_grid();
  }

  for (std::size_t i = 0; i < search_coords.size(); ++i) {
    std::array<float, 3> tmpcoord = {
        static_cast<float>(search_coords[i].x - _lmin[0]),
        static_cast<float>(search_coords[i].y - _lmin[1]),
        static_cast<float>(search_coords[i].z - _lmin[2])};

    std::array<int, 3> cellcoord;
    coord_to_cell_xyz(tmpcoord.data(), cellcoord);

    for (int xi = -1; xi <= 1; ++xi) {
      const int cx = cellcoord[0] + xi;
      for (int yi = -1; yi <= 1; ++yi) {
        const int cy = cellcoord[1] + yi;
        for (int zi = -1; zi <= 1; ++zi) {
          const int cz = cellcoord[2] + zi;

          const int cellid = cell_xyz_to_cell_id(cx, cy, cz);
          if (cellid == END) {
            continue;
          }

          int j = head_id[cellid];
          while (j != END) {
            const float *coord_j = &coords_bbox[3 * j];
            const double dx = static_cast<double>(tmpcoord[0]) - static_cast<double>(coord_j[0]);
            const double dy = static_cast<double>(tmpcoord[1]) - static_cast<double>(coord_j[1]);
            const double dz = static_cast<double>(tmpcoord[2]) - static_cast<double>(coord_j[2]);
            const double d2 = dx * dx + dy * dy + dz * dz;
            if (d2 <= cutoff_sq) {
              results.add_neighbors(static_cast<int>(i), j, static_cast<float>(d2));
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
  for (const auto &c : search_coords) {
    scoords.emplace_back(c[0], c[1], c[2]);
  }
  return search(scoords);
}

void FastNS::pack_grid() {
  head_id.assign(cell_offsets[2] * ncells[2], END);
  next_id.assign(coords_bbox.size() / 3, END);

  for (size_t i = 0; i < coords_bbox.size() / 3; ++i) {
    int j = coord_to_cell_id(&coords_bbox[3 * i]);
    next_id[i] = head_id[j];
    head_id[j] = i;
  }
};

int FastNS::coord_to_cell_id(const float *__restrict coord) const {
  std::array<int, DIMENSIONS> xyz;
  coord_to_cell_xyz(coord, xyz);
  return xyz[0] + xyz[1] * cell_offsets[1] + xyz[2] * cell_offsets[2];
}

void FastNS::coord_to_cell_xyz(const float *__restrict coord, std::array<int, DIMENSIONS> &xyz) const {
  if (ncells[0] <= 0 || ncells[1] <= 0 || ncells[2] <= 0) {
    throw std::runtime_error("FastNS grid dimensions invalid (ncells <= 0).");
  }

  xyz[0] = std::clamp(static_cast<int>(coord[0] / cellsize[0]), 0, ncells[0] - 1);
  xyz[1] = std::clamp(static_cast<int>(coord[1] / cellsize[1]), 0, ncells[1] - 1);
  xyz[2] = std::clamp(static_cast<int>(coord[2] / cellsize[2]), 0, ncells[2] - 1);
}

int FastNS::cell_xyz_to_cell_id(int cx, int cy, int cz) const {
  if (cx < 0 || cx == ncells[0] || cy < 0 || cy == ncells[1] || cz < 0 || cz == ncells[2]) {
    return END;
  }
  return cx + cy * cell_offsets[1] + cz * cell_offsets[2];
}

NSResults FastNS::brute_force_self(double cutoff_sq) const {
  NSResults results;
  results.reserve_space(n_points);

  const int n = static_cast<int>(n_points);
  for (int i = 0; i < n; ++i) {
    const float *ai = &coords_bbox[3 * i];
    for (int j = i + 1; j < n; ++j) {
      const float *aj = &coords_bbox[3 * j];
      const double dx = static_cast<double>(ai[0]) - static_cast<double>(aj[0]);
      const double dy = static_cast<double>(ai[1]) - static_cast<double>(aj[1]);
      const double dz = static_cast<double>(ai[2]) - static_cast<double>(aj[2]);
      const double d2 = dx * dx + dy * dy + dz * dz;
      if (d2 <= cutoff_sq) {
        results.add_neighbors(i, j, static_cast<float>(d2));
      }
    }
  }

  return results;
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
  NSResults filtered;
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

void FastNS::log_occupancy_stats() const {
  // Occupancy stats
  const int ncell_total = ncells[0] * ncells[1] * ncells[2];
  int max_occ = 0;
  long long sum_occ = 0;
  if (ncell_total > 0) {
    for (int cid = 0; cid < ncell_total; ++cid) {
      int cnt = 0;
      int j = head_id[cid];
      while (j != END) { ++cnt; j = next_id[j]; }
      sum_occ += cnt; if (cnt > max_occ) max_occ = cnt;
    }
  }
  const double mean_occ = ncell_total > 0 ? static_cast<double>(sum_occ) / static_cast<double>(ncell_total) : 0.0;

  Logger::get_logger()->debug(
      "FastNS.grid: ncells=({}, {}, {}), cellsize=({:.6f}, {:.6f}, {:.6f}), mean_occ={:.3f}, max_occ={}",
      ncells[0], ncells[1], ncells[2], cellsize[0], cellsize[1], cellsize[2], mean_occ, max_occ);
}

namespace ns_utils {

NSResults remove_adjascent_residueid_pairs(const Luni &luni, NSResults &results, int res_diff) {
  Pairs     f_pairs;
  Distances distances;

  for (const auto &[pair, dist] : results) {
    auto *fatom = luni.get_molecule().getAtomWithIdx(pair.first);
    auto *satom = luni.get_molecule().getAtomWithIdx(pair.second);

    auto *finfo = static_cast<const RDKit::AtomPDBResidueInfo *>(fatom->getMonomerInfo());
    auto *sinfo = static_cast<const RDKit::AtomPDBResidueInfo *>(satom->getMonomerInfo());

    // FIX: this is another filter (filtering out H-H pairs)
    if (fatom->getAtomicNum() == 1 || satom->getAtomicNum() == 1) continue;

    auto f_resid = finfo->getResidueNumber();
    auto s_resid = sinfo->getResidueNumber();

    auto is_either_nonprotein = !definitions::is_polymer(finfo->getResidueName()) || !definitions::is_polymer(sinfo->getResidueName());
    if (std::abs(f_resid - s_resid) > res_diff || is_either_nonprotein) {
      f_pairs.push_back(pair);
      distances.push_back(dist);
    }
  }
  return NSResults(f_pairs, distances);
}

} // namespace ns_utils

} // namespace lahuta
