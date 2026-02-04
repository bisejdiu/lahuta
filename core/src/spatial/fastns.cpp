/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   constexpr std::string_view first = "besian", last = "sejdiu", host = "gmail.com";
 *   return std::string(first) + std::string(last) + "@" + std::string(host);
 * }();
 *
 */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <mutex>
#include <stdexcept>

#include <Eigen/Dense>
#include <rdkit/Geometry/point.h>

#include "fastns.hpp"
#include "logging/logging.hpp"
#include "nsutils.hpp"
#include "spatial/nsresults_tls.hpp"

// clang-format off
namespace lahuta {

namespace {

constexpr double MIN_VAL = std::numeric_limits<double>::lowest();
constexpr double MAX_VAL = std::numeric_limits<double>::max();

constexpr int   END = -1;
constexpr int   MAX_GRID_DIM = 1290;
constexpr float CELL_EPSILON = 1e-3f;
constexpr std::int64_t MAX_TOTAL_CELLS = static_cast<std::int64_t>(64) * 1024 * 1024; // 64M cells

constexpr std::size_t BRUTE_FORCE_THRESHOLD = 500;    // upper cutoff for internal brute-force fallback
constexpr std::size_t DENSE_THRESHOLD       = 20'000; // lower cutoff for dense path in self_search()
constexpr std::size_t CELL_EIGEN_THRESHOLD  = 8;      // min points in cell to use Eigen block GEMM
constexpr std::size_t CROSS_EIGEN_THRESHOLD = 50;
constexpr std::size_t BRUTE_EIGEN_THRESHOLD = 100;

constexpr double COORD_EPSILON = 1e-10;
constexpr double SCALE_THRESHOLD = 1.0 / std::numeric_limits<float>::epsilon();

constexpr std::array<std::array<int, DIMENSIONS>, 13> NEIGHBOR_CELLS = {{
  {1, 0,  0}, {1, 1,  0}, {0, 1,  0}, {-1, 1,  0},
  {1, 0, -1}, {1, 1, -1}, {0, 1, -1}, {-1, 1, -1},
  {1, 0,  1}, {1, 1,  1}, {0, 1,  1}, {-1, 1,  1},
  {0, 0,  1}
}};

static std::once_flag eigen_threads_once; // Eigen global thread setting needs to be done once process-wide to avoid TSAN races

} // namespace

FastNS::FastNS(const RDGeom::POINT3D_VECT &coords) {

  // min/max bounds
  _lmin.resize(DIMENSIONS, MAX_VAL);
  _lmax.resize(DIMENSIONS, MIN_VAL);

  for (const auto &coord : coords) {
    for (int i = 0; i < DIMENSIONS; ++i) {
      const double val = coord[i];
      const double abs_val = std::abs(val);
      _lmax[i] = std::max(_lmax[i], val);
      _lmin[i] = std::min(_lmin[i], val);
      if (abs_val > COORD_EPSILON) {
        if (abs_val < min_abs_coord_) min_abs_coord_ = abs_val;
        if (abs_val > max_abs_coord_) max_abs_coord_ = abs_val;
      }
    }
  }

  if (min_abs_coord_ > 0.0 && std::isfinite(min_abs_coord_) && std::isfinite(max_abs_coord_)) {
    scale_ratio_ = max_abs_coord_ / min_abs_coord_;
    has_mixed_scales_ = scale_ratio_ > SCALE_THRESHOLD;
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

  initialize_eigen_views();
}

FastNS::FastNS(const std::vector<std::vector<double>> &coords) {

  // min/max bounds
  _lmin.resize(DIMENSIONS, MAX_VAL);
  _lmax.resize(DIMENSIONS, MIN_VAL);

  for (const auto &coord : coords) {
    for (int i = 0; i < DIMENSIONS; ++i) {
      const double val = coord[i];
      const double abs_val = std::abs(val);
      _lmax[i] = std::max(_lmax[i], val);
      _lmin[i] = std::min(_lmin[i], val);
      if (abs_val > COORD_EPSILON) {
        if (abs_val < min_abs_coord_) min_abs_coord_ = abs_val;
        if (abs_val > max_abs_coord_) max_abs_coord_ = abs_val;
      }
    }
  }

  if (min_abs_coord_ > 0.0 && std::isfinite(min_abs_coord_) && std::isfinite(max_abs_coord_)) {
    scale_ratio_ = max_abs_coord_ / min_abs_coord_;
    has_mixed_scales_ = scale_ratio_ > SCALE_THRESHOLD;
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

  initialize_eigen_views();
}

FastNS::FastNS(const double *coords_ptr, std::size_t npts) {
  if (!coords_ptr || npts == 0) return;

  _lmin.assign(DIMENSIONS, MAX_VAL);
  _lmax.assign(DIMENSIONS, MIN_VAL);

  const double *p = coords_ptr;
  for (std::size_t i = 0; i < npts; ++i, p += DIMENSIONS) {
    const double x = p[0];
    const double y = p[1];
    const double z = p[2];
    if (x < _lmin[0]) _lmin[0] = x; if (x > _lmax[0]) _lmax[0] = x;
    if (y < _lmin[1]) _lmin[1] = y; if (y > _lmax[1]) _lmax[1] = y;
    if (z < _lmin[2]) _lmin[2] = z; if (z > _lmax[2]) _lmax[2] = z;
    const double ax = std::abs(x);
    const double ay = std::abs(y);
    const double az = std::abs(z);
    if (ax > COORD_EPSILON) { if (ax < min_abs_coord_) min_abs_coord_ = ax; if (ax > max_abs_coord_) max_abs_coord_ = ax; }
    if (ay > COORD_EPSILON) { if (ay < min_abs_coord_) min_abs_coord_ = ay; if (ay > max_abs_coord_) max_abs_coord_ = ay; }
    if (az > COORD_EPSILON) { if (az < min_abs_coord_) min_abs_coord_ = az; if (az > max_abs_coord_) max_abs_coord_ = az; }
  }

  if (min_abs_coord_ > 0.0 && std::isfinite(min_abs_coord_) && std::isfinite(max_abs_coord_)) {
    scale_ratio_ = max_abs_coord_ / min_abs_coord_;
    has_mixed_scales_ = scale_ratio_ > SCALE_THRESHOLD;
  }

  coords_bbox.reserve(npts * 3);
  p = coords_ptr;
  for (std::size_t i = 0; i < npts; ++i, p += DIMENSIONS) {
    coords_bbox.push_back(static_cast<float>(p[0] - _lmin[0]));
    coords_bbox.push_back(static_cast<float>(p[1] - _lmin[1]));
    coords_bbox.push_back(static_cast<float>(p[2] - _lmin[2]));
  }
  n_points = npts;
  assert(coords_bbox.size() == 3 * n_points);

  initialize_eigen_views();
}

void FastNS::initialize_eigen_views() {
  if (n_points == 0 || coords_bbox.empty()) return;
  std::call_once(eigen_threads_once, [](){ Eigen::setNbThreads(1); }); // we handle our own threading
  P = std::make_unique<Eigen::Map<const RowMatX3f>>(coords_bbox.data(), static_cast<Eigen::Index>(n_points), 3); // Eigen Map view of coords_bbox
  P_norm2 = P->rowwise().squaredNorm(); // Precompute per-point squared norms once
}

void FastNS::prepare_scratch(int max_occ) {
  max_occ_ = std::max(1, max_occ);
  A_.resize(max_occ_, 3);
  B_.resize(max_occ_, 3);
  A_norm2_.resize(max_occ_);
  B_norm2_.resize(max_occ_);
  G_.resize(max_occ_, max_occ_);  // allocate a reasonable worst-case Gram block
  IdxA_.reserve(max_occ_);        // preallocate index buffers
  IdxB_.reserve(max_occ_);
}

void FastNS::reset_state() {
  grid_ready = false;
  brute_force_mode = false;
  brute_force_results.clear();
  head_id.clear();
  next_id.clear();
  ncells       = {0, 0, 0};
  cell_offsets = {0, 0, 0};
  cellsize     = {0.0f, 0.0f, 0.0f};
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
      Logger::get_logger()->trace(
          "FastNS.build: n={} cutoff={:.6f} using brute-force fallback ({} pairs)",
          n_points, cutoff, brute_force_results.size());
      return true;
    }

    Logger::get_logger()->warn(
        "FastNS.build: failed to configure grid for n={} cutoff={:.6f} (fallback disabled or too large)",
        n_points, cutoff);
    return false;
  }

  if (has_mixed_scales_) {
    Logger::get_logger()->warn(
        "FastNS.build: mixed-scale coordinates detected (ratio={:.2e}, orders~{:.1f}, range=[{:.3e}, {:.3e}])."
        "Small distances may underflow to 0.0 with float precision.", scale_ratio_, std::log10(scale_ratio_), min_abs_coord_, max_abs_coord_);
  }

  Logger::get_logger()->trace(
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

    // If still above the cap, we'll greedily reduce the largest dimension until acceptable
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

  if (head_id.empty() || next_id.empty()) {
    const_cast<FastNS*>(this)->build_grid();
  }

  //
  // The heuristic below uses absolute point count rather than mean occupancy, as the latter
  // is difficult to get right for heterogeneous molecular systems with our need for filtering
  // based on atom types. - Besian, October 2025
  //
  if (n_points >= DENSE_THRESHOLD) {
    Logger::get_logger()->trace("FastNS.self_search: DENSE path (n={})", n_points);
    return self_search_dense();
  } else {
    Logger::get_logger()->trace("FastNS.self_search: SPARSE path (n={})", n_points);
    return self_search_sparse();
  }
}

// Mostly empty cells, many cells with 0-2 points
NSResults FastNS::self_search_sparse() const {
  TlsResultsScope scope;
  NSResults &results = scope.results();
  results.reserve(ns_utils::estimate_pairs(n_points, static_cast<float>(cutoff), box));

  const double cutoff_sq = cutoff * cutoff;
  const float cutoff_sq_f = static_cast<float>(cutoff_sq);

  // No vector allocations, no batching, just direct linked-list traversal
  for (int cx = 0; cx < ncells[0]; ++cx) {
    for (int cy = 0; cy < ncells[1]; ++cy) {
      for (int cz = 0; cz < ncells[2]; ++cz) {
        int ci = cell_xyz_to_cell_id(cx, cy, cz);

        int i = head_id[ci];
        while (i != END) {
          const float *coord_i = &coords_bbox[3 * i];

          int j = next_id[i];
          while (j != END) {
            const float *coord_j = &coords_bbox[3 * j];
            const float d2 = sqdist3(coord_i, coord_j);
            if (d2 <= cutoff_sq_f) {
              results.add_neighbors(i, j, d2);
            }
            j = next_id[j];
          }

          for (const auto &offset : NEIGHBOR_CELLS) {
            int ox = cx + offset[0];
            int oy = cy + offset[1];
            int oz = cz + offset[2];
            int cj = cell_xyz_to_cell_id(ox, oy, oz);
            if (cj == END) continue;

            int j = head_id[cj];
            while (j != END) {
              const float *coord_j = &coords_bbox[3 * j];
              const float d2 = sqdist3(coord_i, coord_j);
              if (d2 <= cutoff_sq_f) {
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

// Bond computation with all atoms, cells with high occupancy
NSResults FastNS::self_search_dense() const {
  TlsResultsScope scope;
  NSResults &results = scope.results();
  results.reserve(ns_utils::estimate_pairs(n_points, static_cast<float>(cutoff), box));

  const double cutoff_sq = cutoff * cutoff;
  const float cutoff_sq_f = static_cast<float>(cutoff_sq);

#ifdef EIGEN_RUNTIME_NO_MALLOC
  Eigen::internal::set_is_malloc_allowed(false);
#endif

  for (int cx = 0; cx < ncells[0]; ++cx) {
    for (int cy = 0; cy < ncells[1]; ++cy) {
      for (int cz = 0; cz < ncells[2]; ++cz) {
        int ci = cell_xyz_to_cell_id(cx, cy, cz);

        IdxA_.clear();
        int idx = head_id[ci];
        while (idx != END) {
          IdxA_.push_back(idx);
          idx = next_id[idx];
        }

        const int cell_size = static_cast<int>(IdxA_.size());

        bool cell_in_A = false;
        if (cell_size >= CELL_EIGEN_THRESHOLD) {
          for (int r = 0; r < cell_size; ++r) {
            const int pt_idx = IdxA_[r];
            A_.row(r) = Eigen::Map<const Eigen::RowVector3f>(&coords_bbox[3 * pt_idx]);
            A_norm2_[r] = P_norm2[pt_idx];
          }
          cell_in_A = true;
        }

        // Within-cell pairs
        if (cell_in_A) {
          // Use preloaded A_ for within-cell Eigen GEMM
          auto A = A_.topRows(cell_size);
          auto A_norms = A_norm2_.head(cell_size);
          auto G = G_.topLeftCorner(cell_size, cell_size);

          G.noalias() = (-2.0f) * (A * A.transpose());
          G.colwise() += A_norms;
          G.rowwise() += A_norms.transpose();

          for (int r = 0; r < cell_size; ++r) {
            for (int c = r + 1; c < cell_size; ++c) {
              const float d2 = G(r, c);
              if (d2 <= cutoff_sq_f) {
                results.add_neighbors(IdxA_[r], IdxA_[c], d2);
              }
            }
          }
        } else {
          // Scalar path for sparse within-cell
          for (int k1 = 0; k1 < cell_size; ++k1) {
            const int i = IdxA_[k1];
            const float *coord_i = &coords_bbox[3 * i];
            for (int k2 = k1 + 1; k2 < cell_size; ++k2) {
              const int j = IdxA_[k2];
              const float *coord_j = &coords_bbox[3 * j];
              const float d2 = sqdist3(coord_i, coord_j);
              if (d2 <= cutoff_sq_f) {
                results.add_neighbors(i, j, d2);
              }
            }
          }
        }

        // Process neighbor cell pairs
        for (const auto &offset : NEIGHBOR_CELLS) {
          int ox = cx + offset[0];
          int oy = cy + offset[1];
          int oz = cz + offset[2];

          int cj = cell_xyz_to_cell_id(ox, oy, oz);
          if (cj == END) continue;

          IdxB_.clear();
          int j_idx = head_id[cj];
          while (j_idx != END) {
            IdxB_.push_back(j_idx);
            j_idx = next_id[j_idx];
          }

          const int neighbor_size = static_cast<int>(IdxB_.size());

          // Use Eigen block GEMM when both cells are dense enough
          if (cell_in_A && neighbor_size >= CELL_EIGEN_THRESHOLD) {
            for (int c = 0; c < neighbor_size; ++c) {
              const int j = IdxB_[c];
              B_.row(c) = Eigen::Map<const Eigen::RowVector3f>(&coords_bbox[3 * j]);
              B_norm2_[c] = P_norm2[j];
            }

            auto A = A_.topRows(cell_size);
            auto B = B_.topRows(neighbor_size);
            auto A_norms = A_norm2_.head(cell_size);
            auto B_norms = B_norm2_.head(neighbor_size);
            auto G = G_.topLeftCorner(cell_size, neighbor_size);

            G.noalias() = (-2.0f) * (A * B.transpose());
            G.colwise() += A_norms;
            G.rowwise() += B_norms.transpose();

            for (int r = 0; r < cell_size; ++r) {
              for (int c = 0; c < neighbor_size; ++c) {
                const float d2 = G(r, c);
                if (d2 <= cutoff_sq_f) {
                  results.add_neighbors(IdxA_[r], IdxB_[c], d2);
                }
              }
            }
          } else {
            for (int k = 0; k < cell_size; ++k) {
              const int i = IdxA_[k];
              const float *coord_i = &coords_bbox[3 * i];
              for (int n = 0; n < neighbor_size; ++n) {
                const int j = IdxB_[n];
                const float *coord_j = &coords_bbox[3 * j];
                const float d2 = sqdist3(coord_i, coord_j);
                if (d2 <= cutoff_sq_f) {
                  results.add_neighbors(i, j, d2);
                }
              }
            }
          }
        }
      }
    }
  }

#ifdef EIGEN_RUNTIME_NO_MALLOC
  Eigen::internal::set_is_malloc_allowed(true);
#endif

  return results;
}

NSResults FastNS::search(const RDGeom::POINT3D_VECT &search_coords) const {
  if (search_coords.empty()) {
    return {};
  }

  TlsResultsScope scope;
  NSResults &results = scope.results();

  const std::size_t m = search_coords.size();
  results.reserve(ns_utils::estimate_pairs(m, n_points, static_cast<float>(cutoff), box));

  const double cutoff_sq = cutoff * cutoff;

  if (brute_force_mode) {
    const std::size_t m = search_coords.size();
    const std::size_t n = n_points;

    if (m >= CROSS_EIGEN_THRESHOLD || n >= CROSS_EIGEN_THRESHOLD) {
      // We use a tiled approach to avoid materializing the full m x n distance matrix
      constexpr std::size_t TILE_SIZE = 256;

      const auto& coords = *P;  // n x 3 (all reference points)
      const float cutoff_sq_f = static_cast<float>(cutoff_sq);

      for (std::size_t i_start = 0; i_start < m; i_start += TILE_SIZE) {
        const std::size_t i_end = std::min(i_start + TILE_SIZE, m);
        const std::size_t tile_m = i_end - i_start;

        // Gather this tile of search coords
        Eigen::MatrixXf search_tile(tile_m, 3);
        Eigen::VectorXf search_norm2(tile_m);
        for (std::size_t k = 0; k < tile_m; ++k) {
          const std::size_t i = i_start + k;
          search_tile(k, 0) = static_cast<float>(search_coords[i].x - _lmin[0]);
          search_tile(k, 1) = static_cast<float>(search_coords[i].y - _lmin[1]);
          search_tile(k, 2) = static_cast<float>(search_coords[i].z - _lmin[2]);
        }
        search_norm2 = search_tile.rowwise().squaredNorm();

        // Compute tile_m x n distances using GEMM: D = -2*(Q*P^T) + ||Q||^2 + ||P||^2
        Eigen::MatrixXf D(tile_m, n);
        D.noalias() = -2.0f * (search_tile * coords.transpose());
        D.colwise() += search_norm2;
        D.rowwise() += P_norm2.transpose();

        // Stream results from this tile
        for (std::size_t k = 0; k < tile_m; ++k) {
          const std::size_t i = i_start + k;
          for (std::size_t j = 0; j < n; ++j) {
            const float d2 = D(k, j);
            if (d2 <= cutoff_sq_f) {
              results.add_neighbors(static_cast<int>(i), static_cast<int>(j), d2);
            }
          }
        }
      }

    } else {
      const float cutoff_sq_f = static_cast<float>(cutoff_sq);
      for (std::size_t i = 0; i < m; ++i) {
        const float sx = static_cast<float>(search_coords[i].x - _lmin[0]);
        const float sy = static_cast<float>(search_coords[i].y - _lmin[1]);
        const float sz = static_cast<float>(search_coords[i].z - _lmin[2]);
        const float search_coord[3] = {sx, sy, sz};

        for (std::size_t j = 0; j < n; ++j) {
          const float *coord_j = &coords_bbox[3 * j];
          const float d2 = sqdist3(search_coord, coord_j);
          if (d2 <= cutoff_sq_f) {
            results.add_neighbors(static_cast<int>(i), static_cast<int>(j), d2);
          }
        }
      }
    }
    return results;
  }

  if (!grid_ready) throw std::runtime_error("FastNS grid not built. Call build() before search().");

  if (head_id.empty() || next_id.empty()) {
    const_cast<FastNS*>(this)->build_grid();
  }

  const float cutoff_sq_f = static_cast<float>(cutoff_sq);

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
            const float d2 = sqdist3(tmpcoord.data(), coord_j);
            if (d2 <= cutoff_sq_f) {
              results.add_neighbors(static_cast<int>(i), j, d2);
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

  // Calculate maximum cell occupancy for scratch space allocation
  int max_occ = 0;
  for (int cell_id : head_id) {
    if (cell_id != END) {
      int count = 0;
      int current = cell_id;
      while (current != END) {
        count++;
        current = next_id[current];
      }
      max_occ = std::max(max_occ, count);
    }
  }

  prepare_scratch(max_occ); // Prepare scratch space based on maximum occupancy
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
  TlsResultsScope scope;
  NSResults &results = scope.results();
  results.reserve(ns_utils::estimate_pairs(n_points, static_cast<float>(std::sqrt(cutoff_sq)), box));

  const int n = static_cast<int>(n_points);

  if (n_points >= BRUTE_EIGEN_THRESHOLD) {
    constexpr std::size_t TILE_SIZE = 256;

    const auto& coords = *P;
    const float cutoff_sq_f = static_cast<float>(cutoff_sq);

    for (std::size_t i_start = 0; i_start < n_points; i_start += TILE_SIZE) {
      const std::size_t i_end = std::min(i_start + TILE_SIZE, n_points);
      const std::size_t tile_rows = i_end - i_start;

      auto A = coords.middleRows(i_start, tile_rows);
      auto A_norm2 = P_norm2.segment(i_start, tile_rows);

      // Within-tile upper triangle (A vs A)
      Eigen::MatrixXf D_aa(tile_rows, tile_rows);
      D_aa.noalias() = -2.0f * (A * A.transpose());
      D_aa.colwise() += A_norm2;
      D_aa.rowwise() += A_norm2.transpose();

      for (std::size_t r = 0; r < tile_rows; ++r) {
        for (std::size_t c = r + 1; c < tile_rows; ++c) {
          const float d2 = D_aa(r, c);
          if (d2 <= cutoff_sq_f) {
            results.add_neighbors(
              static_cast<int>(i_start + r), 
              static_cast<int>(i_start + c), 
              d2
            );
          }
        }
      }

      // Cross-tile (A vs remaining B)
      if (i_end < n_points) {
        auto B = coords.bottomRows(n_points - i_end);
        auto B_norm2 = P_norm2.tail(n_points - i_end);

        Eigen::MatrixXf D_ab(tile_rows, n_points - i_end);
        D_ab.noalias() = -2.0f * (A * B.transpose());
        D_ab.colwise() += A_norm2;
        D_ab.rowwise() += B_norm2.transpose();

        for (std::size_t r = 0; r < tile_rows; ++r) {
          for (std::size_t c = 0; c < (n_points - i_end); ++c) {
            const float d2 = D_ab(r, c);
            if (d2 <= cutoff_sq_f) {
              results.add_neighbors(
                static_cast<int>(i_start + r),
                static_cast<int>(i_end + c),
                d2
              );
            }
          }
        }
      }
    }

  } else {
    const float cutoff_sq_f = static_cast<float>(cutoff_sq);
    for (int i = 0; i < n; ++i) {
      const float *ai = &coords_bbox[3 * i];
      for (int j = i + 1; j < n; ++j) {
        const float *aj = &coords_bbox[3 * j];
        const float d2 = sqdist3(ai, aj);
        if (d2 <= cutoff_sq_f) {
          results.add_neighbors(i, j, d2);
        }
      }
    }
  }

  return results;
}

void FastNS::log_occupancy_stats() const {
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

} // namespace lahuta
