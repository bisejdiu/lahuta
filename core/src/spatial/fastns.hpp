#ifndef LAHUTA_FASTNS_HPP
#define LAHUTA_FASTNS_HPP

#include <Eigen/Dense>
#include <array>
#include <cassert>
#include <cmath>
#include <memory>
#include <vector>

#include <rdkit/Geometry/point.h>

#include "nsresults.hpp"

// clang-format off
namespace lahuta {

constexpr int DIMENSIONS = 3;

class FastNS {
public:
  FastNS(const RDGeom::POINT3D_VECT &coords);
  FastNS(const std::vector<std::vector<double>> &coords);
  FastNS(const double *coords_ptr, std::size_t n_points); // Used by Python bindings
  ~FastNS() = default;

  bool build(double cutoff, bool brute_force_fallback = true);

  // find all neighbors among the input coordinates within the cutoff distance
  NSResults self_search() const;

  // find all neighbors between provided search coords and the input coordinates within the cutoff distance
  NSResults search(const RDGeom::POINT3D_VECT &search_coords) const;
  NSResults search(const std::vector<std::vector<double>> &search_coords) const;

  double get_cutoff() const { return cutoff; }

  std::array<int,   DIMENSIONS> get_ncells()   const { return ncells; }
  std::array<float, DIMENSIONS> get_cellsize() const { return cellsize; }
  std::array<float, DIMENSIONS> get_box()      const { return box; }
  bool is_grid_ready() const { return grid_ready; }

  double get_scale_ratio()   const { return scale_ratio_; }
  bool   has_mixed_scales()  const { return has_mixed_scales_; }
  double get_min_abs_coord() const { return min_abs_coord_; }
  double get_max_abs_coord() const { return max_abs_coord_; }

  template <typename T>
  static inline T dist_sq(const T* __restrict a, const T* __restrict b) {
      T dx = a[0] - b[0];
      T dy = a[1] - b[1];
      T dz = a[2] - b[2];
      return dx * dx + dy * dy + dz * dz;
  }

  inline static float sqdist3(const float* __restrict a, const float* __restrict b) {
    const float dx = a[0] - b[0];
    const float dy = a[1] - b[1];
    const float dz = a[2] - b[2];
    return std::fmaf(dz, dz, std::fmaf(dy, dy, dx*dx)); // 2 FMAs
  }

private:
  bool configure_grid(const std::array<float, DIMENSIONS> &extents, float padding, float min_cell);
  void reset_state();
  void pack_grid();
  void build_grid();
  void log_occupancy_stats() const;

  // Initialize Eigen Map and precomputed norms
  void initialize_eigen_views();

  // Prepare reusable scratch space for cell-pair operations
  void prepare_scratch(int max_occ);

  NSResults brute_force_self(double cutoff_sq) const;

  // Adaptive search strategies based on cell occupancy
  NSResults self_search_sparse() const;
  NSResults self_search_dense()  const;

  inline int  coord_to_cell_id (const float *__restrict coord) const;
  inline void coord_to_cell_xyz(const float *__restrict coord, std::array<int, 3> &xyz) const;
  inline int  cell_xyz_to_cell_id(int cx, int cy, int cz) const;

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

  double min_abs_coord_ = std::numeric_limits<double>::max();
  double max_abs_coord_ = 0.0;
  double scale_ratio_   = 1.0;
  bool   has_mixed_scales_ = false;

  using MatX3f = Eigen::Matrix<float, Eigen::Dynamic, 3>;  // Column-major by default for optimal GEMM
  using VecXf  = Eigen::VectorXf;

  // Mapped view of coords_bbox
  // Row-major externally, but we'll gather into column-major scratch
  using RowMatX3f = Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>;

  static_assert(std::is_same_v<RowMatX3f::Scalar, float>, "RowMatX3f must use float");
  static_assert(std::is_same_v<Eigen::MatrixXf::Scalar, float>, "MatrixXf must use float");

  std::unique_ptr<Eigen::Map<const RowMatX3f>> P;

  VecXf P_norm2; // Precomputed per-point squared norms

  //
  // Reusable scratch space for cell-pair operations. Marked mutable to allow for
  // use in const search methods.
  //
  // THREAD SAFETY: These buffers are NOT thread-safe for concurrent calls to
  // self_search() or search() on the same FastNS instance. If we add outer-loop
  // parallelism (e.g., OpenMP), use per-thread buffers via TLS
  // (similar to TlsResultsScope) or we'll make them local with reserve(max_occ_).
  //
  mutable int max_occ_ = 0;
  mutable MatX3f A_, B_;
  mutable Eigen::VectorXf A_norm2_, B_norm2_;
  mutable Eigen::MatrixXf G_;          // reuse for Gram blocks
  mutable std::vector<int> IdxA_;      // reusable index buffer for cell_points
  mutable std::vector<int> IdxB_;      // reusable index buffer for neighbor_points

  double cutoff = 0.0;
  bool grid_ready = false;
  bool brute_force_mode = false;
  NSResults brute_force_results;

};

} // namespace lahuta

#endif // LAHUTA_FASTNS_HPP
