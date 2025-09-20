#ifndef LAHUTA_KD_INDEX_HPP
#define LAHUTA_KD_INDEX_HPP

#include <memory>
#include <vector>

#include <rdkit/Geometry/point.h>

#include "kd_tree.hpp"
#include "nsgrid.hpp" // for NSResults

// clang-format off
namespace lahuta {

// Thin wrapper over KDTree3f for cross radius searches.
class KDTreeIndex {
public:
  KDTreeIndex() = default;

  // Build on target points. Converts to contiguous float32 [x,y,z]*n.
  bool build(const RDGeom::POINT3D_VECT &points, int leaf_size = 40);
  // Directly from a contiguous array of f64 coordinates with shape (n,3).
  bool build(const double *coords_f64, std::size_t n_points, int leaf_size = 40);

  bool build_view_f32(const float  *coords_f32, std::size_t n_points, int leaf_size = 40); // Using a non-owning view into external f32 memory with shape (n,3).
  bool build_view_f64(const double *coords_f64, std::size_t n_points, int leaf_size = 40); // Using a non-owning view into external f64 memory with shape (n,3).

  // Optional: set an external owner to extend lifetime of externally referenced memory.
  void set_external_owner(std::shared_ptr<const void> owner) { external_owner_ = std::move(owner); } // for python interop

  bool ready() const { return kd_.ready() || kd64_.ready(); }

  template <typename Callback>
  void radius_search_point(const float p[3], float r2, Callback &&cb) const {
    kd_.radius_search(p, r2, std::forward<Callback>(cb));
  }

  // Batched search. Returns NSResults of (query_index, target_index) with squared distances
  NSResults radius_search(const RDGeom::POINT3D_VECT &queries, double radius) const;

private:
  void reset_all();

  std::vector<float>  coords_  ; // owned target coords (float32)
  std::vector<double> coords64_; // owned target coords (float64)
  KDTree3f kd_;
  KDTree3d kd64_;
  std::shared_ptr<const void> external_owner_{}; // retains external memory (e.g., numpy array)
};

} // namespace lahuta

#endif // LAHUTA_KD_INDEX_HPP
