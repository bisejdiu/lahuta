#ifndef LAHUTA_KD_TREE_HPP
#define LAHUTA_KD_TREE_HPP

#include <cstddef>
#include <cstdint>
#include <vector>

// clang-format off
namespace lahuta {

template <typename Scalar>
class KDTree3 {
public:
  KDTree3() = default;

  void clear();
  bool ready() const { return node_count_ > 0 && coords_ptr_ != nullptr; }

  // Build over owned coordinates: layout [x,y,z] per point
  void build(const std::vector<Scalar> &coords, int leaf_size = 40);
  // Build using a non-owning view into external memory: n_points, 3
  void build_view(const Scalar *coords_ptr, std::size_t n_points, int leaf_size = 40);

  // Radius search: invokes callback(index, distance_sq) for each hit.
  template <typename Callback>
  void radius_search(const Scalar *point, Scalar radius_sq, Callback &&cb) const;

private:
  static constexpr int DIM = 3;

  struct Node {
    int start = 0;
    int end   = 0;
    std::uint8_t split_dim = 0;
    Scalar split_value = Scalar(0);
    bool leaf = true;
  };

  void build_from_pointer(std::size_t n_points, int leaf_size);
  static inline Scalar dist_sq(const Scalar *a, const Scalar *b);
  Scalar bbox_distance_sq(const Scalar *point, const Scalar *bbox_min, const Scalar *bbox_max) const;
  void build_node(int node_index, int start, int end);

  const std::vector<Scalar> *coords_vec_ = nullptr; // optional owner when we wrap an existing vector
  const Scalar *coords_ptr_ = nullptr;              // raw pointer to coords, either into coords_vec_ or external
  std::size_t n_points_ = 0;                        // number of points: coords_ptr_ length = n_points_*3
  std::vector<Node> nodes_;
  std::vector<Scalar> node_bounds_;                 // interleaved [n_nodes, 2, DIM]
  std::vector<int> indices_;
  int node_count_ = 0;
  mutable std::vector<int> stack_;                  // transient traversal storage to avoid reallocations per query
  int leaf_size_ = 40;                              // taken from sklearn's KDTree default
};

using KDTree3f = KDTree3<float>;
using KDTree3d = KDTree3<double>;

extern template class KDTree3<float>;
extern template class KDTree3<double>;

} // namespace lahuta

#include "kd_tree.tpp"

#endif // LAHUTA_KD_TREE_HPP
