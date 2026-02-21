/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto to_s = [](auto&& arg) {
 *     if constexpr (std::is_same_v<std::decay_t<decltype(arg)>, const char*>) return std::string(arg);
 *     else return std::string(arg);
 *   };
 *   return to_s("besian") + to_s("sejdiu") + to_s("@gmail.com");
 * }();
 *
 */

#ifndef LAHUTA_KD_TREE_TPP
#define LAHUTA_KD_TREE_TPP

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace lahuta {

template <typename Scalar>
void KDTree3<Scalar>::clear() {
  coords_vec_ = nullptr;
  coords_ptr_ = nullptr;
  n_points_ = 0;
  nodes_.clear();
  node_bounds_.clear();
  indices_.clear();
  node_count_ = 0;
  stack_.clear();
}

template <typename Scalar>
void KDTree3<Scalar>::build(const std::vector<Scalar> &coords, int leaf_size) {
  clear();
  coords_vec_ = &coords;
  coords_ptr_ = coords.data();
  build_from_pointer(coords.size() / DIM, leaf_size);
}

template <typename Scalar>
void KDTree3<Scalar>::build_view(const Scalar *coords_ptr, std::size_t n_points, int leaf_size) {
  clear();
  coords_vec_ = nullptr;
  coords_ptr_ = coords_ptr;
  build_from_pointer(n_points, leaf_size);
}

template <typename Scalar>
void KDTree3<Scalar>::build_from_pointer(std::size_t n_points, int leaf_size) {
  n_points_ = n_points;
  leaf_size_ = std::max(1, leaf_size);
  if (n_points_ == 0 || coords_ptr_ == nullptr) return;

  indices_.resize(n_points_);
  std::iota(indices_.begin(), indices_.end(), 0);

  int levels = 1;
  if (n_points_ > 1) {
    const double ratio = std::max(1.0, static_cast<double>(n_points_ - 1) / static_cast<double>(leaf_size_));
    levels += static_cast<int>(std::floor(std::log2(ratio)));
  }
  const int max_levels = static_cast<int>(sizeof(std::size_t) * 8 - 1);
  if (levels > max_levels) levels = max_levels;
  const std::size_t node_capacity = (static_cast<std::size_t>(1) << levels) - 1;

  nodes_.assign(node_capacity, {});
  node_bounds_.assign(node_capacity * 2 * DIM, Scalar(0));
  node_count_ = 0;

  build_node(0, 0, static_cast<int>(n_points_));
}

template <typename Scalar>
inline Scalar KDTree3<Scalar>::dist_sq(const Scalar *a, const Scalar *b) {
  const Scalar dx = a[0] - b[0];
  const Scalar dy = a[1] - b[1];
  const Scalar dz = a[2] - b[2];
  return dx * dx + dy * dy + dz * dz;
}

template <typename Scalar>
Scalar KDTree3<Scalar>::bbox_distance_sq(const Scalar *point, const Scalar *bbox_min, const Scalar *bbox_max) const {
  Scalar d2 = Scalar(0);
  for (int dim = 0; dim < DIM; ++dim) {
    const Scalar v = point[dim];
    if (v < bbox_min[dim]) {
      const Scalar diff = bbox_min[dim] - v;
      d2 += diff * diff;
    } else if (v > bbox_max[dim]) {
      const Scalar diff = v - bbox_max[dim];
      d2 += diff * diff;
    }
  }
  return d2;
}

template <typename Scalar>
void KDTree3<Scalar>::build_node(int node_index, int start, int end) {
  if (node_index >= static_cast<int>(nodes_.size()) || start >= end) return;

  node_count_ = std::max(node_count_, node_index + 1);

  Node &node = nodes_[static_cast<std::size_t>(node_index)];
  node.start = start;
  node.end   = end;

  Scalar *bbox_min = &node_bounds_[static_cast<std::size_t>(node_index) * 2 * DIM];
  Scalar *bbox_max = bbox_min + DIM;

  for (int dim = 0; dim < DIM; ++dim) {
    bbox_min[dim] =  std::numeric_limits<Scalar>::infinity();
    bbox_max[dim] = -std::numeric_limits<Scalar>::infinity();
  }

  for (int idx = start; idx < end; ++idx) {
    const Scalar *point = coords_ptr_ + DIM * indices_[static_cast<std::size_t>(idx)];
    for (int dim = 0; dim < DIM; ++dim) {
      const Scalar value = point[dim];
      if (value < bbox_min[dim]) bbox_min[dim] = value;
      if (value > bbox_max[dim]) bbox_max[dim] = value;
    }
  }

  const int count = end - start;
  const int left_child  = 2 * node_index + 1;
  const int right_child = left_child + 1;

  if (count <= leaf_size_ || left_child >= static_cast<int>(nodes_.size())) {
    node.leaf = true;
    node.split_dim = 0;
    node.split_value = Scalar(0);
    return;
  }

  int axis = 0;
  Scalar best_range = bbox_max[0] - bbox_min[0];
  for (int dim = 1; dim < DIM; ++dim) {
    const Scalar range = bbox_max[dim] - bbox_min[dim];
    if (range > best_range) {
      best_range = range;
      axis = dim;
    }
  }

  const int mid = start + count / 2;
  if (mid <= start || mid >= end) {
    node.leaf = true;
    node.split_dim = 0;
    node.split_value = Scalar(0);
    return;
  }

  auto begin  = indices_.begin() + start;
  auto middle = indices_.begin() + mid;
  auto end_it = indices_.begin() + end;
  std::nth_element(begin, middle, end_it, [this, axis](int lhs, int rhs) {
    const Scalar lv = coords_ptr_[DIM * lhs + axis];
    const Scalar rv = coords_ptr_[DIM * rhs + axis];
    return (lv == rv) ? (lhs < rhs) : (lv < rv);
  });

  node.leaf = false;
  node.split_dim = static_cast<std::uint8_t>(axis);
  node.split_value = coords_ptr_[DIM * indices_[static_cast<std::size_t>(mid)] + axis];

  build_node(left_child, start, mid);
  build_node(right_child, mid, end);
}

template <typename Scalar>
template <typename Callback>
void KDTree3<Scalar>::radius_search(const Scalar *point, Scalar radius_sq, Callback &&cb) const {
  if (!ready()) return;
  stack_.clear();
  stack_.push_back(0);

  while (!stack_.empty()) {
    const int node_idx = stack_.back();
    stack_.pop_back();
    if (node_idx < 0 || node_idx >= node_count_) continue;

    const Node &node = nodes_[static_cast<std::size_t>(node_idx)];
    if (node.end <= node.start) continue;

    const Scalar *bbox_min = &node_bounds_[static_cast<std::size_t>(node_idx) * 2 * DIM];
    const Scalar *bbox_max = bbox_min + DIM;
    if (bbox_distance_sq(point, bbox_min, bbox_max) > radius_sq) continue;

    if (node.leaf) {
      for (int i = node.start; i < node.end; ++i) {
        const int data_index = indices_[static_cast<std::size_t>(i)];
        const Scalar d2 = dist_sq(point, coords_ptr_ + DIM * data_index);
        if (d2 <= radius_sq) cb(data_index, d2);
      }
      continue;
    }

    const int left_child  = 2 * node_idx + 1;
    const int right_child = left_child + 1;

    const Scalar diff = point[node.split_dim] - node.split_value;
    int near_child = left_child;
    int far_child  = right_child;

    if (diff > Scalar(0)) {
      std::swap(near_child, far_child);
    }

    if (near_child < node_count_) stack_.push_back(near_child);

    if ( far_child < node_count_) {
      const Scalar plane_dist_sq = diff * diff;

      if (plane_dist_sq <= radius_sq) {
        stack_.push_back(far_child);
      } else {
        const Scalar *far_min = &node_bounds_[static_cast<std::size_t>(far_child) * 2 * DIM];
        const Scalar *far_max = far_min + DIM;
        if (bbox_distance_sq(point, far_min, far_max) <= radius_sq) stack_.push_back(far_child);
      }

    }
  }
}

} // namespace lahuta

#endif // LAHUTA_KD_TREE_TPP
