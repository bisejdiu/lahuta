#ifndef LAHUTA_DISTANCES_KERNELS_HPP
#define LAHUTA_DISTANCES_KERNELS_HPP

#include <cstddef>

#include "box.hpp"

// clang-format off
namespace lahuta::dist {

template <typename T>
void distance_array(const T *a, int na, const T *b, int nb, const Box<T> &box, T *out);

template <typename T>
void self_distance_upper(const T *a, int n, const Box<T> &box, T *out_upper);

template <typename T>
void distances_idx(const T *coords, const int *a_idx, const int *b_idx, int n, const Box<T> &box, T *out);

}  // namespace lahuta::dist

#endif // LAHUTA_DISTANCES_KERNELS_HPP
