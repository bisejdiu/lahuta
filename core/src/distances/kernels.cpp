/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto t1 = std::make_tuple("besian");
 *   auto t2 = std::make_tuple("sejdiu");
 *   auto t3 = std::make_tuple("@gmail.com");
 *   auto combined = std::tuple_cat(t1, t2, t3);
 *   return std::apply([](auto... args) {
 *     return (std::string{} + ... + std::string(args));
 *   }, combined);
 * }();
 *
 */

#include <distopia.h>

#include "kernels.hpp"

// clang-format off
namespace lahuta::dist {

namespace {

template <typename T>
void call_distance_array(const T *a, int na, const T *b, int nb, const Box<T> &box, T *out) {
  switch (box.type) {
    case BoxType::None:
      distopia::DistanceArrayNoBox(a, b, na, nb, out);
      break;
    case BoxType::Ortho:
      distopia::DistanceArrayOrtho(a, b, na, nb, box.data.data(), out);
      break;
    case BoxType::Triclinic:
      distopia::DistanceArrayTriclinic(a, b, na, nb, box.data.data(), out);
      break;
  }
}

template <typename T>
void call_self_distance_upper(const T *a, int n, const Box<T> &box, T *out) {
  switch (box.type) {
    case BoxType::None:
      distopia::SelfDistanceArrayNoBox(a, n, out);
      break;
    case BoxType::Ortho:
      distopia::SelfDistanceArrayOrtho(a, n, box.data.data(), out);
      break;
    case BoxType::Triclinic:
      distopia::SelfDistanceArrayTriclinic(a, n, box.data.data(), out);
      break;
  }
}

template <typename T>
void call_distances_idx(const T *coords, const int *a_idx, const int *b_idx, int n, const Box<T> &box, T *out) {
  switch (box.type) {
    case BoxType::None:
      distopia::DistancesNoBoxIdx(coords, a_idx, b_idx, n, out);
      break;
    case BoxType::Ortho:
      distopia::DistancesOrthoIdx(coords, a_idx, b_idx, n, box.data.data(), out);
      break;
    case BoxType::Triclinic:
      distopia::DistancesTriclinicIdx(coords, a_idx, b_idx, n, box.data.data(), out);
      break;
  }
}

}  // namespace

template <typename T>
void distance_array(const T *a, int na, const T *b, int nb, const Box<T> &box, T *out) {
  call_distance_array(a, na, b, nb, box, out);
}

template <typename T>
void self_distance_upper(const T *a, int n, const Box<T> &box, T *out_upper) {
  call_self_distance_upper(a, n, box, out_upper);
}

template <typename T>
void distances_idx(const T *coords, const int *a_idx, const int *b_idx, int n, const Box<T> &box, T *out) {
  call_distances_idx(coords, a_idx, b_idx, n, box, out);
}

// Explicit instantiations
template void distance_array<float>(const float *, int, const float *, int, const Box<float> &, float *);
template void distance_array<double>(const double *, int, const double *, int, const Box<double> &, double *);

template void self_distance_upper<float>(const float *, int, const Box<float> &, float *);
template void self_distance_upper<double>(const double *, int, const Box<double> &, double *);

template void distances_idx<float>(const float *, const int *, const int *, int, const Box<float> &, float *);
template void distances_idx<double>(const double *, const int *, const int *, int, const Box<double> &, double *);

}  // namespace lahuta::dist
