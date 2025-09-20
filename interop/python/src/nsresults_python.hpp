#ifndef LAHUTA_INTEROP_NSRESULTS_PYTHON_HPP
#define LAHUTA_INTEROP_NSRESULTS_PYTHON_HPP

#include <algorithm>
#include <cmath>
#include <cstring>
#include <tuple>
#include <utility>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "nsgrid.hpp"

// clang-format off
namespace py = pybind11;

namespace lahuta::bindings::neighbors {

namespace detail {

using Pair = std::pair<int, int>;

// Counts & prefix sums
inline std::vector<int> counts_cross(const std::vector<Pair>& pairs, int n_queries) {
  std::vector<int> counts(static_cast<size_t>(n_queries), 0);
  for (const auto& p : pairs) counts[static_cast<size_t>(p.first)]++;
  return counts;
}

inline std::vector<int> counts_self(const std::vector<Pair>& pairs, int n) {
  std::vector<int> counts(static_cast<size_t>(n), 0);
  for (const auto& p : pairs) {
    counts[static_cast<size_t>(p.first)]++;
    counts[static_cast<size_t>(p.second)]++;
  }
  return counts;
}

inline std::vector<int> make_indptr(const std::vector<int>& counts) {
  std::vector<int> indptr(counts.size() + 1, 0);
  for (size_t i = 0; i < counts.size(); ++i) {
    indptr[i + 1] = indptr[i] + counts[i];
  }
  return indptr;
}

// CSR allocation, sorting
struct CSR {
  py::array_t<double> distances;  // could be empty if not requested
  py::array_t<int>    indices;
  py::array_t<int>    indptr_py;
  int*                idx_ptr = nullptr;
  double*             dst_ptr = nullptr;
};

inline CSR allocate_csr(const std::vector<int>& indptr, bool with_distance) {
  const int total = indptr.empty() ? 0 : indptr.back();

  py::array_t<int> indices_arr(total);
  int* indices = total ? static_cast<int*>(indices_arr.mutable_data()) : nullptr;

  py::array_t<double> distances_arr;
  double* dst = nullptr;
  if (with_distance) {
    distances_arr = py::array_t<double>(total);
    dst = total ? static_cast<double*>(distances_arr.mutable_data()) : nullptr;
  }

  py::array_t<int> indptr_arr(static_cast<py::ssize_t>(indptr.size()));
  if (!indptr.empty())
    std::memcpy(indptr_arr.mutable_data(), indptr.data(), indptr.size() * sizeof(int));

  return CSR{std::move(distances_arr), std::move(indices_arr), std::move(indptr_arr), indices, dst};
}

inline void sort_rows_csr(int n_rows, int* indices, double* dst, const std::vector<int>& indptr, bool by_distance) {
  if (!indices) return;
  for (int i = 0; i < n_rows; ++i) {
    const int start = indptr[static_cast<size_t>(i)];
    const int end   = indptr[static_cast<size_t>(i + 1)];
    if (end - start <= 1) continue;

    if (by_distance && dst) {
      std::vector<std::pair<double, int>> tmp(static_cast<size_t>(end - start));
      for (int t = 0; t < end - start; ++t) {
        tmp[static_cast<size_t>(t)] = {dst[start + t], indices[start + t]};
      }
      std::sort(tmp.begin(), tmp.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
      for (int t = 0; t < end - start; ++t) {
        dst[start + t]     = tmp[static_cast<size_t>(t)].first;
        indices[start + t] = tmp[static_cast<size_t>(t)].second;
      }
    } else {
      std::sort(indices + start, indices + end);
    }
  }
}

// CSR filling
template <class Distances>
inline void fill_csr_cross(const std::vector<Pair>& pairs, const Distances& d2,
                           const std::vector<int>& indptr, int* indices, double* dst) {
  const int n_rows = static_cast<int>(indptr.size()) - 1;
  std::vector<int> cursor(static_cast<size_t>(n_rows), 0);

  for (size_t k = 0; k < pairs.size(); ++k) {
    const int qi  = pairs[k].first;
    const int pos = indptr[static_cast<size_t>(qi)] + cursor[static_cast<size_t>(qi)]++;
    if (indices) indices[pos] = pairs[k].second;
    if (dst)      dst[pos]    = std::sqrt(static_cast<double>(d2[k]));
  }
}

template <class Distances>
inline void fill_csr_self(const std::vector<Pair>& pairs, const Distances& d2,
                          const std::vector<int>& indptr, int* indices, double* dst) {
  const int n_rows = static_cast<int>(indptr.size()) - 1;
  std::vector<int> cursor(static_cast<size_t>(n_rows), 0);

  for (size_t k = 0; k < pairs.size(); ++k) {
    const int i = pairs[k].first;
    const int j = pairs[k].second;
    const double dist = std::sqrt(static_cast<double>(d2[k]));

    const int pos_i = indptr[static_cast<size_t>(i)] + cursor[static_cast<size_t>(i)]++;
    if (indices) indices[pos_i] = j;
    if (dst)      dst[pos_i]    = dist;

    const int pos_j = indptr[static_cast<size_t>(j)] + cursor[static_cast<size_t>(j)]++;
    if (indices) indices[pos_j] = i;
    if (dst)      dst[pos_j]    = dist;
  }
}

// ragged allocation, sorting
struct Ragged {
  std::vector<py::array_t<int>>    idx_arrays;
  std::vector<py::array_t<double>> dst_arrays; // may be empty (size n but default py::array_t when not requested)
  std::vector<int*>                idx_ptrs;
  std::vector<double*>             dst_ptrs;
};

inline Ragged allocate_ragged(const std::vector<int>& counts, bool with_distance) {
  Ragged r;
  const int n = static_cast<int>(counts.size());
  r.idx_arrays.resize(static_cast<size_t>(n));
  r.dst_arrays.resize(static_cast<size_t>(n));
  r.idx_ptrs  .assign(static_cast<size_t>(n), nullptr);
  r.dst_ptrs  .assign(static_cast<size_t>(n), nullptr);

  for (int i = 0; i < n; ++i) {
    const int c = counts[static_cast<size_t>(i)];
    r.idx_arrays[static_cast<size_t>(i)] = py::array_t<int>(static_cast<py::ssize_t>(c));
    r.idx_ptrs  [static_cast<size_t>(i)] = c ? static_cast<int*>(r.idx_arrays[static_cast<size_t>(i)].mutable_data()) : nullptr;

    if (with_distance) {
      r.dst_arrays[static_cast<size_t>(i)] = py::array_t<double>(static_cast<py::ssize_t>(c));
      r.dst_ptrs  [static_cast<size_t>(i)] = c ? static_cast<double*>(r.dst_arrays[static_cast<size_t>(i)].mutable_data()) : nullptr;
    }
  }
  return r;
}

inline void sort_ragged_row(int c, int* ip, double* dp, bool by_distance) {
  if (c <= 1 || !ip) return;
  if (by_distance && dp) {
    std::vector<std::pair<double, int>> tmp(static_cast<size_t>(c));
    for (int t = 0; t < c; ++t) {
      tmp[static_cast<size_t>(t)] = {dp[t], ip[t]};
    }
    std::sort(tmp.begin(), tmp.end(), [](const auto& a, const auto& b){ return a.first < b.first; });
    for (int t = 0; t < c; ++t) {
      dp[t] = tmp[static_cast<size_t>(t)].first; ip[t] = tmp[static_cast<size_t>(t)].second;
    }
  } else {
    std::sort(ip, ip + c);
  }
}

} // namespace detail

inline std::tuple<py::array_t<double>, py::array_t<int>, py::array_t<int>>
flatten_cross(const NSResults& res, int n_queries, bool return_distance, bool sort_results) {
  const auto& pairs = res.get_pairs();
  const auto& d2    = res.get_distances();

  auto counts = detail::counts_cross(pairs, n_queries);
  auto indptr = detail::make_indptr(counts);
  auto csr    = detail::allocate_csr(indptr, return_distance);

  detail::fill_csr_cross(pairs, d2, indptr, csr.idx_ptr, csr.dst_ptr);
  if (sort_results) detail::sort_rows_csr(n_queries, csr.idx_ptr, csr.dst_ptr, indptr, return_distance);

  return {csr.distances, csr.indices, csr.indptr_py};
}

inline std::tuple<py::array_t<double>, py::array_t<int>, py::array_t<int>>
flatten_self(const NSResults& res, int n, bool return_distance, bool sort_results) {
  const auto& pairs = res.get_pairs();
  const auto& d2    = res.get_distances();

  auto counts = detail::counts_self(pairs, n);
  auto indptr = detail::make_indptr(counts);
  auto csr    = detail::allocate_csr(indptr, return_distance);

  detail::fill_csr_self(pairs, d2, indptr, csr.idx_ptr, csr.dst_ptr);
  if (sort_results) detail::sort_rows_csr(n, csr.idx_ptr, csr.dst_ptr, indptr, return_distance);

  return {csr.distances, csr.indices, csr.indptr_py};
}

inline std::pair<py::list, py::list>
build_ragged_cross(const NSResults& res, int n_queries, bool return_distance, bool sort_results) {
  const auto& pairs = res.get_pairs();
  const auto& d2    = res.get_distances();

  auto counts = detail::counts_cross(pairs, n_queries);
  auto rag    = detail::allocate_ragged(counts, return_distance);

  std::vector<int> cursor(static_cast<size_t>(n_queries), 0);
  for (size_t k = 0; k < pairs.size(); ++k) {
    const int qi  = pairs[k].first;
    const int pos = cursor[static_cast<size_t>(qi)]++;
    if (auto* ip = rag.idx_ptrs[static_cast<size_t>(qi)]) ip[pos] = pairs[k].second;
    if (return_distance) {
      if (auto* dp = rag.dst_ptrs[static_cast<size_t>(qi)]) dp[pos] = std::sqrt(static_cast<double>(d2[k]));
    }
  }

  py::list indices_list(n_queries);
  py::list distances_list(n_queries);

  for (int i = 0; i < n_queries; ++i) {
    if (sort_results) {
      detail::sort_ragged_row(counts[static_cast<size_t>(i)],
                              rag.idx_ptrs[static_cast<size_t>(i)],
                              rag.dst_ptrs[static_cast<size_t>(i)],
                              return_distance);
    }
    indices_list[i] = rag.idx_arrays[static_cast<size_t>(i)];
    if (return_distance) {
      distances_list[i] = rag.dst_arrays[static_cast<size_t>(i)];
    }
  }

  return {indices_list, distances_list};
}

inline std::pair<py::list, py::list>
build_ragged_self(const NSResults& res, int n, bool return_distance, bool sort_results) {
  const auto& pairs = res.get_pairs();
  const auto& d2    = res.get_distances();

  auto counts = detail::counts_self(pairs, n);
  auto rag    = detail::allocate_ragged(counts, return_distance);

  std::vector<int> cursor(static_cast<size_t>(n), 0);
  for (size_t k = 0; k < pairs.size(); ++k) {
    const int i = pairs[k].first;
    const int j = pairs[k].second;
    const double dist = std::sqrt(static_cast<double>(d2[k]));

    int pos_i = cursor[static_cast<size_t>(i)]++;
    if (auto* ip = rag.idx_ptrs[static_cast<size_t>(i)]) ip[pos_i] = j;
    if (return_distance) {
      if (auto* dp = rag.dst_ptrs[static_cast<size_t>(i)]) dp[pos_i] = dist;
    }

    int pos_j = cursor[static_cast<size_t>(j)]++;
    if (auto* ip = rag.idx_ptrs[static_cast<size_t>(j)]) ip[pos_j] = i;
    if (return_distance) {
      if (auto* dp = rag.dst_ptrs[static_cast<size_t>(j)]) dp[pos_j] = dist;
    }
  }

  py::list indices_list(n);
  py::list distances_list(n);

  for (int i = 0; i < n; ++i) {
    if (sort_results) {
      detail::sort_ragged_row(counts[static_cast<size_t>(i)],
                              rag.idx_ptrs[static_cast<size_t>(i)],
                              rag.dst_ptrs[static_cast<size_t>(i)],
                              return_distance);
    }
    indices_list[i] = rag.idx_arrays[static_cast<size_t>(i)];
    if (return_distance) {
      distances_list[i] = rag.dst_arrays[static_cast<size_t>(i)];
    }
  }

  return {indices_list, distances_list};
}

} // namespace lahuta::bindings::neighbors

#endif // LAHUTA_INTEROP_NSRESULTS_PYTHON_HPP
