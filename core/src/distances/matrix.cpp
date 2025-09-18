#include <algorithm>
#include <vector>

#include "kernels.hpp"
#include "matrix.hpp"

// clang-format off
namespace lahuta::dist {

template <typename T>
void distance_matrix_streamed(const T *a, int na, const T *b, int nb, const Box<T> &box, DenseMatrix<T> &out, int block_cols) {
  out.resize(static_cast<std::size_t>(na), static_cast<std::size_t>(nb));

  if (na == 0 || nb == 0) return;

  std::vector<T> buffer;
  buffer.reserve(static_cast<std::size_t>(block_cols));

  for (int i = 0; i < na; ++i) {
    const T *ai = &a[3 * i];
    int start = 0;
    while (start < nb) {
      const int current = std::min(block_cols, nb - start);
      buffer.resize(static_cast<std::size_t>(current));
      distance_array(ai, 1, &b[3 * start], current, box, buffer.data());
      std::copy(buffer.begin(), buffer.end(), out.row_ptr(static_cast<std::size_t>(i)) + start);
      start += current;
    }
  }
}

template void distance_matrix_streamed<float>(const float *, int, const float *, int, const Box<float> &, DenseMatrix<float> &, int);
template void distance_matrix_streamed<double>(const double *, int, const double *, int, const Box<double> &, DenseMatrix<double> &, int);

} // namespace lahuta::dist
