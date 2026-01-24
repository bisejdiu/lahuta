#ifndef LAHUTA_RADIX_SORT_HPP
#define LAHUTA_RADIX_SORT_HPP

/*
 * Stable Least Significant Bit (LSB) radix sort for sorting an array of key-value pairs. It is
 * specialized for sorting the MMseqs2 score matrix. It has not been tested for general use.
 */

#include <algorithm>
#include <cstddef>

namespace lahuta {

inline unsigned short extract_key(unsigned short value, int shift) { return (value >> shift) & 0xFF; }

using Pair = std::pair<unsigned short, unsigned int>;
inline void stable_radix_sort(Pair *arr, size_t size, Pair *buffer) {
  constexpr int bits = 16;
  constexpr int radix = 256; // 8 bits per pass
  constexpr int passes = bits / 8;
  unsigned int count[radix] = {};

  Pair *__restrict input = arr;
  Pair *__restrict output = buffer;

  for (int pass = 0; pass < passes; ++pass) {
    std::fill_n(count, radix, 0);

    for (size_t i = 0; i < size; ++i) {
      ++count[(input[i].first >> (pass * 8)) & 0xFF];
    }

    for (int i = 1; i < radix; ++i) {
      count[i] += count[i - 1];
    }

    for (ptrdiff_t i = size - 1; i >= 0; --i) {
      unsigned short key = extract_key(input[i].first, pass * 8);
      output[--count[key]] = input[i];
    }

    std::swap(input, output);
  }

  if (input != arr) {
    std::copy(buffer, buffer + size, arr);
  }

  std::reverse(arr, arr + size);
}

} // namespace lahuta

#endif // LAHUTA_RADIX_SORT_HPP
