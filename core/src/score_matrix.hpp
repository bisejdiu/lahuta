// implementation based on mmseqs/src/prefiltering/ExtendedSubstitutionMatrix.cpp
#ifndef LAHUTA_SCORE_MATRIX_HPP
#define LAHUTA_SCORE_MATRIX_HPP

#include <array>
#include <atomic>
#include <cstddef>
#include <thread>

#include "ScoreMatrix.h"
#include "ctpl/ctpl.h"
#include "radix_sort.hpp"
#include "substitution_matrix.hpp"
#include "utils.hpp"

namespace lahuta {

template <std::size_t AlphabetSize, std::size_t KmerSize> //
class ConstExprIndexer {
public:
  static constexpr std::size_t PermutationSize = __pow__(AlphabetSize, KmerSize);

  constexpr ConstExprIndexer(
      const std::array<std::array<unsigned char, KmerSize>, PermutationSize> &permutations)
      : index_table{} {
    for (std::size_t i = 0; i < PermutationSize; ++i) {
      index_table[i] = compute_index(permutations[i]);
    }
  }

  constexpr std::size_t get_index(std::size_t permutation_index) const {
    return index_table[permutation_index];
  }

private:
  std::array<std::size_t, PermutationSize> index_table;

  static constexpr std::size_t compute_index(const std::array<unsigned char, KmerSize> &kmer) {
    std::size_t index = 0;
    std::size_t power = 1;
    for (std::size_t i = 0; i < KmerSize; ++i) {
      index += kmer[i] * power;
      power *= AlphabetSize;
    }
    return index;
  }
};

template <size_t KmerSize, size_t AlphabetSize>
constexpr std::array<std::array<unsigned char, AlphabetSize>, KmerSize> build_input() {
  std::array<std::array<unsigned char, AlphabetSize>, KmerSize> dimension_vector{};

  for (size_t i = 0; i < KmerSize; i++) {
    for (size_t j = 0; j < AlphabetSize; j++) {
      dimension_vector[i][j] = static_cast<unsigned char>(j);
    }
  }
  return dimension_vector;
}

template <size_t KmerSize, size_t AlphabetSize>
constexpr void cartesian_product(
    const std::array<std::array<unsigned char, AlphabetSize>, KmerSize> &input,
    std::array<std::array<unsigned char, KmerSize>, __pow__(AlphabetSize, KmerSize)> &output,
    std::array<unsigned char, KmerSize> &current_result, size_t depth, size_t &result_index) {

  if (depth == KmerSize) {
    output[result_index++] = current_result;
    return;
  }

  for (const auto &element : input[depth]) {
    current_result[depth] = element;
    cartesian_product(input, output, current_result, depth + 1, result_index);
  }
}

template <size_t KmerSize, size_t AlphabetSize>
inline short calc_score(
    const std::array<unsigned char, KmerSize> &i_seq, const std::array<unsigned char, KmerSize> &j_seq,
    const std::array<std::array<short, AlphabetSize>, AlphabetSize> &subMatrix) {
  short score = 0;
  for (size_t i = 0; i < KmerSize; i++) {
    score += subMatrix[i_seq[i]][j_seq[i]];
  }
  return score;
}

template <size_t KmerSize, size_t AlphabetSize>
constexpr std::array<std::array<unsigned char, KmerSize>, __pow__(AlphabetSize, KmerSize)>
generate_permutations() {

  constexpr auto permutation_size = __pow__(AlphabetSize, KmerSize);

  std::array<std::array<unsigned char, KmerSize>, permutation_size> permutation{};
  std::array<unsigned char, KmerSize> outputTemp{};
  size_t result_index = 0;
  constexpr auto kmer_input = build_input<KmerSize, AlphabetSize>();

  cartesian_product<KmerSize, AlphabetSize>(kmer_input, permutation, outputTemp, 0, result_index);
  return permutation;
}

template <size_t KmerSize, size_t AlphabetSize>
ScoreMatrix calculate_score_matrix(const __SubMatrix__ &matrix) {
  constexpr static size_t size = __pow__(AlphabetSize, KmerSize);
  constexpr static size_t row_size = ((size / MAX_ALIGN_INT) + 1) * MAX_ALIGN_INT;
  constexpr static auto permutation = generate_permutations<KmerSize, AlphabetSize>();
  constexpr static ConstExprIndexer<AlphabetSize, KmerSize> indexer(permutation);

  auto *score = (short *)mem_align(MAX_ALIGN_INT, (size * (row_size)) * sizeof(short));
  auto *index = (unsigned int *)mem_align(MAX_ALIGN_INT, (size * (row_size)) * sizeof(unsigned int));

  // chunking would give an additional 3-5% speedup (not worth the complexity)
  std::atomic<size_t> current_index(0);
  auto worker = [&](int) {
    auto *score_matrix = new std::pair<unsigned short, unsigned int>[size];
    auto *buffer = new std::pair<unsigned short, unsigned int>[size];

    while (true) {
      size_t i = current_index.fetch_add(1);
      if (i >= permutation.size()) break;

      const unsigned int i_index = indexer.get_index(i);
      for (size_t j = 0; j < permutation.size(); j++) {
        const unsigned int j_index = indexer.get_index(j);
        const auto score_val = calc_score<>(permutation[i], permutation[j], matrix.subMatrix);
        unsigned short score_val_shifted = static_cast<unsigned>(score_val + 1E5);
        score_matrix[j] = {score_val_shifted, j_index};
      }

      stable_radix_sort(score_matrix, size, buffer);

      for (size_t z = 0; z < size; z++) {
        score[(i_index * row_size) + z] = static_cast<short>(score_matrix[z].first - 1E5);
        index[(i_index * row_size) + z] = score_matrix[z].second;
      }

      for (size_t z = size; z < row_size; z++) {
        score[(i_index * row_size) + z] = -255;
        index[(i_index * row_size) + z] = 0;
      }
    }

    delete[] score_matrix;
    delete[] buffer;
  };

  const size_t num_threads = std::thread::hardware_concurrency();
  ctpl::thread_pool pool(num_threads);
  std::vector<std::future<void>> futures;

  for (size_t t = 0; t < num_threads; t++) {
    futures.emplace_back(pool.push(worker));
  }

  for (auto &future : futures) {
    future.get();
  }

  return ScoreMatrix{score, index, size, row_size};
}

} // namespace lahuta

#endif // LAHUTA_SCORE_MATRIX_HPP
