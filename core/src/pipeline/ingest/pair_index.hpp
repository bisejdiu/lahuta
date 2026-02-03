/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
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

#ifndef LAHUTA_PIPELINE_INGEST_PAIR_INDEX_HPP
#define LAHUTA_PIPELINE_INGEST_PAIR_INDEX_HPP

#include <string>
#include <vector>

#include "fseek/seq.hpp"

namespace lahuta::pipeline {

struct PairRef {
  const SeqData *a;
  const SeqData *b;
};

class PairIndex {
public:
  using value_type = PairRef;

  PairIndex(const SeqCollection &s, bool include_same_file = false)
      : seqs_(s.get_data()), allow_same_file_(include_same_file) {}

  std::optional<value_type> next() {
    while (i_ < seqs_.size()) {
      if (++j_ >= seqs_.size()) {
        ++i_;
        j_ = i_;
        continue;
      }
      if (!allow_same_file_ && seqs_[i_].file_name == seqs_[j_].file_name) continue;
      return value_type{&seqs_[i_], &seqs_[j_]};
    }
    return std::nullopt;
  }

  void reset() {
    i_ = 0;
    j_ = 0;
  }

private:
  const std::vector<SeqData> &seqs_;
  std::size_t i_ = 0, j_ = 0;
  bool allow_same_file_;
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_INGEST_PAIR_INDEX_HPP
