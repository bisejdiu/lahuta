#ifndef LAHUTA_PIPELINE_PAIR_INDEX_SOURCE_HPP
#define LAHUTA_PIPELINE_PAIR_INDEX_SOURCE_HPP

#include <seq.hpp>
#include <string>
#include <vector>

namespace lahuta::sources {

struct PairRef {
  const SeqData *a;
  const SeqData *b;
};

class PairIndexSource {
public:
  using value_type = PairRef;

  PairIndexSource(const SeqCollection &s, bool include_same_file = false)
      : seqs_(s.get_data()), allow_same_file_(include_same_file) {}

  std::optional<value_type> next() {
    while (i_ < seqs_.size()) {
      if (++j_ >= seqs_.size()) {
        ++i_;
        j_ = i_; // FIX: there is an important subtlety here. Document!
        continue;
      }
      if (!allow_same_file_ && seqs_[i_].file_name == seqs_[j_].file_name) continue;
      return value_type{&seqs_[i_], &seqs_[j_]};
    }
    return std::nullopt;
  }

private:
  const std::vector<SeqData> &seqs_;
  std::size_t i_ = 0, j_ = 0;
  bool allow_same_file_;
};

} // namespace lahuta::sources

#endif // LAHUTA_PIPELINE_PAIR_INDEX_SOURCE_HPP
