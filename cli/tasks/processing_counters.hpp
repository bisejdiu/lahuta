#ifndef LAHUTA_CLI_TASKS_PROCESSING_COUNTERS_HPP
#define LAHUTA_CLI_TASKS_PROCESSING_COUNTERS_HPP

#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <mutex>
#include <unordered_map>
#include <vector>

namespace lahuta::cli {

template <std::size_t ThresholdCount>
class ProcessingCounters {
public:
  struct Shard {
    std::uint64_t missing_sequence{0};
    std::uint64_t missing_dssp{0};
    std::uint64_t missing_positions{0};
    std::uint64_t length_mismatch{0};
    std::uint64_t atom_mismatch{0};
    std::uint64_t trimmed_empty{0};
    std::uint64_t prefix_errors{0};
    std::array<std::uint64_t, ThresholdCount> pass_by_threshold{};
    std::uint64_t total_kept_residues{0};
    std::uint64_t total_removed_residues{0};
    std::uint64_t total_trimmed_n{0};
    std::uint64_t total_trimmed_c{0};

    void add_pass(std::size_t index) {
      if (index < pass_by_threshold.size()) {
        ++pass_by_threshold[index];
      }
    }
  };

  struct Totals {
    std::uint64_t processed{0};
    std::uint64_t written{0};
    std::uint64_t missing_sequence{0};
    std::uint64_t missing_dssp{0};
    std::uint64_t missing_positions{0};
    std::uint64_t length_mismatch{0};
    std::uint64_t atom_mismatch{0};
    std::uint64_t trimmed_empty{0};
    std::uint64_t failed_confidence{0};
    std::uint64_t prefix_errors{0};
    std::array<std::uint64_t, ThresholdCount> pass_by_threshold{};
    std::uint64_t total_kept_residues{0};
    std::uint64_t total_removed_residues{0};
    std::uint64_t total_trimmed_n{0};
    std::uint64_t total_trimmed_c{0};

    void merge_from(const Shard &shard) {
      missing_sequence  += shard.missing_sequence;
      missing_dssp      += shard.missing_dssp;
      missing_positions += shard.missing_positions;
      length_mismatch   += shard.length_mismatch;
      atom_mismatch     += shard.atom_mismatch;
      trimmed_empty     += shard.trimmed_empty;
      prefix_errors     += shard.prefix_errors;
      for (std::size_t i = 0; i < ThresholdCount; ++i) {
        pass_by_threshold[i] += shard.pass_by_threshold[i];
      }
      total_kept_residues    += shard.total_kept_residues;
      total_removed_residues += shard.total_removed_residues;
      total_trimmed_n        += shard.total_trimmed_n;
      total_trimmed_c        += shard.total_trimmed_c;
    }
  };

  ProcessingCounters() : id_(next_id_.fetch_add(1, std::memory_order_relaxed)) {}

  std::uint64_t bump_processed() { return processed_.fetch_add(1, std::memory_order_relaxed) + 1; }

  void bump_written() { written_.fetch_add(1, std::memory_order_relaxed); }

  void bump_failed_confidence() { failed_confidence_.fetch_add(1, std::memory_order_relaxed); }

  std::uint64_t written_total() const { return written_.load(std::memory_order_relaxed); }

  std::uint64_t failed_confidence_total() const { return failed_confidence_.load(std::memory_order_relaxed); }

  Shard &local() {
    thread_local std::unordered_map<std::size_t, Shard *> shard_cache;
    auto it = shard_cache.find(id_);
    if (it != shard_cache.end() && it->second) {
      return *it->second;
    }
    Shard &shard = register_shard();
    shard_cache.emplace(id_, &shard);
    return shard;
  }

  Totals snapshot() const {
    Totals totals;
    totals.processed         = processed_.load(std::memory_order_relaxed);
    totals.written           = written_.load(std::memory_order_relaxed);
    totals.failed_confidence = failed_confidence_.load(std::memory_order_relaxed);

    std::lock_guard<std::mutex> lock(mutex_);
    for (const auto &shard : shards_) {
      if (shard) {
        totals.merge_from(*shard);
      }
    }
    return totals;
  }

private:
  Shard &register_shard() {
    std::lock_guard<std::mutex> lock(mutex_);
    shards_.push_back(std::make_unique<Shard>());
    return *shards_.back();
  }

  static inline std::atomic<std::size_t> next_id_{0};
  const std::size_t id_;
  mutable std::mutex mutex_;
  mutable std::vector<std::unique_ptr<Shard>> shards_;
  std::atomic<std::uint64_t> processed_{0};
  std::atomic<std::uint64_t> written_{0};
  std::atomic<std::uint64_t> failed_confidence_{0};
};

} // namespace lahuta::cli

#endif // LAHUTA_CLI_TASKS_PROCESSING_COUNTERS_HPP
