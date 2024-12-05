#ifndef LAHUTA_PROCESSOR_HPP
#define LAHUTA_PROCESSOR_HPP

#include "aligner.hpp"
#include "file_system.hpp"
#include <memory>
#include <queue>
#include <string>
#include <vector>

namespace lahuta {

class LahutaProcessor {
public:
  struct ProcessingConfig {
    size_t query_chunk_size;
    size_t target_chunk_size;
    bool symmetric_mode = false;
    bool allow_query_self_operations = false;
    std::string file_extension;
    bool recursive_search = false;
  };

  struct ChunkRange {
    size_t begin;
    size_t end;

    size_t size() const { return end - begin; }
  };

  struct ProcessingChunk {
    ChunkRange query_range;
    ChunkRange target_range;
    bool is_diagonal_chunk = false;
  };

  using FileList = std::vector<std::string>;
  using ChunkQueue = std::queue<ProcessingChunk>; // can we make it lock free when parallelizing?

  LahutaProcessor(ProcessingConfig config) : cfg(std::move(config)) {}


  void process_files(const FileList &query_files, const FileList &target_files) {
    const FileList &target_files_ref = target_files.empty() ? query_files : target_files;

    if (cfg.symmetric_mode || &query_files == &target_files) {
      cfg.symmetric_mode = true;
      orchestrate_symmetric(query_files.size());
    } else {
      orchestrate(query_files.size(), target_files_ref.size());
    }

    /*orchestrate(query_files.size(), target_files_ref.size());*/
    execute(query_files, target_files_ref);
  }

  void process_files(const FileList &files) {
    process_files(files, files);
  }

  /// Takes ownership of the runner
  void set_runner(std::unique_ptr<Runner> runner) { this->runner_ = std::move(runner); }

private:
  void orchestrate_symmetric(size_t num_files) {
    for (size_t i = 0; i < num_files; i += cfg.query_chunk_size) {
      ChunkRange range_i{i, std::min(i + cfg.query_chunk_size, num_files)};

      // diagonal chunk (i,i)
      processing_queue_.push({range_i, range_i, true});

      // Processing only upper triangular chunks (i,j) where j > i
      for (size_t j = range_i.end; j < num_files; j += cfg.target_chunk_size) {
        ChunkRange range_j{j, std::min(j + cfg.target_chunk_size, num_files)};
        processing_queue_.push({range_i, range_j, false});
      }
    }
  }

  LahutaAligner get_default_runner() { return LahutaAligner(); }
  void orchestrate(size_t num_queries, size_t num_targets) {
    for (size_t q_start = 0; q_start < num_queries; q_start += cfg.query_chunk_size) {
      ChunkRange query_range{q_start, std::min(q_start + cfg.query_chunk_size, num_queries)};

      for (size_t t_start = 0; t_start < num_targets; t_start += cfg.target_chunk_size) {
        ChunkRange target_range{t_start, std::min(t_start + cfg.target_chunk_size, num_targets)};

        processing_queue_.push({query_range, target_range});
      }
    }
  }

  void execute(const FileList &query_files, const FileList &target_files) {
    while (!processing_queue_.empty()) {
      auto chunk = processing_queue_.front();
      processing_queue_.pop();

      process_chunk(chunk, query_files, target_files);
    }
  }

  void
  process_chunk(const ProcessingChunk &chunk, const FileList &query_files, const FileList &target_files) {
    // Extract file lists for current chunk
    FileList chunk_queries(
        query_files.begin() + chunk.query_range.begin,
        query_files.begin() + chunk.query_range.end);

    FileList chunk_targets(
        target_files.begin() + chunk.target_range.begin,
        target_files.begin() + chunk.target_range.end);

    spdlog::warn(
        "Processing chunk: Q[{}:{}] x T[{}:{}]",
        chunk.query_range.begin,
        chunk.query_range.end,
        chunk.target_range.begin,
        chunk.target_range.end);

    /*get_runner().run(chunk_queries, chunk_targets);*/

    if (cfg.symmetric_mode && chunk.is_diagonal_chunk) {
      // For diagonal chunks, only process pairs where j > i
      process_diagonal_chunk(chunk_queries);
    } else {
      get_runner().run(chunk_queries, chunk_targets);
    }
  }

  void process_diagonal_chunk(const FileList &files) {
    // Process only pairs where j > i within the diagonal chunk
    for (size_t i = 0; i < files.size(); ++i) {
      FileList single_query{files[i]};

      FileList remaining_targets(
          files.begin() + i + static_cast<int>(!cfg.allow_query_self_operations),
          files.end());

      if (!remaining_targets.empty()) {
        get_runner().run(single_query, remaining_targets);
      }
    }
  }

  std::unique_ptr<Runner> set_default_runner() { return std::make_unique<LahutaAligner>(); }
  Runner &get_runner() {
    if (!runner_) runner_ = set_default_runner();
    return *runner_;
  }

  ProcessingConfig cfg;
  ChunkQueue processing_queue_;
  std::unique_ptr<Runner> runner_;
};

} // namespace lahuta

#endif // LAHUTA_PROCESSOR_HPP
