#ifndef LAHUTA_SEARCH_AND_ALIGN_HPP
#define LAHUTA_SEARCH_AND_ALIGN_HPP

#include "aligner.hpp"
#include "chunking/chunk_defs.hpp"
#include "chunking/standard.hpp"
#include "chunking/symmetric.hpp"

// clang-format off

namespace lahuta {

class LahutaProcessor {
public:
  explicit LahutaProcessor(std::shared_ptr<LahutaAligner> aligner, const ProcessingConfig &config)
    : aligner_(aligner), config_(config) {}

  void process(const FileList &query_files, const FileList &target_files) {
    auto chunking_fn    = std::make_shared<StandardChunking>(config_.query_chunk_size, config_.target_chunk_size);
    auto chunk_proc_fn  = std::make_shared<StandardChunkProcessor>();

    ChunkQueue chunks = chunking_fn->create_chunks(query_files.size(), target_files.size());

    while (!chunks.empty()) {
      const auto &chunk = chunks.front();
      auto chunk_files  = chunk_proc_fn->process_chunk(chunk, query_files, target_files);

      for (auto &files : chunk_files) {
        // apply runner
        aligner_->run(files.query_files, files.target_files);
      }
      chunks.pop();
    }
  }

  void process(const FileList &query_files) {
    auto chunking_fn    = std::make_shared<SymmetricChunking>(config_.query_chunk_size);
    auto chunk_proc_fn  = std::make_shared<SymmetricChunkProcessor>(config_.allow_self_ops);

    ChunkQueue chunks = chunking_fn->create_chunks(query_files.size());

    while (!chunks.empty()) {
      const auto &chunk = chunks.front();
      auto chunk_files  = chunk_proc_fn->process_chunk(chunk, query_files);

      for (auto &files : chunk_files) {
        // apply runner
        aligner_->run(files.query_files, files.target_files);
      }
      chunks.pop();
    }
  }

private:
  std::shared_ptr<LahutaAligner> aligner_;
  ProcessingConfig config_;
};

} // namespace lahuta

#endif // LAHUTA_SEARCH_AND_ALIGN_HPP

