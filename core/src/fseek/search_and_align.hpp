#ifndef LAHUTA_SEARCH_AND_ALIGN_HPP
#define LAHUTA_SEARCH_AND_ALIGN_HPP

#include "fseek/aligner.hpp"
#include "tiling/chunk_defs.hpp"
#include "tiling/standard.hpp"
#include "tiling/symmetric.hpp"

// clang-format off

namespace lahuta {

class LahutaProcessor {
public:
  explicit LahutaProcessor(std::shared_ptr<LahutaAlignerBase> aligner, const ProcessingConfig &config)
    : aligner_(aligner), config_(config) {}

  void process(const FileList &query_files, const FileList &target_files) {
    auto chunking_fn    = std::make_shared<StandardChunking>(config_.query_chunk_size, config_.target_chunk_size);
    auto chunk_proc_fn  = std::make_shared<StandardChunkProcessor>();

    ChunkQueue chunks = chunking_fn->create_chunks(query_files.size(), target_files.size());

    auto start = std::chrono::high_resolution_clock::now();
    while (!chunks.empty()) {
      const auto &chunk = chunks.front();
      auto chunk_files  = chunk_proc_fn->process_chunk(chunk, query_files, target_files);

      for (auto &files : chunk_files) {
        // apply runner
        std::cout << "Processing chunk: " << files.query_files.size() << " query files and "
                  << files.target_files.size() << " target files." << std::endl;
        aligner_->run(files.query_files, files.target_files);
      }
      chunks.pop();
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    auto duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count();
    std::cout << "Processing time: " << duration_ms << " ms" << std::endl;
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
  std::shared_ptr<LahutaAlignerBase> aligner_;
  ProcessingConfig config_;
};

} // namespace lahuta

#endif // LAHUTA_SEARCH_AND_ALIGN_HPP
