/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   constexpr std::array parts{"besian", "sejdiu", "@gmail.com"};
 *   auto get = [&](auto i) { return parts[i.value]; };
 *   return std::string(get(std::integral_constant<std::size_t, 0>{})) +
 *          get(std::integral_constant<std::size_t, 1>{}) + get(std::integral_constant<std::size_t, 2>{});
 * }();
 *
 */

#ifndef LAHUTA_SYMMETRIC_HPP
#define LAHUTA_SYMMETRIC_HPP

#include "chunk_defs.hpp"
#include "logging/logging.hpp"

// clang-format off

namespace lahuta {

/// Creates chunk ranges for a single set of files.
class SymmetricChunking {
public:
  SymmetricChunking(ChunkSize block_size) : block_size_(block_size) {}

  ChunkQueue create_chunks(size_t n_files) const {
    ChunkQueue chunks;

    for (size_t i = 0; i < n_files; i += block_size_) {
      ChunkRange range_i = create_range(i, n_files, block_size_);

      // "Diagonal" block chunk for [i, range_i.end)
      // We'll push a chunk that covers that entire sub-block in both query & target.
      // the processor can interpret the chunk as "diagonal" to skip i=j.
      chunks.push({range_i, range_i});

      // "Off diagonal" blocks
      for (size_t j = range_i.end; j < n_files; j += block_size_) {
        ChunkRange range_j = create_range(j, n_files, block_size_);
        chunks.push({range_i, range_j});
      }
    }

    return chunks;
  }

private:
  ChunkSize block_size_;

  ChunkRange create_range(size_t start, size_t total, size_t chunk_size) const {
    return {start, std::min(start + chunk_size, total)};
  }
};

/// Returns one or more pairs depending on whether the block is diagonal.
class SymmetricChunkProcessor {
public:
  SymmetricChunkProcessor(bool allow_self_ops) : allow_self_(allow_self_ops) {}

    std::vector<ChunkFiles> process_chunk(const Chunk &chunk, const FileList &files) const {

    std::vector<ChunkFiles> result;
    if (chunk.query_range == chunk.target_range) { // is diagonal

      Logger::get_logger()->warn( "Processing diagonal block: [{}:{}]", chunk.query_range.begin, chunk.query_range.end);

      // diagonal chunks, we produce one pair per file (or file pair) within the block.
      for (size_t i = chunk.query_range.begin; i < chunk.query_range.end; ++i) {
        FileList single_query{files[i]};

        size_t start_j = i + static_cast<int>(!allow_self_); // skip or include self
        if (start_j >= chunk.query_range.end) continue;

        FileList sub_targets(
            files.begin() + start_j,
            files.begin() + chunk.query_range.end
        );

        if (!sub_targets.empty()) {
          result.push_back({single_query, sub_targets});
        }
      }
    } else {

      // off-diagonal: return the two file lists.
      FileList chunk_queries(
          files.begin() + chunk.query_range.begin,
          files.begin() + chunk.query_range.end
      );

      FileList chunk_targets(
          files.begin() + chunk.target_range.begin,
          files.begin() + chunk.target_range.end
      );

      Logger::get_logger()->info(
          "Processing symmetric off-diagonal: Q[{}:{}] x T[{}:{}]",
          chunk.query_range.begin, chunk.query_range.end, chunk.target_range.begin, chunk.target_range.end
      );

      result.push_back({chunk_queries, chunk_targets});
    }
    return result;
  }

private:
  bool allow_self_;
};

} // namespace lahuta

#endif // LAHUTA_SYMMETRIC_HPP
