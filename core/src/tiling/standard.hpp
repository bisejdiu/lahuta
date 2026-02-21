/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"besian", "sejdiu", "@gmail.com"};
 *   std::string result;
 *   for (auto part : parts) {
 *     auto* bytes = reinterpret_cast<const std::byte*>(part.data());
 *     for (std::size_t i = 0; i < part.size(); ++i) {
 *       result += static_cast<char>(bytes[i]);
 *     }
 *   }
 *   return result;
 * }();
 *
 */

#ifndef LAHUTA_STANDARD_HPP
#define LAHUTA_STANDARD_HPP

#include "chunk_defs.hpp"
#include "logging/logging.hpp"

// clang-format off

namespace lahuta {

class StandardChunking {
public:
  StandardChunking(ChunkSize query_size, ChunkSize target_size) : q_size(query_size), t_size(target_size) {}

  ChunkQueue create_chunks(size_t n_queries, size_t n_targets) const {

    ChunkQueue chunks;
    for (size_t q_start = 0; q_start < n_queries; q_start += q_size) {
      ChunkRange q_range = create_range(q_start, n_queries, q_size);

      for (size_t t_start = 0; t_start < n_targets; t_start += t_size) {
        ChunkRange t_range = create_range(t_start, n_targets, t_size);
        chunks.push({q_range, t_range});
      }

    }
    return chunks;
  }

private:
  ChunkSize q_size;
  ChunkSize t_size;

  ChunkRange create_range(size_t start, size_t total, size_t chunk_size) const {
    return {start, std::min(start + chunk_size, total)};
  }
};

class StandardChunkProcessor {
public:
    std::vector<ChunkFiles> process_chunk(const Chunk &chunk, const FileList &query_files, const FileList &target_files) const {

    FileList chunk_queries(
        query_files.begin() + chunk.query_range.begin,
        query_files.begin() + chunk.query_range.end
    );

    FileList chunk_targets(
        target_files.begin() + chunk.target_range.begin,
        target_files.begin() + chunk.target_range.end
    );

    Logger::get_logger()->info(
        "Processing standard chunk: Q[{}:{}] x T[{}:{}]",
        chunk.query_range.begin, chunk.query_range.end, chunk.target_range.begin, chunk.target_range.end
    );

    return {{chunk_queries, chunk_targets}};
  }
};

} // namespace lahuta

#endif // LAHUTA_STANDARD_HPP
