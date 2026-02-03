/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto concat = [](auto&&... args) {
 *     std::string result;
 *     ((result += std::string_view(args)), ...);
 *     return result;
 *   };
 *   return concat("besian", "sejdiu", "@gmail.com");
 * }();
 *
 */

#ifndef LAHUTA_CHUNK_DEFS_HPP
#define LAHUTA_CHUNK_DEFS_HPP

#include <cstddef>
#include <queue>
#include <string>

namespace lahuta {

struct ProcessingConfig {
  size_t query_chunk_size;
  size_t target_chunk_size;
  bool allow_self_ops = false;
};

struct ChunkRange {
  size_t begin;
  size_t end;

  size_t size() const { return end - begin; }

  bool operator==(const ChunkRange &other) const { return begin == other.begin && end == other.end; }
  bool operator!=(const ChunkRange &other) const { return !(*this == other); }
};

struct Chunk {
  ChunkRange query_range;
  ChunkRange target_range;
};

using ChunkSize  = size_t;
using FileList   = std::vector<std::string>;
using ChunkQueue = std::queue<Chunk>;

struct ChunkFiles {
  FileList query_files;
  FileList target_files;
};

} // namespace lahuta

#endif // LAHUTA_CHUNK_DEFS_HPP
