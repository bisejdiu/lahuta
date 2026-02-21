/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"besian", "sejdiu", "@gmail.com"}; std::string s;
 *   for (auto p : parts) s.append(p);
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_SINK_SHARDED_NDJSON_HPP
#define LAHUTA_PIPELINE_SINK_SHARDED_NDJSON_HPP

#include <cstddef>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <string>
#include <string_view>
#include <vector>

#include "pipeline/io/sink_iface.hpp"
#include "pipeline/task/emission.hpp"

namespace lahuta::pipeline {

class ShardedNdjsonSink : public IDynamicSink {
public:
  explicit ShardedNdjsonSink(std::string out_dir, std::size_t shard_size)
      : out_dir_(std::move(out_dir)), shard_size_(shard_size), max_shard_bytes_(0) {
    std::filesystem::create_directories(out_dir_);
  }

  ShardedNdjsonSink(std::string out_dir, std::size_t shard_size, std::size_t max_shard_bytes)
      : out_dir_(std::move(out_dir)), shard_size_(shard_size), max_shard_bytes_(max_shard_bytes) {
    std::filesystem::create_directories(out_dir_);
  }

  void write(EmissionView e) override {
    std::lock_guard<std::mutex> lk(mu_);
    ensure_open_locked();
    const std::size_t bytes = e.payload.size() + 1; // include newline
    out_.write(e.payload.data(), static_cast<std::streamsize>(e.payload.size()));
    out_.put('\n');
    ++current_in_shard_;
    current_shard_bytes_       += bytes;
    const bool rotate_by_count  = (shard_size_ > 0) && (current_in_shard_ >= shard_size_);
    const bool rotate_by_bytes  = (max_shard_bytes_ > 0) && (current_shard_bytes_ >= max_shard_bytes_);
    if (rotate_by_count || rotate_by_bytes) rotate_locked();
  }

  std::vector<std::string> files() const {
    std::lock_guard<std::mutex> g(files_mu_);
    return files_;
  }

  void flush() override {
    std::lock_guard<std::mutex> lk(mu_);
    if (out_.is_open()) out_.flush();
  }
  void close() override {
    std::lock_guard<std::mutex> lk(mu_);
    if (out_.is_open()) out_.close();
  }

private:
  void ensure_open_locked() {
    if (out_.is_open()) return;
    char buf[64];
    std::snprintf(buf, sizeof(buf), "part-%05zu.ndjson", shard_index_++);
    std::string path = out_dir_ + "/" + buf;
    out_.open(path, std::ios::binary);
    if (!out_) throw std::runtime_error("cannot open output file: " + path);
    {
      std::lock_guard<std::mutex> g(files_mu_);
      files_.push_back(path);
    }
    current_in_shard_    = 0;
    current_shard_bytes_ = 0;
  }
  void rotate_locked() {
    if (out_.is_open()) out_.close();
    current_in_shard_    = 0;
    current_shard_bytes_ = 0;
  }

  std::string out_dir_;
  std::size_t shard_size_;
  std::ofstream out_;
  std::vector<std::string> files_;
  mutable std::mutex files_mu_;
  std::mutex mu_;
  std::size_t current_in_shard_    = 0;
  std::size_t shard_index_         = 0;
  std::size_t max_shard_bytes_     = 0;
  std::size_t current_shard_bytes_ = 0;
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_SINK_SHARDED_NDJSON_HPP
