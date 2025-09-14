#ifndef LAHUTA_PIPELINE_DYNAMIC_SINK_SHARDED_NDJSON_HPP
#define LAHUTA_PIPELINE_DYNAMIC_SINK_SHARDED_NDJSON_HPP

#include <cstddef>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <string>
#include <string_view>
#include <vector>

#include "pipeline/dynamic/sink_iface.hpp"
#include "pipeline/dynamic/types.hpp"

// clang-format off
namespace lahuta::pipeline::dynamic {

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
    ensure_open();
    const std::size_t bytes = e.payload.size() + 1; // include newline
    out_.write(e.payload.data(), static_cast<std::streamsize>(e.payload.size()));
    out_.put('\n');
    ++current_in_shard_;
    current_shard_bytes_ += bytes;
    const bool rotate_by_count = (shard_size_ > 0) && (current_in_shard_ >= shard_size_);
    const bool rotate_by_bytes = (max_shard_bytes_ > 0) && (current_shard_bytes_ >= max_shard_bytes_);
    if (rotate_by_count || rotate_by_bytes) rotate();
  }

  std::vector<std::string> files() const {
    std::lock_guard<std::mutex> g(files_mu_);
    return files_;
  }

  void flush() override { if (out_.is_open()) out_.flush(); }
  void close() override { if (out_.is_open()) out_.close(); }

private:
  void ensure_open() {
    if (out_.is_open()) return;
    char buf[64];
    std::snprintf(buf, sizeof(buf), "part-%05zu.ndjson", shard_index_++);
    std::string path = out_dir_ + "/" + buf;
    out_.open(path, std::ios::binary);
    if (!out_) throw std::runtime_error("cannot open output file: " + path);
    { std::lock_guard<std::mutex> g(files_mu_); files_.push_back(path); }
    current_in_shard_    = 0;
    current_shard_bytes_ = 0;
  }
  void rotate() {
    if (out_.is_open()) out_.close();
    current_in_shard_    = 0;
    current_shard_bytes_ = 0;
  }

  std::string out_dir_;
  std::size_t shard_size_;
  std::ofstream out_;
  std::vector<std::string> files_;
  mutable std::mutex files_mu_;
  std::size_t current_in_shard_    = 0;
  std::size_t shard_index_         = 0;
  std::size_t max_shard_bytes_     = 0;
  std::size_t current_shard_bytes_ = 0;
};

} // namespace lahuta::pipeline::dynamic

#endif // LAHUTA_PIPELINE_DYNAMIC_SINK_SHARDED_NDJSON_HPP
