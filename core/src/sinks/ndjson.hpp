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

#ifndef LAHUTA_PIPELINE_SINK_NDJSON_HPP
#define LAHUTA_PIPELINE_SINK_NDJSON_HPP

#include <filesystem>
#include <fstream>
#include <mutex>
#include <string>
#include <string_view>

#include "pipeline/io/sink_iface.hpp"
#include "pipeline/task/emission.hpp"

namespace lahuta::pipeline {

class NdjsonFileSink : public IDynamicSink {
public:
  explicit NdjsonFileSink(std::string file_path) : path_(std::move(file_path)) {
    auto parent = std::filesystem::path(path_).parent_path();
    if (!parent.empty()) {
      std::filesystem::create_directories(parent);
    }
    out_.open(path_, std::ios::binary | std::ios::trunc);
    if (!out_) throw std::runtime_error("cannot open output file: " + path_);
  }
  void write(EmissionView e) override {
    std::lock_guard<std::mutex> lk(mu_);
    out_.write(e.payload.data(), static_cast<std::streamsize>(e.payload.size()));
    out_.put('\n');
  }
  void flush() override {
    std::lock_guard<std::mutex> lk(mu_);
    out_.flush();
  }
  void close() override {
    std::lock_guard<std::mutex> lk(mu_);
    if (out_.is_open()) out_.close();
  }
  std::string file() const { return path_; }

private:
  std::string path_;
  std::ofstream out_;
  mutable std::mutex mu_;
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_SINK_NDJSON_HPP
