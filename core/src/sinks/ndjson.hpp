#ifndef LAHUTA_PIPELINE_DYNAMIC_SINK_NDJSON_HPP
#define LAHUTA_PIPELINE_DYNAMIC_SINK_NDJSON_HPP

#include <filesystem>
#include <fstream>
#include <mutex>
#include <string>
#include <string_view>

#include "pipeline/dynamic/sink_iface.hpp"
#include "pipeline/dynamic/types.hpp"

// clang-format off
namespace lahuta::pipeline::dynamic {

class NdjsonFileSink : public IDynamicSink {
public:
  explicit NdjsonFileSink(std::string file_path)
    : path_(std::move(file_path)) {
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

} // namespace lahuta::pipeline::dynamic

#endif // LAHUTA_PIPELINE_DYNAMIC_SINK_NDJSON_HPP
