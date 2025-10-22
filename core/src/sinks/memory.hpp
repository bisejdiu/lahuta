#ifndef LAHUTA_PIPELINE_DYNAMIC_SINK_MEMORY_HPP
#define LAHUTA_PIPELINE_DYNAMIC_SINK_MEMORY_HPP

#include <mutex>
#include <string>
#include <vector>

#include "pipeline/dynamic/sink_iface.hpp"
#include "pipeline/dynamic/types.hpp"

// clang-format off
namespace lahuta::pipeline::dynamic {

class MemorySink : public IDynamicSink {
public:
  void write(EmissionView e) override {
    std::lock_guard<std::mutex> lk(m_);
    out_.emplace_back(e.payload);
  }
  std::vector<std::string> result() const {
    std::lock_guard<std::mutex> lk(m_);
    return out_;
  }
  std::vector<std::string> result_bytes() const {
    std::lock_guard<std::mutex> lk(m_);
    return out_;
  }
  void clear() {
    std::lock_guard<std::mutex> lk(m_);
    out_.clear();
  }
  void flush() override {}
  void close() override {}

private:
  mutable std::mutex m_;
  std::vector<std::string> out_;
};

} // namespace lahuta::pipeline::dynamic

#endif // LAHUTA_PIPELINE_DYNAMIC_SINK_MEMORY_HPP
