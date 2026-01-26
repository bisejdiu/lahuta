#ifndef LAHUTA_PIPELINE_SINK_MEMORY_HPP
#define LAHUTA_PIPELINE_SINK_MEMORY_HPP

#include <mutex>
#include <string>
#include <vector>

#include "pipeline/io/sink_iface.hpp"
#include "pipeline/task/emission.hpp"

namespace lahuta::pipeline {

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

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_SINK_MEMORY_HPP
