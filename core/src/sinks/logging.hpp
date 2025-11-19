#ifndef LAHUTA_PIPELINE_DYNAMIC_SINK_LOGGING_HPP
#define LAHUTA_PIPELINE_DYNAMIC_SINK_LOGGING_HPP

#include <cstdio>
#include <mutex>
#include <string_view>

#include "pipeline/dynamic/sink_iface.hpp"
#include "pipeline/dynamic/types.hpp"

// clang-format off
namespace lahuta::pipeline::dynamic {

class LoggingSink : public IDynamicSink {
public:
  void write(EmissionView e) override {
    std::lock_guard<std::mutex> lk(mu_);
    fwrite(e.payload.data(), 1, e.payload.size(), stdout);
    fputc('\n', stdout);
    fflush(stdout);
  }
  void flush() override {
    std::lock_guard<std::mutex> lk(mu_);
    fflush(stdout);
  }
  void close() override {}

private:
  std::mutex mu_;
};

} // namespace lahuta::pipeline::dynamic

#endif // LAHUTA_PIPELINE_DYNAMIC_SINK_LOGGING_HPP
