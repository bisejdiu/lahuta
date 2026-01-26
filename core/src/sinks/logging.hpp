#ifndef LAHUTA_PIPELINE_SINK_LOGGING_HPP
#define LAHUTA_PIPELINE_SINK_LOGGING_HPP

#include <cstdio>
#include <mutex>
#include <string_view>

#include "pipeline/io/sink_iface.hpp"
#include "pipeline/task/emission.hpp"

namespace lahuta::pipeline {

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

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_SINK_LOGGING_HPP
