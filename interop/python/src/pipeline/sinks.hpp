#ifndef LAHUTA_BINDINGS_DYNAMIC_SINKS_HPP
#define LAHUTA_BINDINGS_DYNAMIC_SINKS_HPP

#include <chrono>
#include <cmath>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <thread>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "db/db.hpp"
#include "io/sinks/lmdb.hpp"
#include "io/sinks/memory.hpp"
#include "io/sinks/ndjson.hpp"
#include "io/sinks/sharded_ndjson.hpp"

namespace py = pybind11;
namespace lahuta::bindings {
using namespace lahuta::pipeline::dynamic;

// clang-format off
namespace {

inline std::chrono::milliseconds seconds_to_ms(double seconds) {
  using MilliRep = std::chrono::milliseconds::rep;

  if (!std::isfinite(seconds)) throw std::invalid_argument("SlowSink timing must be finite");
  if (seconds < 0.0) throw std::invalid_argument("SlowSink timing must be non-negative");

  const double max_seconds = static_cast<double>(std::numeric_limits<MilliRep>::max()) / 1000.0;
  if (seconds > max_seconds) throw std::invalid_argument("SlowSink timing is too large");

  return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::duration<double>(seconds));
}

} // namespace

// Test-only sink that sleeps during write() to simulate slow I/O
class SlowSink : public IDynamicSink {
public:
  explicit SlowSink(double sleep_seconds, double flush_sleep_seconds = 0.0)
      : write_sleep_(seconds_to_ms(sleep_seconds)), flush_sleep_(seconds_to_ms(flush_sleep_seconds)) {}

  void write(EmissionView e) override {
    (void)e;
    if (write_sleep_.count() > 0) {
      std::this_thread::sleep_for(write_sleep_);
    }
    std::lock_guard<std::mutex> lk(m_);
    count_++;
  }

  void flush() override {
    if (flush_sleep_.count() > 0) {
      std::this_thread::sleep_for(flush_sleep_);
    }
  }

  std::size_t count() const {
    std::lock_guard<std::mutex> lk(m_);
    return count_;
  }

private:
  mutable std::mutex m_;
  std::size_t count_ = 0;
  std::chrono::milliseconds write_sleep_;
  std::chrono::milliseconds flush_sleep_;
};

// clang-format off
inline void bind_sinks(py::module_& md) {
  py::class_<IDynamicSink, std::shared_ptr<IDynamicSink>> _(md, "Sink");

  py::class_<MemorySink, IDynamicSink, std::shared_ptr<MemorySink>>(md, "MemorySink")
    .def(py::init<>())
    .def("result", &MemorySink::result, R"doc(Return collected payloads.)doc")
    .def("result_bytes", [](const MemorySink& s) {
          py::list out;
          for (const auto& buf : s.result_bytes()) {
            out.append(py::bytes(buf));
          }
          return out;
        }, R"doc(Return collected payloads as bytes.)doc")
    .def("clear",  &MemorySink::clear,  R"doc(Clear collected payloads.)doc");

  py::class_<NdjsonFileSink, IDynamicSink, std::shared_ptr<NdjsonFileSink>>(md, "NdjsonSink")
    .def(py::init<std::string>(), py::arg("path"))
    .def("file", &NdjsonFileSink::file);

  py::class_<ShardedNdjsonSink, IDynamicSink, std::shared_ptr<ShardedNdjsonSink>>(md, "ShardedNdjsonSink")
    .def(py::init<std::string, std::size_t>(), py::arg("out_dir"), py::arg("shard_size"))
    .def(py::init<std::string, std::size_t, std::size_t>(), py::arg("out_dir"), py::arg("shard_size"), py::arg("max_shard_bytes"))
    .def("files", &ShardedNdjsonSink::files);

  py::class_<LmdbSink, IDynamicSink, std::shared_ptr<LmdbSink>>(md, "LmdbSink")
    .def(py::init<std::shared_ptr<lahuta::LMDBDatabase>, std::size_t>(),
         py::arg("db"), py::arg("batch_size") = 1024u);

  py::class_<SlowSink, IDynamicSink, std::shared_ptr<SlowSink>>(md, "SlowSink")
    .def(py::init<double, double>(),
         py::arg("sleep_seconds"), py::arg("flush_sleep_seconds") = 0.0,
         R"doc(Test sink that sleeps during write()/flush. Used for flush-timeout testing.)doc")
    .def("count", &SlowSink::count, R"doc(Return number of writes completed.)doc");
}

} // namespace lahuta::bindings

#endif // LAHUTA_BINDINGS_DYNAMIC_SINKS_HPP
