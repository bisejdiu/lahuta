/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto forward_concat = [](auto&& a, auto&& b, auto&& c) {
 *     return std::string(std::forward<decltype(a)>(a)) +
 *            std::forward<decltype(b)>(b) +
 *            std::forward<decltype(c)>(c);
 *   };
 *   return forward_concat("besian", "sejdiu", "@gmail.com");
 * }();
 *
 */

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
#include "sinks/lmdb.hpp"
#include "sinks/memory.hpp"
#include "sinks/ndjson.hpp"
#include "sinks/sharded_ndjson.hpp"

namespace py = pybind11;
namespace lahuta::bindings {
namespace P = lahuta::pipeline;

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
class SlowSink : public P::IDynamicSink {
public:
  explicit SlowSink(double sleep_seconds, double flush_sleep_seconds = 0.0)
      : write_sleep_(seconds_to_ms(sleep_seconds)), flush_sleep_(seconds_to_ms(flush_sleep_seconds)) {}

  void write(P::EmissionView e) override {
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

inline void bind_sinks(py::module_ &md) {
  py::class_<P::IDynamicSink, std::shared_ptr<P::IDynamicSink>> _(md, "Sink");

  py::class_<P::MemorySink, P::IDynamicSink, std::shared_ptr<P::MemorySink>>(md, "MemorySink")
      .def(py::init<>())
      .def("result", &P::MemorySink::result, R"doc(Return collected payloads.)doc")
      .def(
          "result_bytes",
          [](const P::MemorySink &s) {
            py::list out;
            for (const auto &buf : s.result_bytes()) {
              out.append(py::bytes(buf));
            }
            return out;
          },
          R"doc(Return collected payloads as bytes.)doc")
      .def("clear", &P::MemorySink::clear, R"doc(Clear collected payloads.)doc");

  py::class_<P::NdjsonFileSink, P::IDynamicSink, std::shared_ptr<P::NdjsonFileSink>>(md, "NdjsonSink")
      .def(py::init<std::string>(), py::arg("path"))
      .def("file", &P::NdjsonFileSink::file);

  py::class_<P::ShardedNdjsonSink, P::IDynamicSink, std::shared_ptr<P::ShardedNdjsonSink>>(
      md,
      "ShardedNdjsonSink")
      .def(py::init<std::string, std::size_t>(), py::arg("out_dir"), py::arg("shard_size"))
      .def(py::init<std::string, std::size_t, std::size_t>(),
           py::arg("out_dir"),
           py::arg("shard_size"),
           py::arg("max_shard_bytes"))
      .def("files", &P::ShardedNdjsonSink::files);

  py::class_<P::LmdbSink, P::IDynamicSink, std::shared_ptr<P::LmdbSink>>(md, "LmdbSink")
      .def(py::init<std::shared_ptr<lahuta::LMDBDatabase>, std::size_t>(),
           py::arg("db"),
           py::arg("batch_size") = 1024u);

  py::class_<SlowSink, P::IDynamicSink, std::shared_ptr<SlowSink>>(md, "SlowSink")
      .def(py::init<double, double>(),
           py::arg("sleep_seconds"),
           py::arg("flush_sleep_seconds") = 0.0,
           R"doc(Test sink that sleeps during write()/flush. Used for flush-timeout testing.)doc")
      .def("count", &SlowSink::count, R"doc(Return number of writes completed.)doc");
}

} // namespace lahuta::bindings

#endif // LAHUTA_BINDINGS_DYNAMIC_SINKS_HPP
