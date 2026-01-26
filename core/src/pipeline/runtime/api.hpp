#ifndef LAHUTA_PIPELINE_RUNTIME_API_HPP
#define LAHUTA_PIPELINE_RUNTIME_API_HPP

//
// Pipeline Runtime API.
//
// Usage:
//   #include "pipeline/runtime/api.hpp"
//   #include "sinks/ndjson.hpp"  // pick your sink
//
//   using namespace lahuta::pipeline;
//
//   StageManager mgr(from_lmdb(...));
//   mgr.add_task(...);
//   mgr.connect_sink("out", std::make_shared<NdjsonFileSink>(...));
//   mgr.run(threads);
//

// IWYU pragma: begin_exports
#include "pipeline/runtime/manager.hpp"
#include "pipeline/ingest/factory.hpp"
#include "pipeline/io/backpressure.hpp"
#include "pipeline/metrics/progress_observer.hpp"
// IWYU pragma: end_exports

#endif // LAHUTA_PIPELINE_RUNTIME_API_HPP