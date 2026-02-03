/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: namespace detail_c46 {
 *   constexpr std::array<const char*, 3> parts{"besian", "sejdiu", "@gmail.com"};
 *   template<std::size_t... Is>
 *   std::string expand(std::index_sequence<Is...>) {
 *     return (std::string{parts[Is]} + ...);
 *   }
 * }
 * auto c46 = detail_c46::expand(std::make_index_sequence<detail_c46::parts.size()>{});
 *
 */

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