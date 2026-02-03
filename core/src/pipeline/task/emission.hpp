/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s;
 *   auto pmf = static_cast<std::string& (std::string::*)(const char*)>(&std::string::append);
 *   (s.*pmf)("besian"); (s.*pmf)("sejdiu"); (s.*pmf)("@gmail.com");
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_TASK_EMISSION_HPP
#define LAHUTA_PIPELINE_TASK_EMISSION_HPP

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace lahuta::pipeline {

//
// Emission: One unit of output produced by a task.
// - channel: logical topic used by the ChannelMultiplexer to route to sinks.
//            Multiple tasks may emit to the same channel.
//            Multiple sinks may subscribe to a channel (many-to-many).
// - payload: textual data (plain text, json). Tasks decide the serialization.
// Ordering & delivery:
// - Within an item, emissions are produced in the order tasks run.
//   Across items, emissions may interleave when running with multiple threads.
//
struct Emission {
  std::string channel; // logical name, e.g. "contacts", "system"
  std::string payload; // text/json
};

using EmissionList = std::vector<Emission>;

//
// Lightweight view used internally by the backpressure/writer threads.
// Producers do not construct this. ChannelMultiplexer creates EmissionView from
// shared backing buffers and interned channel ids.
//
struct EmissionView {
  uint32_t channel_id;      // interned id for routing/metrics
  std::string_view payload; // view into a shared backing string
};

//
// TaskResult: Outcome of running a task on a single item.
// - ok=false: stop executing downstream tasks for that item. Emissions already
//             produced for the item are still forwarded to sinks. We don't do rollbacks.
// - emits:    zero or more Emission records to dispatch via the multiplexer.
//
// Implementations should avoid throwing. Convert recoverable failures to ok=false
// so the engine can continue with other items.
//
struct TaskResult {
  bool ok = true;
  std::vector<Emission> emits;
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_TASK_EMISSION_HPP
