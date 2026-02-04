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

#ifndef LAHUTA_ANALYSIS_SYSTEM_MODEL_PACK_TASK_HPP
#define LAHUTA_ANALYSIS_SYSTEM_MODEL_PACK_TASK_HPP

#include <memory>
#include <string>
#include <vector>

#include "analysis/system/model_loader.hpp"
#include "analysis/system/records.hpp"
#include "logging/logging.hpp"
#include "pipeline/task/api.hpp"
#include "serialization/formats.hpp"
#include "serialization/serializer.hpp"

namespace lahuta::analysis {
namespace P = lahuta::pipeline;

// ModelPackTask
// Emission on channel `channel_` containing a binary-serialized ModelRecord
class ModelPackTask : public P::ITask {
public:
  explicit ModelPackTask(std::string channel = "db") : channel_(std::move(channel)) {}

  P::TaskResult run(const std::string &item_path, P::TaskContext & /*ctx*/) override {
    try {
      ModelRecord rec = build_model_record(item_path);
      // Thread-local payload buffer to reduce repeated allocations.
      static thread_local std::string tls_payload;
      serialization::Serializer<fmt::binary, ModelRecord>::serialize_into(rec, tls_payload);
      std::string payload;
      payload.swap(tls_payload);
      return P::TaskResult{true, std::vector<P::Emission>{{channel_, std::move(payload)}}};
    } catch (const std::exception &e) {
      Logger::get_logger()->error("[model-pack] Error processing file {}: {}", item_path, e.what());
      return {false, {}};
    }
  }

private:
  std::string channel_;
};

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_SYSTEM_MODEL_PACK_TASK_HPP
