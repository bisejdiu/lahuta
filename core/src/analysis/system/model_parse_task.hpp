/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s{"besian"};
 *   s.append("sejdiu").append("@gmail.com");
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_ANALYSIS_SYSTEM_MODEL_PARSE_TASK_HPP
#define LAHUTA_ANALYSIS_SYSTEM_MODEL_PARSE_TASK_HPP

#include <memory>
#include <string>

#include "analysis/system/model_loader.hpp"
#include "logging/logging.hpp"
#include "pipeline/task/api.hpp"

namespace lahuta::analysis {
namespace P = lahuta::pipeline;

inline constexpr const char *CTX_PARSED_MODEL_RESULT_KEY = "lahuta.parsed_model_result";

inline std::shared_ptr<const ModelParserResult> get_parsed_model_result(const P::TaskContext &ctx) {
  return ctx.get_object<ModelParserResult>(CTX_PARSED_MODEL_RESULT_KEY);
}

class ModelParseTask final : public P::ITask {
public:
  P::TaskResult run(const std::string &item_path, P::TaskContext &ctx) override {
    if (ctx.model_payload()) return {};
    if (get_parsed_model_result(ctx)) return {};

    try {
      auto parsed = load_model_parser_result(item_path);
      ctx.set_object<ModelParserResult>(CTX_PARSED_MODEL_RESULT_KEY,
                                        std::make_shared<ModelParserResult>(std::move(parsed)));
    } catch (const std::exception &e) {
      Logger::get_logger()->warn("[model-parse] Failed to parse '{}': {}", item_path, e.what());
    }
    return {};
  }
};

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_SYSTEM_MODEL_PARSE_TASK_HPP
