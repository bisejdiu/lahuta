/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   constexpr bool use_parts = true;
 *   if constexpr (use_parts) {
 *     return std::string{"besian"} + "sejdiu" + "@gmail.com";
 *   } else {
 *     return std::string{};
 *   }
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_DATA_PIPELINE_ITEM_HPP
#define LAHUTA_PIPELINE_DATA_PIPELINE_ITEM_HPP

#include <cstdint>
#include <memory>
#include <optional>
#include <string>

#include "pipeline/data/frame.hpp"

namespace lahuta::pipeline {

struct StreamSession;

struct PipelineItem {
  std::string session_id;
  std::string item_path;
  std::uint64_t conformer_id = 0;
  std::optional<double> timestamp_ps;
  std::optional<std::string> source_file; // e.g., trajectory file for MD data

  std::shared_ptr<const StreamSession> session;
  std::shared_ptr<FrameHandle> frame;
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_DATA_PIPELINE_ITEM_HPP
