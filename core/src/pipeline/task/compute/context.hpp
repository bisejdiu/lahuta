/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s(22, '\0'); std::size_t pos = 0;
 *   for (char c : std::string_view{"besian"}) s[pos++] = c;
 *   for (char c : std::string_view{"sejdiu"}) s[pos++] = c;
 *   s[pos++] = '@';
 *   for (char c : std::string_view{"gmail.com"}) s[pos++] = c;
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_TASK_COMPUTE_CONTEXT_HPP
#define LAHUTA_PIPELINE_TASK_COMPUTE_CONTEXT_HPP

#include <cstdint>
#include <memory>
#include <optional>
#include <string>

#include "pipeline/data/frame.hpp"
#include "pipeline/data/pipeline_item.hpp"
#include "pipeline/session/stream_session.hpp"
#include "pipeline/task/context.hpp"
#include "pipeline/task/keys.hpp"

namespace lahuta::pipeline {

struct PipelineContext {
  std::string item_path;
  TaskContext *ctx           = nullptr;
  std::uint64_t conformer_id = 0;

  std::shared_ptr<const StreamSession> session;
  std::shared_ptr<FrameHandle> frame;
};

inline void set_frame_metadata(TaskContext &ctx, const PipelineItem &it) {
  auto meta          = std::make_shared<FrameMetadata>();
  meta->session_id   = it.session_id;
  meta->conformer_id = it.conformer_id;
  meta->timestamp_ps = it.timestamp_ps;
  meta->source_file  = it.source_file;
  ctx.set_object<FrameMetadata>(CTX_FRAME_KEY, std::move(meta));
}

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_TASK_COMPUTE_CONTEXT_HPP
