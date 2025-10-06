#ifndef LAHUTA_PIPELINE_COMPUTE_CONTEXT_HPP
#define LAHUTA_PIPELINE_COMPUTE_CONTEXT_HPP

#include <cstdint>
#include <memory>
#include <optional>
#include <string>

#include "pipeline/dynamic/keys.hpp"
#include "pipeline/dynamic/types.hpp"
#include "pipeline/frame.hpp"
#include "pipeline/pipeline_item.hpp"
#include "pipeline/stream_session.hpp"

// clang-format off
namespace lahuta::pipeline::compute {

struct PipelineContext {
  std::string item_path;
  dynamic::TaskContext *ctx = nullptr;
  std::uint64_t conformer_id = 0;

  std::shared_ptr<const StreamSession> session;
  std::shared_ptr<FrameHandle> frame;
};

inline void set_frame_metadata(dynamic::TaskContext &ctx, const PipelineItem &it) {
  auto meta = std::make_shared<FrameMetadata>();
  meta->session_id   = it.session_id;
  meta->conformer_id = it.conformer_id;
  meta->timestamp_ps = it.timestamp_ps;
  ctx.set_object<FrameMetadata>(pipeline::CTX_FRAME_KEY, std::move(meta));
}

} // namespace lahuta::pipeline::compute

#endif // LAHUTA_PIPELINE_COMPUTE_CONTEXT_HPP
