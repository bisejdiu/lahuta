#ifndef LAHUTA_PIPELINE_PIPELINE_ITEM_HPP
#define LAHUTA_PIPELINE_PIPELINE_ITEM_HPP

#include <cstdint>
#include <memory>
#include <optional>
#include <string>

#include "pipeline/frame.hpp"

namespace lahuta {

struct StreamSession;

struct PipelineItem {
  std::string session_id;
  std::string item_path;
  std::uint64_t conformer_id = 0;
  std::optional<double> timestamp_ps;

  std::shared_ptr<const StreamSession> session;
  std::shared_ptr<FrameHandle> frame;
};

} // namespace lahuta

#endif // LAHUTA_PIPELINE_PIPELINE_ITEM_HPP
