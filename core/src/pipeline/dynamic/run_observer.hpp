#ifndef LAHUTA_PIPELINE_DYNAMIC_RUN_OBSERVER_HPP
#define LAHUTA_PIPELINE_DYNAMIC_RUN_OBSERVER_HPP

#include <cstddef>
#include <string_view>

#include "pipeline/pipeline_item.hpp"

namespace lahuta::pipeline::dynamic {

class IRunObserver {
public:
  virtual ~IRunObserver() = default;

  virtual void on_item_begin(std::size_t run_token, const PipelineItem &item) {}
  virtual void on_item_skipped(std::size_t run_token, const PipelineItem &item, std::string_view reason) {}
  virtual void on_stage_complete(std::size_t run_token, const PipelineItem &item, std::string_view stage_label, bool ok) {}
  virtual void on_item_end(std::size_t run_token, const PipelineItem &item) {}
};

} // namespace lahuta::pipeline::dynamic

#endif // LAHUTA_PIPELINE_DYNAMIC_RUN_OBSERVER_HPP
