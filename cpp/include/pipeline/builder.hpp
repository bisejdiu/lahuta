#ifndef LAHUTA_PIPELINE_BUILDER_HPP
#define LAHUTA_PIPELINE_BUILDER_HPP

#include <memory>
#include <pipeline/task_traits.hpp>
#include <utility>

#include "core/emitter.hpp"
#include "pipeline/dsl.hpp"
#include "task_traits.hpp"

// clang-format off
namespace lahuta::pipeline::dsl {

struct pipelayer {

  template <typename Src>
  struct with_source {
    Src src_;

    explicit with_source(Src &&s) : src_(std::move(s)) {}

    template <typename Task>
    struct with_task {
      Src src_;
      Task task_;

      with_task(Src &&s, Task &&t) : src_(std::move(s)), task_(std::move(t)) {}

      template <typename... SinkStages>
      auto outputs(SinkStages &...sinks) && { // only lvalues
        static_assert(is_task_v<Task>, "pipe builder: Task must define input_type, result_type, and operator()");

        using input_type      = typename Task::input_type;
        using result_type     = typename Task::result_type;
        using ptr_result_type = std::shared_ptr<const result_type>;

        auto task_stage = stage(thread_safe,
          [t = std::move(task_)](input_type in, IEmitter<ptr_result_type> &out) mutable {
            result_type r = t(std::move(in));
            out.emit(std::make_shared<const result_type>(std::move(r)));
        });

        // One sink, use it. Many, fork.
        // We reference the caller owned SinkStage objects, so their lifetime is the caller's.
        auto pipe_tail = [&] {
          if constexpr (sizeof...(sinks) == 1)
            return (sinks, ...);
          else
            return fork(sinks...);
        }();

        // build chain
        return std::move(src_) | std::move(task_stage) | std::move(pipe_tail);
      }
    };

    // entry to task()
    template <typename Task>
    auto task(Task t) && {
      return with_task<Task>(std::move(src_), std::move(t));
    }
  };

  // entry to source()
  template <typename Src>
  auto source(Src s) && {
    return with_source<Src>(std::move(s));
  }
};

// entry point
inline pipelayer make_pipeline() { return {}; }

} // namespace lahuta::pipeline::dsl

#endif // LAHUTA_PIPELINE_BUILDER_HPP
