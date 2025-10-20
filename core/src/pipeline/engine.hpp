#ifndef LAHUTA_PIPELINE_ENGINE_HPP
#define LAHUTA_PIPELINE_ENGINE_HPP

#include <ctpl.h>

#include "core/emitter.hpp"

// clang-format off
namespace lahuta::pipeline {

class PipelineEngine {
public:
  explicit PipelineEngine(size_t threads=1) : pool_(threads), in_flight_(0) {}

  template<typename Source, typename Stage, typename Em>
  void run(Source& src, Stage& st, Em& em){
    static_assert(std::is_base_of_v<IEmitter<typename Stage::output_type>, Em>);

    if (!st.thread_safe() || pool_.size()==1)
      serial(src, st, em);
    else
      parallel(src, st, em);
  }

private:
  template<typename Src, typename St, typename Em>
  static void serial(Src& s, St& st, Em& em){
    while(auto v = s.next()) st.process(std::move(*v), em);
  }

  template<typename Src, typename St, typename Em>
  void parallel(Src& s, St& st, Em& em){
    while(auto v = s.next()){
      in_flight_.fetch_add(1, std::memory_order_relaxed);
        pool_.push([this, &st, &em, val = std::move(*v)](int) mutable {
        try {
          st.process(std::move(val), em);
        } catch(...) {
          //
          // TODO: How do we handle exceptions?
          //
        }
        if (in_flight_.fetch_sub(1) == 1) {
          std::scoped_lock lk(m_);
          cv_.notify_one();
        }
      });
    }
    std::unique_lock lk(m_);
    cv_.wait(lk, [&]{ return in_flight_.load() == 0; });
  }

private:
  std::atomic_size_t      in_flight_;
  std::mutex              m_;
  std::condition_variable cv_;
  ctpl::thread_pool       pool_;
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_ENGINE_HPP
