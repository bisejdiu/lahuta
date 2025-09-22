#ifndef LAHUTA_PIPELINE_RUN_HPP
#define LAHUTA_PIPELINE_RUN_HPP

#include "pipe.hpp"
#include "pipeline/engine.hpp"
#include "sources.hpp"
#include <type_traits>

// clang-format off
namespace lahuta::pipeline::dsl {

// gather stages in order (tuple of references)
template<typename P> struct stages_tuple;

template <typename P> struct stages_tuple;                                             // primary
template <typename P> struct stages_tuple<P &>     : stages_tuple<std::decay_t<P>> {}; // strip ref
template <typename P> struct stages_tuple<const P> : stages_tuple<std::decay_t<P>> {}; // strip cv

// leaf stage
template<typename S>
struct stages_tuple {
  using type = std::tuple<S&>;
  static type get(S& s) { return {s}; }
};

// source_t (no stages)
template<typename Src>
struct stages_tuple< source_t<Src> > {
  using type = std::tuple<>;
  static type get(source_t<Src>&) { return {}; }
};

// pipe_t : concat LHS and RHS
template<typename L, typename R>
struct stages_tuple< pipe_t<L,R> > {
  using type = decltype( std::tuple_cat(
    stages_tuple<L>::get( std::declval<L&>() ),
    stages_tuple<R>::get( std::declval<R&>() )
  ));
  static type get(pipe_t<L,R>& p) {
    return std::tuple_cat(
      stages_tuple<L>::get(p.lhs),
      stages_tuple<R>::get(p.rhs)
    );
  }
};

// unwrap std::reference_wrapper if we are passed a collect(emitter)
template<typename T> auto& unwrap_sink(T& s)                        { return s;       }
template<typename U> auto& unwrap_sink(std::reference_wrapper<U> w) { return w.get(); }

template<typename... Ss>
class Chain {
public:
  using input_type  = typename std::tuple_element_t<0, std::tuple<Ss...>>::input_type;
  using output_type = typename std::tuple_element_t<sizeof...(Ss)-1,std::tuple<Ss...>>::output_type;

  constexpr bool thread_safe() const noexcept { return true; }

  explicit Chain(std::tuple<Ss&...> t) : st_(t) {}
  void process(input_type ptr, IEmitter<output_type>& out) { step<0>(std::move(ptr), out); }

private:
  template<std::size_t I, typename P>
  void step(P&& p, IEmitter<output_type>& out) {

    if constexpr (I+1 == sizeof...(Ss)) {
      std::get<I>(st_).process(std::forward<P>(p), out);

    } else {

      using Mid = typename std::tuple_element_t<I, std::tuple<Ss...>>::output_type;
      struct Forward : IEmitter<Mid> {
        Chain& ch; IEmitter<output_type>& out;
        Forward(Chain& c, IEmitter<output_type>& o) : ch(c), out(o) {}
        void emit(Mid q) override { ch.template step<I+1>(std::move(q), out); }
      } fwd{*this,out};

      std::get<I>(st_).process(std::forward<P>(p), fwd);

    }
  }

private:
  std::tuple<Ss&...> st_;
};

/// Run the pipeline
template<typename Pipeline, typename FinalSink>
void run(pipe_t<Pipeline, FinalSink>& pl, std::size_t threads = 1) {
  auto& src    = find_source(pl.lhs);
  auto  stages = stages_tuple<Pipeline>::get(pl.lhs);
  Chain chain{stages};
  PipelineEngine{threads}.run(src, chain, unwrap_sink(pl.rhs));
}

} // namespace lahuta::pipeline::dsl

#endif // LAHUTA_PIPELINE_RUN_HPP
