#ifndef LAHUTA_PIPELINE_SINK_ADAPTER_HPP
#define LAHUTA_PIPELINE_SINK_ADAPTER_HPP

#include "serialization/serializer.hpp"
#include <type_traits>
#include <utility>

// clang-format off
namespace lahuta::pipeline {

// adapter so the pipeline fanout will treat as IEmitter<Ptr<Result>>
template <
    typename FormatTag,
    typename Payload,
    template <typename> typename Backend,
    typename... BackendCtorArgs>
class SinkAdapter {
  static_assert(
      std::is_pointer_v<typename Payload::element_type> == false,
      "Pass shared_ptr<const Result>, not raw pointer");

  using Result      = typename Payload::element_type;
  using SerializerT = serialization::Serializer<FormatTag, Result>;

public:
  using input_type = Payload;
  template <class... A>
  explicit SinkAdapter(A &&...a)
    : backend_(std::forward<A>(a)...) {}

  void emit(Payload &&p) { backend_.write(SerializerT::serialize(*p)); }

private:
  Backend<SerializerT> backend_;
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_SINK_ADAPTER_HPP
