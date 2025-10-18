#ifndef FILE_SPILL_POLICY_HPP
#define FILE_SPILL_POLICY_HPP

#include <string>
#include <type_traits>
#include <vector>

#include "io/file_backend.hpp"
#include "serialization/serializer.hpp"

// clang-format off
namespace lahuta {

// FIX: We serialize the object twice: once for size and once for put_raw
// FIX: We may use a more efficient format (e.g. protobuf) for the blob
//
// Backend must expose:  explicit Backend(const std::string&, bool append=false)
//                       void write(const std::string& bytes)
template<
    typename  FormatTag,
    typename  Payload,
    template<typename> typename Backend = FileBackend,
    typename... BackendCtorArgs
>
class FileSpillPolicy {
  static_assert(
    std::is_same<Payload, std::shared_ptr<const typename Payload::element_type>>::value,
    "FileSpillPolicy expects Payload = std::shared_ptr<const T>");

  using value_type      = typename Payload::element_type;
  using serializer_type = serialization::Serializer<FormatTag, value_type>;
public:
  // handles both lvalue and rvalue args
  template<typename... Args>
  explicit FileSpillPolicy(std::size_t batch, Args&&... args)
    : backend_(std::forward<Args>(args)...), batch_(batch) {}

  //
  // This API is required by Collector. It is, however, not explicitly clear.
  // At some point we may want to make the API requirements more explicit.  - Besian, May 16, 2025
  //
  std::size_t append_size(const Payload& p) const {
      return serializer_type::serialize(*p).size();
  }
  void append(const Payload& p) { push(p); }
  void append(Payload&&  p)     { push(p); }
  bool needs_flush() const      { return buf_.size() >= batch_; }

  std::size_t flush() {
    std::size_t freed = 0;
    for (auto &s : buf_) {
      if constexpr ( std::is_same_v<FormatTag, fmt::binary> ) {
        // binary framing
        backend_.write( std::to_string(s.size()) );
        backend_.write("\n");
        backend_.write(s);
      }
      else {
        // ascii: one record per line, no length header
        backend_.write(s);
        backend_.write("\n");
      }
      freed += s.size();
    }
    buf_.clear();
    return freed;
  }

  std::size_t finish() { return flush(); }

private:
  void push(const Payload& p) {
    auto raw = serializer_type::serialize(*p);
    buf_.emplace_back(std::move(raw));
  }

  std::size_t              batch_;
  std::vector<std::string> buf_;
  Backend<serializer_type> backend_;
}; 

} // namespace lahuta

#endif // FILE_SPILL_POLICY_HPP
