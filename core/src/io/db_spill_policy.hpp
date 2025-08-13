#ifndef LAHUTA_DB_SPILL_POLICY_HPP
#define LAHUTA_DB_SPILL_POLICY_HPP

#include "io/lmdb_backend.hpp"
#include "serialization/serializer.hpp"
#include <vector>

// clang-format off
namespace lahuta {

//
// We serialize the entire Task::result_type object, rather than the corresponding data member.
// This leads to redundancy in handling strings/keys, but it is more convenient to have the pairing
// of key and value in the same object. Below I have commented out the code that would
// serialize only the data member. The issue then is that we need an additional step in the pipeline
// to pair the key and value together.  - Besian, May 16, 2025
//
template<
    typename  FormatTag,
    typename  Payload,
    template<typename> typename Backend = LMDBBackend,
    typename... BackendCtorArgs
>
class DBSpillPolicy {
    static_assert(
      std::is_same_v<Payload, std::shared_ptr<const typename Payload::element_type>>,
      "DBSpillPolicy expects Payload = shared_ptr<const T>");

    using value_type      = typename Payload::element_type;
    using serializer_type = serialization::Serializer<FormatTag, value_type>;

    // using record_type     = typename Payload::element_type;
    // using data_type       = std::remove_cv_t<std::remove_reference_t<
    //                          decltype(std::declval<record_type>().data)
    //                        >>;
    // using serializer_type = serialization::Serializer<FormatTag, data_type>;

    Backend<serializer_type> backend_;
    const std::size_t    batch_;
    std::vector<std::pair<std::string,std::string>> buf_;

public:
    template<typename... Args>
    explicit DBSpillPolicy(std::size_t batch, Args&&... xs)
      : backend_(std::forward<Args>(xs)...), batch_(batch) {}

    bool needs_flush() const noexcept {
      return buf_.size() >= batch_;
    }

    std::size_t append_size(const Payload& p) const {
      // key + serialized-value
      auto const& rec = *p;
      auto const  val = serializer_type::serialize(rec);
      return rec.file_path.size() + val.size();
    }

    void append(const Payload& p) { push(p); }
    void append(Payload&&      p) { push(p); }

    std::size_t flush() {
      backend_.begin_txn();

      std::size_t freed = 0;
      for (auto& [k,v] : buf_) {
        backend_.write(k, v);
        freed += k.size() + v.size();
      }
      backend_.commit_txn();
      buf_.clear();
      return freed;
    }

    std::size_t finish() { return flush(); }

private:
    void push(const Payload& p) {
      auto const& rec = *p;
      auto val = serializer_type::serialize(rec);
      buf_.emplace_back(rec.file_path, std::move(val));
    }
};

} // namespace lahuta

#endif // LAHUTA_DB_SPILL_POLICY_HPP
