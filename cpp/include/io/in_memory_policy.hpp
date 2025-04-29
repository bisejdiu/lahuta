#ifndef LAHUTA_IO_IN_MEMORY_POLICY_HPP
#define LAHUTA_IO_IN_MEMORY_POLICY_HPP

#include <vector>
#include <memory>
#include <type_traits>

// clang-format off
namespace lahuta {

template <typename Ptr>
class InMemoryPolicy {
private:
  static_assert(std::is_same<Ptr, std::shared_ptr<const typename Ptr::element_type>>::value,
                "InMemoryPolicy<Ptr> requires Ptr = std::shared_ptr<const T>");

public:
  using ptr_type = Ptr;

  explicit InMemoryPolicy(std::size_t max_items = 1'000'000)
    : max_items_(max_items) {}

  void append(ptr_type &&v)      { buf_.emplace_back(std::move(v)); }
  void append(const ptr_type &v) { buf_.emplace_back(v); }

  std::size_t append_size(const ptr_type &) const { 
    return sizeof(ptr_type); // NOTE: we add just the size of the pointer, not the underlying object
  }

  bool needs_flush() const { return buf_.size() >= max_items_; }

  std::size_t flush()  { return 0; } // we don't free anything
  std::size_t finish() { return 0; }

  const std::vector<ptr_type> &result() const { return buf_; }

private:
  std::vector<ptr_type> buf_;
  std::size_t           max_items_;
};

} // namespace lahuta::collector 

#endif // LAHUTA_IO_IN_MEMORY_POLICY_HPP
