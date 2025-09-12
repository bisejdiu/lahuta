#ifndef LAHUTA_SPAN_HPP
#define LAHUTA_SPAN_HPP

#include <cstddef>

namespace lahuta {

// A very simple, leightweight, span class. All members are constexpr. Class can be templated
template <typename T>
class span {
  const T    *data_;
  std::size_t size_;

public:
  constexpr span(const T *data, std::size_t size) noexcept
    : data_(data), size_(size) {}

  template <typename C>
  constexpr span(const C &data_container) noexcept
    : data_(data_container.data()), size_(data_container.size()) {}

  constexpr const T &operator[](std::size_t idx) const noexcept { return data_[idx]; }
  constexpr std::size_t size() const noexcept { return size_; }
  constexpr const T    *data() const noexcept { return data_; }

  constexpr const T *begin() const noexcept { return data_; }
  constexpr const T *end()   const noexcept { return data_ + size_; }
};

} // namespace lahuta

#endif // LAHUTA_SPAN_HPP
