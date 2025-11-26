#ifndef LAHUTA_SPAN_HPP
#define LAHUTA_SPAN_HPP

#include <cstddef>
#include <iterator>

// clang-format off
namespace lahuta {

template <class T>
class span {
  T*           data_ = nullptr;
  std::size_t  size_ = 0;

public:
  constexpr span() noexcept = default;
  constexpr span(T* p, std::size_t n) noexcept : data_(p), size_(n) {}

  template <class C> constexpr span(C& c) noexcept : data_(std::data(c)), size_(std::size(c)) {}
  template <class C> constexpr span(const C& c) noexcept : data_(std::data(c)), size_(std::size(c)) {}

  template <class C>
  span(C&&) = delete; // no temporaries

  constexpr T& operator[](std::size_t i) noexcept { return data_[i]; }
  constexpr const T& operator[](std::size_t i) const noexcept { return data_[i]; }

  [[nodiscard]] constexpr std::size_t size() const noexcept { return size_; }
  [[nodiscard]] constexpr bool empty() const noexcept { return size_ == 0; }
  [[nodiscard]] constexpr T* data() noexcept { return data_; }
  [[nodiscard]] constexpr const T* data() const noexcept { return data_; }
  [[nodiscard]] constexpr std::size_t size_bytes() const noexcept { return size_ * sizeof(T); }

  constexpr T* begin() noexcept { return data_; }
  constexpr T* end()   noexcept { return data_ + size_; }
  constexpr const T* begin() const noexcept { return data_; }
  constexpr const T* end()   const noexcept { return data_ + size_; }
};

} // namespace lahuta

#endif // LAHUTA_SPAN_HPP
