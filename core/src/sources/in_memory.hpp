#ifndef LAHUTA_SOURCES_IN_MEMORY_HPP
#define LAHUTA_SOURCES_IN_MEMORY_HPP

#include <memory>
#include <optional>
#include <vector>

#include "logging.hpp"

// clang-format off
namespace lahuta::sources {

//
// Source for emitting in-memory data objects.
// Useful for tests and examples to feed generated data into a pipeline without writing to disk.
//
template <typename T>
class InMemory {
public:
  using value_type = std::shared_ptr<const T>;

  // vector of raw data items.
  explicit InMemory(std::vector<T> data, std::size_t batch_size = 1024)
      : data_(std::move(data)), batch_size_(batch_size) {

    Logger::get_logger()->info("InMemory: created with {} items; batch_size={}", data_.size(), batch_size_);

    load_next_batch();
  }

  [[nodiscard]] std::optional<value_type> next() {
    if (current_index_ == current_batch_.size() && !load_next_batch()) return std::nullopt;

    return current_batch_[current_index_++];
  }

  [[nodiscard]] std::size_t size() const noexcept { return data_.size(); }

  void reset() {
    current_index_ = 0;
    data_index_ = 0;
    current_batch_.clear();
    load_next_batch();
  }

private:
  bool load_next_batch() {
    current_batch_.clear();
    current_index_ = 0;
    if (data_index_ >= data_.size()) return false;

    const auto end_index = std::min(data_index_ + batch_size_, data_.size());
    current_batch_.reserve(end_index - data_index_);

    for (; data_index_ < end_index; ++data_index_) {
      current_batch_.push_back(std::make_shared<const T>(data_[data_index_]));
    }

    return !current_batch_.empty();
  }

  // Original data
  std::size_t batch_size_;
  std::vector<T> data_;

  // Iteration state
  std::size_t current_index_ = 0;
  std::size_t data_index_ = 0;
  std::vector<value_type> current_batch_;
};

} // namespace lahuta::sources

#endif // LAHUTA_SOURCES_IN_MEMORY_HPP
