#ifndef LAHUTA_SOURCES_VECTOR_SOURCE_HPP
#define LAHUTA_SOURCES_VECTOR_SOURCE_HPP

#include <string>
#include <vector>

namespace lahuta::sources {

// old code
struct VectorSource {
  using value_type = std::string;

  VectorSource(std::vector<std::string> items) : items_(std::move(items)), idx_(0) {}

  std::optional<value_type> next() {
    if (idx_ < items_.size()) return items_[idx_++];
    return std::nullopt;
  }

private:
  std::vector<std::string> items_;
  size_t idx_;
};

} // namespace lahuta::sources

#endif // LAHUTA_SOURCES_VECTOR_SOURCE_HPP
