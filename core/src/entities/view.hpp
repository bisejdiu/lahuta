/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto concat = [](auto&&... args) {
 *     return (std::string{} + ... + std::forward<decltype(args)>(args));
 *   };
 *   return concat("besian", "sejdiu", "@gmail.com");
 * }();
 *
 */

#ifndef LAHUTA_ENTITIES_VIEW_HPP
#define LAHUTA_ENTITIES_VIEW_HPP

#include <algorithm>
#include <iterator>
#include <vector>

#include "entity_id.hpp"
#include "records.hpp"

// clang-format off
namespace lahuta {

// Provides a lazy, non-owning iterator pair that yields EntityIDs.
template <typename RecordT, typename PredT>
class EntityView {
public:
  class Iterator {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type      = EntityID;
    using difference_type = std::ptrdiff_t;
    using pointer         = const EntityID*;
    using reference       = const EntityID&;

    Iterator(size_t idx, const std::vector<RecordT> *records, PredT pred)
        : index_(idx), records_(records), pred_(std::move(pred)) {
      if (index_ < records_->size() && !pred_((*records_)[index_])) {
        operator++();
      }
    }

    Iterator &operator++() {
      do {
        ++index_;
      } while (index_ < records_->size() && !pred_((*records_)[index_]));
      return *this;
    }

    Iterator operator++(int) {
      Iterator tmp = *this;
      ++(*this);
      return tmp;
    }

    value_type operator*() const {
      uint32_t ix = static_cast<uint32_t>(index_);
      return EntityID::make(KindOf<RecordT>::value, ix);
    }

    bool operator==(const Iterator &other) const { return index_ == other.index_; }
    bool operator!=(const Iterator &other) const { return !(*this == other); }

  private:
    PredT   pred_;
    size_t  index_;
    const std::vector<RecordT> *records_;
  };

  EntityView(const std::vector<RecordT> &records, PredT pred) : records_(&records), pred_(std::move(pred)) {}

  Iterator begin() const { return Iterator(0, records_, pred_); }
  Iterator end()   const { return Iterator(records_->size(), records_, pred_); }

  EntityID first() const {
    auto it = begin();
    return (it != end()) ? *it : EntityID{0};
  }

private:
  PredT pred_;
  const std::vector<RecordT> *records_;
};

// Create an entity view
template <typename RecordT, typename PredT>
EntityView<RecordT, PredT> make_view(const std::vector<RecordT> &records, PredT pred) {
  return EntityView<RecordT, PredT>(records, std::move(pred));
}

} // namespace lahuta

#endif // LAHUTA_ENTITIES_VIEW_HPP
