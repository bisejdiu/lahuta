#ifndef LAHUTA_FACTORIZE_HPP
#define LAHUTA_FACTORIZE_HPP

#include <cstddef>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace lahuta::common {

inline std::vector<int> factorize(const std::vector<std::string> &labels) {
  std::vector<int> ids(labels.size());

  // hash map from labels to ids
  std::unordered_map<std::string_view, int> label_to_id;
  label_to_id.reserve(labels.size());

  int current_id = 0;
  for (size_t i = 0; i < labels.size(); ++i) {
    std::string_view label = labels[i];
    auto it = label_to_id.find(label);
    if (it == label_to_id.end()) {
      label_to_id[label] = current_id;
      ids[i] = current_id;
      ++current_id;
    } else {
      ids[i] = it->second;
    }
  }

  return ids;
}

} // namespace lahuta::common

#endif // LAHUTA_FACTORIZE_HPP