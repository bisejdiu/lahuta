#ifndef LAHUTA_MAPPING_DISJOINT_SET_HPP
#define LAHUTA_MAPPING_DISJOINT_SET_HPP

#include "_defs.hpp"
#include <mutex>
#include <vector>

// clang-format off
namespace lahuta::mapping {

class DynamicEquivalenceForest { // DSU
private:
  mutable std::mutex m_;
  mutable std::vector<u32> parent_; // mutable for path compression
  std::vector<u32> rank_;           // union‑by‑rank
  std::vector<u32> head_;           // head of singly‑linked member list
  std::vector<u32> next_;           // next pointer for each element
  std::vector<u32> size_;           // size of each root

  u32 _find_nl(u32 x) const {
    if (parent_[x] != x) {
      parent_[x] = _find_nl(parent_[x]);
    }
    return parent_[x];
  }

public:
  u32 make() {
    std::lock_guard<std::mutex> lock(m_);
    u32 id = static_cast<u32>(parent_.size());

    parent_.push_back(id);
    rank_  .push_back(0);
    head_  .push_back(id);
    next_  .push_back(static_cast<u32>(-1));
    size_  .push_back(1);
    return id;
  }

  u32 find(u32 x) const {
    std::lock_guard<std::mutex> lock(m_);
    return _find_nl(x);
  }

  void join(u32 x, u32 y) {
    std::lock_guard<std::mutex> lock(m_);
    u32 rx = _find_nl(x);
    u32 ry = _find_nl(y);
    if (rx == ry) return;

    // keep rx the new root (rank / size heuristic)
    if (rank_[rx] < rank_[ry] || (rank_[rx] == rank_[ry] && size_[rx] < size_[ry])) {
      std::swap(rx, ry);
    }

    // attach ry under rx
    parent_[ry] = rx;
    if (rank_[rx] == rank_[ry]) {
      ++rank_[rx];
    }

    // splice member lists: [head_rx] + [head_ry]
    u32 tail = head_[ry];
    while (next_[tail] != static_cast<u32>(-1)) {
      tail = next_[tail];
    }

    next_[tail] = head_[rx];
    head_[rx]   = head_[ry];
    size_[rx]  += size_[ry];
  }

  /// Get the root of the set containing x
  template <typename T = u32>
  std::vector<T> get_members(u32 x) const {
    std::lock_guard<std::mutex> lock(m_);
    u32 root = _find_nl(x);

    std::vector<T> out;
    out.reserve(size_[root]);

    for (u32 v = head_[root]; v != static_cast<u32>(-1); v = next_[v])
      out.push_back(static_cast<T>(v));

    return out; // this should be O(|S|)
  }

  /// Check if two elements are in the same set
  bool is_same(u32 x, u32 y) const {
    std::lock_guard<std::mutex> lock(m_);
    return _find_nl(x) == _find_nl(y);
  }

  void clear() {
    std::lock_guard<std::mutex> lock(m_);
    parent_.clear();
    rank_  .clear();
    head_  .clear();
    next_  .clear();
    size_  .clear();
  }

  size_t size() const {
    std::lock_guard<std::mutex> lock(m_);
    return parent_.size();
  }
};

} // namespace lahuta::mapping

#endif // LAHUTA_MAPPING_DISJOINT_SET_HPP
