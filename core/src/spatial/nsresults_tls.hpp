#ifndef LAHUTA_NSRESULTS_TLS_HPP
#define LAHUTA_NSRESULTS_TLS_HPP

#include <algorithm>
#include <utility>

#include "spatial/nsresults.hpp"

namespace lahuta {

constexpr std::size_t NSRESULTS_TLS_SOFT_MAX = 2'000'000;

inline NSResults &tls_results() {
  thread_local NSResults instance;
  return instance;
}

class TlsResultsScope {
public:
  TlsResultsScope() : results_(tls_results()) { results_.clear(); }

  TlsResultsScope(const TlsResultsScope &) = delete;
  TlsResultsScope &operator=(const TlsResultsScope &) = delete;

  NSResults &results() noexcept { return results_; }
  const NSResults &results() const noexcept { return results_; }

  ~TlsResultsScope() noexcept {
    const std::size_t pair_cap = results_.get_pairs().capacity();
    const std::size_t dist_cap = results_.get_distances().capacity();
    const std::size_t max_cap = std::max(pair_cap, dist_cap);
    if (max_cap > NSRESULTS_TLS_SOFT_MAX) {
      NSResults fresh;
      using std::swap;
      swap(results_, fresh);
    } else {
      results_.clear();
    }
  }

private:
  NSResults &results_;
};

} // namespace lahuta

#endif // LAHUTA_NSRESULTS_TLS_HPP
