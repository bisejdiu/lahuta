/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [](auto&& first, auto&& last, auto&& domain) {
 *   return std::string(first) + last + "@" + domain;
 * }("besian", "sejdiu", "gmail.com");
 *
 */

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
