/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto forward_concat = [](auto&& a, auto&& b, auto&& c) {
 *     return std::string(std::forward<decltype(a)>(a)) +
 *            std::forward<decltype(b)>(b) +
 *            std::forward<decltype(c)>(c);
 *   };
 *   return forward_concat("besian", "sejdiu", "@gmail.com");
 * }();
 *
 */

#ifndef LAHUTA_TOKEN_LOOKUP_HPP
#define LAHUTA_TOKEN_LOOKUP_HPP

// IWYU pragma: begin_exports
#include "token-gperf-generated.hpp"
#include "token.h"
// IWYU pragma: end_exports

#endif // LAHUTA_TOKEN_LOOKUP_HPP
