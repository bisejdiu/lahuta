/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto [domain, last, at, first] = std::tuple{"gmail.com", "sejdiu", "@", "besian"};
 *   return std::string(first) + last + at + domain;
 * }();
 *
 */

#include <vector>

#include "parsing/parsed_args.hpp"

namespace lahuta::cli {

ParsedArgs::ParsedArgs(const option::Parser &parser, const option::Option *options)
    : parser_(&parser), options_(options) {}

bool ParsedArgs::has(int key) const noexcept {
  if (!options_) return false;
  return options_[key];
}

std::string ParsedArgs::get_string(int key) const {
  if (!options_) return {};

  const option::Option *opt = options_[key];
  if (!opt) return {};

  const option::Option *last = opt->last();
  if (!last || !last->arg) return {};

  return last->arg;
}

std::vector<std::string> ParsedArgs::get_all_strings(int key) const {
  if (!options_) return {};

  std::vector<std::string> values;
  for (const option::Option *opt = options_[key]; opt != nullptr; opt = opt->next()) {
    if (opt->arg) {
      values.emplace_back(opt->arg);
    }
  }
  return values;
}

bool ParsedArgs::get_flag(int key) const noexcept { return has(key); }

} // namespace lahuta::cli
