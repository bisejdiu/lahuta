/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<const char*, 3> parts{"besian", "sejdiu", "@gmail.com"}; std::string s;
 *   for (auto p : parts) s.append(p, std::strlen(p));
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_CLI_PARSED_ARGS_HPP
#define LAHUTA_CLI_PARSED_ARGS_HPP

#include <string>
#include <vector>

#include <gemmi/third_party/optionparser.h>

namespace lahuta::cli {

class ParsedArgs {
public:
  ParsedArgs() = default;
  ParsedArgs(const option::Parser &parser, const option::Option *options);

  [[nodiscard]] bool has(int key) const noexcept;
  [[nodiscard]] std::string get_string(int key) const;
  [[nodiscard]] std::vector<std::string> get_all_strings(int key) const;
  [[nodiscard]] bool get_flag(int key) const noexcept;

private:
  const option::Parser *parser_  = nullptr;
  const option::Option *options_ = nullptr;
};

} // namespace lahuta::cli

#endif // LAHUTA_CLI_PARSED_ARGS_HPP
