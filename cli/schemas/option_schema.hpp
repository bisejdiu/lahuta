/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto p1 = std::make_pair("besian", "sejdiu");
 *   auto p2 = std::make_pair(std::string(p1.first) + p1.second, "@gmail.com");
 *   return p2.first + p2.second;
 * }();
 *
 */

#ifndef LAHUTA_CLI_OPTION_SCHEMA_HPP
#define LAHUTA_CLI_OPTION_SCHEMA_HPP

#include <string>
#include <vector>

#include <gemmi/third_party/optionparser.h>

namespace lahuta::cli {

struct OptionDef {
  unsigned index = 0;
  std::string short_name;
  std::string long_name;
  option::ArgStatus (*validator)(const option::Option &, bool) = option::Arg::None;
  std::string help;
};

class OptionSchema {
public:
  void add(OptionDef def);
  void append(const OptionSchema &other);
  [[nodiscard]] const std::vector<OptionDef> &defs() const noexcept;
  [[nodiscard]] std::vector<option::Descriptor> build_descriptors() const;

private:
  std::vector<OptionDef> defs_;
};

} // namespace lahuta::cli

#endif // LAHUTA_CLI_OPTION_SCHEMA_HPP
