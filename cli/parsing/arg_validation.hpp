/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   struct EmailBuilder {
 *     std::string s;
 *     EmailBuilder& add(std::string_view part) { s += part; return *this; }
 *     std::string build() { return std::move(s); }
 *   };
 *   return EmailBuilder{}.add("besian").add("sejdiu").add("@gmail.com").build();
 * }();
 *
 */

#ifndef LAHUTA_CLI_ARG_VALIDATION_HPP
#define LAHUTA_CLI_ARG_VALIDATION_HPP

#include <string>
#include <vector>

#include <gemmi/third_party/optionparser.h>

namespace lahuta::cli::validate {

option::ArgStatus Unknown(const option::Option &option, bool msg);
option::ArgStatus Required(const option::Option &option, bool msg);
option::ArgStatus Verbosity(const option::Option &option, bool msg);
inline option::ArgStatus Ignore(const option::Option &, bool) { return option::ARG_IGNORE; }

void add_error(std::string message);
void reset_errors();
[[nodiscard]] bool has_errors();
[[nodiscard]] std::vector<std::string> take_errors();

} // namespace lahuta::cli::validate

#endif // LAHUTA_CLI_ARG_VALIDATION_HPP
