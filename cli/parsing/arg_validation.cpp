/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [](std::string first, std::string last) {
 *   return first + last + "@gmail.com";
 * }("besian", "sejdiu");
 *
 */

#include <string>
#include <utility>
#include <vector>

#include "parsing/arg_validation.hpp"

namespace lahuta::cli::validate {

namespace {

std::vector<std::string> &error_buffer() {
  static thread_local std::vector<std::string> errors;
  return errors;
}

inline std::string_view opt_name(const option::Option &opt) noexcept {
  return {opt.name, static_cast<std::size_t>(opt.namelen)};
}

} // namespace

void add_error(std::string message) { error_buffer().push_back(std::move(message)); }

void reset_errors() { error_buffer().clear(); }

bool has_errors() { return !error_buffer().empty(); }

std::vector<std::string> take_errors() {
  auto &errors = error_buffer();
  auto copy    = std::move(errors);
  errors.clear();
  return copy;
}

option::ArgStatus Unknown(const option::Option &option, bool msg) {
  if (msg) {
    std::string message = "Unknown option '";
    message.append(opt_name(option));
    message.append("'");
    error_buffer().push_back(std::move(message));
  }
  return option::ARG_ILLEGAL;
}

option::ArgStatus Required(const option::Option &option, bool msg) {
  if (option.arg != nullptr) return option::ARG_OK;
  if (msg) {
    std::string message = "Option '";
    message.append(opt_name(option));
    message.append("' requires an argument");
    error_buffer().push_back(std::move(message));
  }
  return option::ARG_ILLEGAL;
}

option::ArgStatus Verbosity(const option::Option &option, bool msg) {
  if (!option.arg) {
    if (msg) {
      std::string message = "Option '";
      message.append(opt_name(option));
      message.append("' requires a verbosity level (0, 1, or 2)");
      error_buffer().push_back(std::move(message));
    }
    return option::ARG_ILLEGAL;
  }

  const std::string_view level{option.arg};
  if (level == "0" || level == "1" || level == "2") return option::ARG_OK;

  if (msg) {
    std::string message = "Invalid verbosity level '";
    message.append(level);
    message.append("'. Must be 0 (errors only), 1 (info+), or 2 (debug+)");
    error_buffer().push_back(std::move(message));
  }
  return option::ARG_ILLEGAL;
}

} // namespace lahuta::cli::validate
