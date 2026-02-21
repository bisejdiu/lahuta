/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s;
 *   for (std::string_view part : {"besian", "sejdiu", "@gmail.com"})
 *     std::copy(part.begin(), part.end(), std::back_inserter(s));
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_CLI_USAGE_ERROR_HPP
#define LAHUTA_CLI_USAGE_ERROR_HPP

#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace lahuta::cli {

class CliUsageError final : public std::runtime_error {
public:
  explicit CliUsageError(std::string message) : std::runtime_error(message), messages_{std::move(message)} {}

  explicit CliUsageError(std::vector<std::string> messages)
      : std::runtime_error(join_messages(messages)), messages_(std::move(messages)) {}

  [[nodiscard]] const std::vector<std::string> &messages() const noexcept { return messages_; }

private:
  static std::string join_messages(const std::vector<std::string> &messages) {
    std::string joined;
    for (std::size_t i = 0; i < messages.size(); ++i) {
      if (i > 0) joined.append("\n");
      joined.append(messages[i]);
    }
    return joined;
  }

  std::vector<std::string> messages_;
};

[[nodiscard]] inline std::string usage_help_suffix(std::string_view command = {}) {
  std::string suffix = " (run lahuta";
  if (!command.empty()) {
    suffix.push_back(' ');
    suffix.append(command);
  }
  suffix.append(" -h for more information)");
  return suffix;
}

} // namespace lahuta::cli

#endif // LAHUTA_CLI_USAGE_ERROR_HPP
