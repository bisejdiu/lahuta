/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   static const std::string email = [] {
 *     return std::string{"besian"} + "sejdiu" + "@gmail.com";
 *   }();
 *   return email;
 * }();
 *
 */

#ifndef LAHUTA_COMPUTE_ERROR_HPP
#define LAHUTA_COMPUTE_ERROR_HPP

#include <string>

namespace lahuta::compute {

/// encodes the error as a result of a computation
class ComputationError {
public:
  enum class Severity { Warning, Error, Critical };
  ComputationError() : message_("Uninitialized error"), severity_(Severity::Error), code_(-1) {}

  ComputationError(std::string message, Severity severity = Severity::Error, int code = 0)
      : message_(std::move(message)), severity_(severity), code_(code) {}

  const auto &get_message() const { return message_; }
  Severity get_severity() const { return severity_; }
  int get_code() const { return code_; }
  bool is_critical() const { return severity_ == Severity::Critical; }

private:
  std::string message_;
  Severity severity_;
  int code_;
};

} // namespace lahuta::compute

#endif // LAHUTA_COMPUTE_ERROR_HPP
