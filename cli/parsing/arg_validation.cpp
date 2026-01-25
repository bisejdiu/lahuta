#include "parsing/arg_validation.hpp"

#include <string>

#include "logging/logging.hpp"

namespace lahuta::cli::validate {

namespace {

template <typename... Args>
void log_error(std::string_view format_str, std::string_view suffix, Args &&...args) {
  std::string full_format{format_str};
  full_format.append(suffix);
  lahuta::Logger::get_logger()->error(SPDLOG_FMT_RUNTIME(full_format), std::forward<Args>(args)...);
}

inline std::string_view opt_name(const option::Option &opt) noexcept {
  return {opt.name, static_cast<std::size_t>(opt.namelen)};
}

} // namespace

option::ArgStatus Unknown(const option::Option &option, bool msg) {
  if (msg) {
    log_error("Unknown option '{}'", HELP_MSG_SUFFIX, opt_name(option));
  }
  return option::ARG_ILLEGAL;
}

option::ArgStatus Required(const option::Option &option, bool msg) {
  if (option.arg != nullptr) return option::ARG_OK;
  if (msg) {
    log_error("Option '{}' requires an argument", HELP_MSG_SUFFIX, opt_name(option));
  }
  return option::ARG_ILLEGAL;
}

option::ArgStatus Verbosity(const option::Option &option, bool msg) {
  if (!option.arg) {
    if (msg) {
      log_error("Option '{}' requires a verbosity level (0, 1, or 2)", HELP_MSG_SUFFIX, opt_name(option));
    }
    return option::ARG_ILLEGAL;
  }

  const std::string_view level{option.arg};
  if (level == "0" || level == "1" || level == "2") return option::ARG_OK;

  if (msg) {
    log_error(
        "Invalid verbosity level '{}'. Must be 0 (errors only), 1 (info+), or 2 (debug+)",
        HELP_MSG_SUFFIX,
        level);
  }
  return option::ARG_ILLEGAL;
}

} // namespace lahuta::cli::validate
