#include <string>

#include "analysis/contacts/provider.hpp"
#include "cli/arg_validation.hpp"
#include "entities/interaction_types.hpp"
#include "logging/logging.hpp"

// clang-format off
namespace lahuta::cli::validate {

namespace {
  template<typename... Args>
  void log_error(std::string_view format_str, std::string_view suffix, Args&&... args) {
    std::string full_format{format_str};
    full_format.append(suffix);
    lahuta::Logger::get_logger()->error(SPDLOG_FMT_RUNTIME(full_format), std::forward<Args>(args)...);
  }

  inline std::string_view opt_name(const option::Option& opt) noexcept {
    return {opt.name, static_cast<std::size_t>(opt.namelen) };
  }
}

option::ArgStatus Unknown(const option::Option& option, bool msg) {
  if (msg) log_error("Unknown option '{}'", HELP_MSG_SUFFIX, opt_name(option));
  return option::ARG_ILLEGAL;
}

option::ArgStatus Required(const option::Option& option, bool msg) {
  if (option.arg != nullptr) return option::ARG_OK;
  if (msg) log_error("Option '{}' requires an argument", HELP_MSG_SUFFIX, opt_name(option));
  return option::ARG_ILLEGAL;
}

option::ArgStatus Provider(const option::Option& option, bool msg) {
  if (!option.arg) {
    if (msg) log_error("Option '{}' requires a provider", HELP_MSG_SUFFIX, opt_name(option));
    return option::ARG_ILLEGAL;
  }

  const std::string_view provider{option.arg};
  if (analysis::contacts::contact_provider_from_string(provider).has_value()) {
    return option::ARG_OK;
  }

  if (msg) log_error("Invalid provider '{}'. Must be 'arpeggio', 'molstar', or 'getcontacts'", HELP_MSG_SUFFIX, provider);
  return option::ARG_ILLEGAL;
}

option::ArgStatus ContactType(const option::Option& option, bool msg) {
  if (!option.arg) {
    if (msg) log_error("Option '{}' requires a contact type", HELP_MSG_SUFFIX, opt_name(option));
    return option::ARG_ILLEGAL;
  }

  const std::string_view type{option.arg};
  if (parse_interaction_type_sequence(type, ',') || parse_interaction_type_sequence(type, '|')) {
    return option::ARG_OK;
  }

  if (msg) log_error("Invalid contact type '{}'", HELP_MSG_SUFFIX, type);
  return option::ARG_ILLEGAL;
}

option::ArgStatus Verbosity(const option::Option& option, bool msg) {
  if (!option.arg) {
    if (msg) log_error("Option '{}' requires a verbosity level (0, 1, or 2)", HELP_MSG_SUFFIX, opt_name(option));
    return option::ARG_ILLEGAL;
  }

  const std::string_view level{option.arg};
  if (level == "0" || level == "1" || level == "2") return option::ARG_OK;

  if (msg) log_error("Invalid verbosity level '{}'. Must be 0 (errors only), 1 (info+), or 2 (debug+)", HELP_MSG_SUFFIX, level);
  return option::ARG_ILLEGAL;
}

[[nodiscard]] const std::unordered_set<std::string>& get_valid_contact_types() noexcept {
  static const std::unordered_set<std::string> valid_types = {
    "all", "generic", "none",
    "hbond", "hydrogenbond", "weak_hbond", "weakhydrogenbond",
    "polar_hbond", "polarhydrogenbond", "weak_polar_hbond", "weakpolarhydrogenbond",
    "hydrophobic", "ionic", "halogen", "metalic", "metalcoordination",
    "cationpi", "aromatic", "pistacking", "pistackingp", "pistackingt",
    "carbonyl", "vdw", "vanderwaals", "donor_pi", "donorpi",
    "sulphur_pi", "sulphurpi", "carbon_pi", "carbonpi"
  };
  return valid_types;
}

} // namespace lahuta::cli::validate 
