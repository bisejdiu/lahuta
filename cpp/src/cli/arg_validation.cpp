#include "cli/arg_validation.hpp"
#include "logging.hpp"
#include <string>

// clang-format off
namespace lahuta::cli::validate {

namespace {
  template<typename... Args>
  void log_error(const std::string& format_str, std::string_view suffix, Args&&... args) {
    const auto full_format = format_str + std::string{suffix};
    lahuta::Logger::get_logger()->error(full_format, std::forward<Args>(args)...);
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
    if (msg) log_error("Option '{}' requires a provider (arpeggio or molstar)", HELP_MSG_SUFFIX, opt_name(option));
    return option::ARG_ILLEGAL;
  }

  const std::string_view provider{option.arg};
  if (provider == "arpeggio" || provider == "molstar") return option::ARG_OK;

  if (msg) log_error("Invalid provider '{}'. Must be 'arpeggio' or 'molstar'", HELP_MSG_SUFFIX, provider);
  return option::ARG_ILLEGAL;
}


option::ArgStatus ContactType(const option::Option& option, bool msg) {
  if (!option.arg) {
    if (msg) log_error("Option '{}' requires a contact type", HELP_MSG_SUFFIX, opt_name(option));
    return option::ARG_ILLEGAL;
  }

  const std::string_view type{option.arg};
  const auto& valid_types = get_valid_contact_types();

  if (const auto it = valid_types.find(std::string{type}); it != valid_types.end()) return option::ARG_OK;

  if (msg) log_error("Invalid contact type '{}'", HELP_MSG_SUFFIX, type);
  return option::ARG_ILLEGAL;
}

[[nodiscard]] const std::unordered_set<std::string>& get_valid_contact_types() noexcept {
  static const std::unordered_set<std::string> valid_types = {
    "hbond", "hydrophobic", "ionic", // these three are common
    "weak_hbond", "halogen", "metalic", "cationpi", "pistacking",
    "polar_hbond", "weak_polar_hbond", "aromatic", "carbonyl",
    "vdw", "donor_pi", "sulphur_pi", "carbon_pi"
  };
  return valid_types;
}

} // namespace lahuta::cli::validate 
