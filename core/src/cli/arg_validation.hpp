#ifndef LAHUTA_CLI_ARG_VALIDATION_HPP
#define LAHUTA_CLI_ARG_VALIDATION_HPP

#include "gemmi/third_party/optionparser.h"
#include <unordered_set>
#include <string>
#include <string_view>

namespace lahuta::cli::validate {

inline constexpr std::string_view HELP_MSG_SUFFIX = " (run lahuta -h for more information)";
inline constexpr std::string_view CONTACTS_HELP_MSG_SUFFIX = " (run lahuta contacts -h for more information)";

option::ArgStatus Unknown    (const option::Option& option, bool msg); // unknwn options
option::ArgStatus Required   (const option::Option& option, bool msg); // options that require an argument
option::ArgStatus Provider   (const option::Option& option, bool msg); // provider options
option::ArgStatus ContactType(const option::Option& option, bool msg); // contact type options
option::ArgStatus Verbosity  (const option::Option& option, bool msg); // verbosity level options
inline option::ArgStatus Ignore(const option::Option&, bool) { return option::ARG_IGNORE; }

// Get the set of valid contact types for validation.
[[nodiscard]] const std::unordered_set<std::string>& get_valid_contact_types() noexcept;

} // namespace lahuta::cli::validate

#endif // LAHUTA_CLI_ARG_VALIDATION_HPP
