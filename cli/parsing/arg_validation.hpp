#ifndef LAHUTA_CLI_ARG_VALIDATION_HPP
#define LAHUTA_CLI_ARG_VALIDATION_HPP

#include <string_view>

#include <gemmi/third_party/optionparser.h>

namespace lahuta::cli::validate {

inline constexpr std::string_view HELP_MSG_SUFFIX = " (run lahuta -h for more information)";

option::ArgStatus Unknown(const option::Option &option, bool msg);
option::ArgStatus Required(const option::Option &option, bool msg);
option::ArgStatus Verbosity(const option::Option &option, bool msg);
inline option::ArgStatus Ignore(const option::Option &, bool) { return option::ARG_IGNORE; }

} // namespace lahuta::cli::validate

#endif // LAHUTA_CLI_ARG_VALIDATION_HPP
