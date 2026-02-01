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

void reset_errors();
[[nodiscard]] bool has_errors();
[[nodiscard]] std::vector<std::string> take_errors();

} // namespace lahuta::cli::validate

#endif // LAHUTA_CLI_ARG_VALIDATION_HPP
