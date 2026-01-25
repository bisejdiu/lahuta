#ifndef LAHUTA_CLI_PARSED_ARGS_HPP
#define LAHUTA_CLI_PARSED_ARGS_HPP

#include <string>
#include <vector>

#include <gemmi/third_party/optionparser.h>

namespace lahuta::cli {

class ParsedArgs {
public:
  ParsedArgs() = default;
  ParsedArgs(const option::Parser &parser, const option::Option *options);

  [[nodiscard]] bool has(int key) const noexcept;
  [[nodiscard]] std::string get_string(int key) const;
  [[nodiscard]] std::vector<std::string> get_all_strings(int key) const;
  [[nodiscard]] bool get_flag(int key) const noexcept;

private:
  const option::Parser *parser_  = nullptr;
  const option::Option *options_ = nullptr;
};

} // namespace lahuta::cli

#endif // LAHUTA_CLI_PARSED_ARGS_HPP
