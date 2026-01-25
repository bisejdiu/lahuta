#ifndef LAHUTA_CLI_OPTION_SCHEMA_HPP
#define LAHUTA_CLI_OPTION_SCHEMA_HPP

#include <string>
#include <vector>

#include <gemmi/third_party/optionparser.h>

namespace lahuta::cli {

struct OptionDef {
  unsigned index = 0;
  std::string short_name;
  std::string long_name;
  option::ArgStatus (*validator)(const option::Option &, bool) = option::Arg::None;
  std::string help;
};

class OptionSchema {
public:
  void add(OptionDef def);
  void append(const OptionSchema &other);
  [[nodiscard]] const std::vector<OptionDef> &defs() const noexcept;
  [[nodiscard]] std::vector<option::Descriptor> build_descriptors() const;

private:
  std::vector<OptionDef> defs_;
};

} // namespace lahuta::cli

#endif // LAHUTA_CLI_OPTION_SCHEMA_HPP
