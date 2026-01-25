#ifndef LAHUTA_CLI_COMMAND_SPEC_HPP
#define LAHUTA_CLI_COMMAND_SPEC_HPP

#include <any>
#include <map>
#include <string>
#include <string_view>

#include "parsing/parsed_args.hpp"
#include "runner/pipeline_plan.hpp"
#include "schemas/option_schema.hpp"

namespace lahuta::cli {

class CommandSpec {
public:
  virtual ~CommandSpec() = default;

  [[nodiscard]] virtual std::string_view name() const      = 0;
  [[nodiscard]] virtual std::string_view summary() const   = 0;
  [[nodiscard]] virtual const OptionSchema &schema() const = 0;

  [[nodiscard]] virtual std::any parse_config(const ParsedArgs &args) const   = 0;
  [[nodiscard]] virtual PipelinePlan build_plan(const std::any &config) const = 0;
};

using CommandSpecRegistry = std::map<std::string, const CommandSpec *>;

[[nodiscard]] const CommandSpec &get_contacts_spec() noexcept;
[[nodiscard]] const CommandSpec &get_createdb_spec() noexcept;
[[nodiscard]] const CommandSpec &get_extract_spec() noexcept;
[[nodiscard]] const CommandSpec &get_positions_spec() noexcept;
[[nodiscard]] const CommandSpec &get_quality_metrics_spec() noexcept;

[[nodiscard]] inline const CommandSpecRegistry &get_command_specs() noexcept {
  static const CommandSpecRegistry registry = []() {
    CommandSpecRegistry map;
    const auto &createdb = get_createdb_spec();
    map.emplace(std::string(createdb.name()), &createdb);
    const auto &extract = get_extract_spec();
    map.emplace(std::string(extract.name()), &extract);
    const auto &positions = get_positions_spec();
    map.emplace(std::string(positions.name()), &positions);
    const auto &quality_metrics = get_quality_metrics_spec();
    map.emplace(std::string(quality_metrics.name()), &quality_metrics);
    const auto &contacts = get_contacts_spec();
    map.emplace(std::string(contacts.name()), &contacts);
    return map;
  }();
  return registry;
}

} // namespace lahuta::cli

#endif // LAHUTA_CLI_COMMAND_SPEC_HPP
