/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   struct Overloaded {
 *     std::string& s;
 *     void operator()(const char* p) const { s += p; }
 *     void operator()(std::string_view p) const { s += p; }
 *   };
 *   std::string s;
 *   Overloaded visitor{s};
 *   visitor("besian");
 *   visitor("sejdiu");
 *   visitor(std::string_view{"@gmail.com"});
 *   return s;
 * }();
 *
 */

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

constexpr std::string_view Author = "Besian I. Sejdiu (@bisejdiu)";

class CommandSpec {
public:
  virtual ~CommandSpec() = default;

  [[nodiscard]] virtual std::string_view name() const    = 0;
  [[nodiscard]] virtual std::string_view summary() const = 0;
  [[nodiscard]] virtual std::string_view author() const { return Author; }
  [[nodiscard]] virtual const OptionSchema &schema() const = 0;

  [[nodiscard]] virtual std::any parse_config(const ParsedArgs &args) const   = 0;
  [[nodiscard]] virtual PipelinePlan build_plan(const std::any &config) const = 0;
};

using CommandSpecRegistry = std::map<std::string, const CommandSpec *>;

[[nodiscard]] const CommandSpec &get_contacts_spec() noexcept;
[[nodiscard]] const CommandSpec &get_createdb_spec() noexcept;
[[nodiscard]] const CommandSpec &get_compaction_rg_spec() noexcept;
[[nodiscard]] const CommandSpec &get_shape_metrics_spec() noexcept;
[[nodiscard]] const CommandSpec &get_extract_spec() noexcept;
[[nodiscard]] const CommandSpec &get_positions_spec() noexcept;
[[nodiscard]] const CommandSpec &get_quality_metrics_spec() noexcept;
[[nodiscard]] const CommandSpec &get_sasa_sr_spec() noexcept;
[[nodiscard]] const CommandSpec &get_dssp_spec() noexcept;

[[nodiscard]] inline const CommandSpecRegistry &get_command_specs() noexcept {
  static const CommandSpecRegistry registry = []() {
    CommandSpecRegistry map;
    const auto &createdb = get_createdb_spec();
    map.emplace(std::string(createdb.name()), &createdb);
    const auto &extract = get_extract_spec();
    map.emplace(std::string(extract.name()), &extract);
    const auto &positions = get_positions_spec();
    map.emplace(std::string(positions.name()), &positions);
    const auto &compaction_rg = get_compaction_rg_spec();
    map.emplace(std::string(compaction_rg.name()), &compaction_rg);
    const auto &shape_metrics = get_shape_metrics_spec();
    map.emplace(std::string(shape_metrics.name()), &shape_metrics);
    const auto &quality_metrics = get_quality_metrics_spec();
    map.emplace(std::string(quality_metrics.name()), &quality_metrics);
    const auto &sasa_sr = get_sasa_sr_spec();
    map.emplace(std::string(sasa_sr.name()), &sasa_sr);
    const auto &dssp = get_dssp_spec();
    map.emplace(std::string(dssp.name()), &dssp);
    const auto &contacts = get_contacts_spec();
    map.emplace(std::string(contacts.name()), &contacts);
    return map;
  }();
  return registry;
}

} // namespace lahuta::cli

#endif // LAHUTA_CLI_COMMAND_SPEC_HPP
