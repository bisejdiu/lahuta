/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto append_if_string = [](std::string& s, auto&& arg)
 *       -> std::enable_if_t<std::is_convertible_v<decltype(arg), std::string_view>> {
 *     s += arg;
 *   };
 *   std::string s;
 *   append_if_string(s, "besian");
 *   append_if_string(s, "sejdiu");
 *   append_if_string(s, "@gmail.com");
 *   return s;
 * }();
 *
 */

#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "analysis/contacts/computation.hpp"
#include "analysis/contacts/provider.hpp"
#include "entities/interaction_types.hpp"
#include "logging/logging.hpp"
#include "parsing/arg_validation.hpp"
#include "parsing/extension_utils.hpp"
#include "parsing/usage_error.hpp"
#include "pipeline/runtime/api.hpp"
#include "pipeline/task/api.hpp"
#include "schemas/shared_options.hpp"
#include "sinks/logging.hpp"
#include "sinks/ndjson.hpp"
#include "specs/command_spec.hpp"
#include "tasks/contacts_md_source.hpp"

namespace lahuta::cli {
namespace A = lahuta::analysis;
namespace P = lahuta::pipeline;
namespace {

constexpr std::string_view Summary = "Compute inter-atomic contacts.";

// Contacts-specific argument validators
namespace contacts_validate {

inline std::string_view opt_name(const option::Option &opt) noexcept {
  return {opt.name, static_cast<std::size_t>(opt.namelen)};
}

option::ArgStatus Provider(const option::Option &option, bool msg) {
  if (!option.arg) {
    if (msg) {
      std::string message = "Option '";
      message.append(opt_name(option));
      message.append("' requires a provider");
      validate::add_error(std::move(message));
    }
    return option::ARG_ILLEGAL;
  }

  const std::string_view provider{option.arg};
  if (A::contact_provider_from_string(provider).has_value()) {
    return option::ARG_OK;
  }

  if (msg) {
    std::string message = "Invalid provider '";
    message.append(provider);
    message.append("'. Must be 'arpeggio', 'molstar', or 'getcontacts'");
    validate::add_error(std::move(message));
  }
  return option::ARG_ILLEGAL;
}

option::ArgStatus ContactType(const option::Option &option, bool msg) {
  if (!option.arg) {
    if (msg) {
      std::string message = "Option '";
      message.append(opt_name(option));
      message.append("' requires a contact type");
      validate::add_error(std::move(message));
    }
    return option::ARG_ILLEGAL;
  }

  const std::string_view type{option.arg};
  if (parse_interaction_type_sequence(type, ',') || parse_interaction_type_sequence(type, '|')) {
    return option::ARG_OK;
  }

  if (msg) {
    std::string message = "Invalid contact type '";
    message.append(type);
    message.append("'");
    validate::add_error(std::move(message));
  }
  return option::ARG_ILLEGAL;
}

} // namespace contacts_validate

namespace contacts_opts {
constexpr unsigned BaseIndex = 200;
enum : unsigned { //
  SourceMD = BaseIndex,
  Provider,
  InteractionType,
  MolstarWaterWaterHbonds,
  OutputPath,
  OutputStdout,
  OutputCompact
};
} // namespace contacts_opts

enum class ContactsSourceMode { Directory, Vector, FileList, Database, MD };

struct ContactsConfig {
  ContactsSourceMode source_mode = ContactsSourceMode::Directory;
  SourceConfig source;
  contacts::MdInputs md_inputs;
  RuntimeConfig runtime;
  ReportConfig report;
  A::ContactProvider provider          = A::ContactProvider::MolStar;
  InteractionTypeSet interaction_types = InteractionTypeSet::all();
  bool molstar_include_water_water_hbonds = false;
  bool is_af2_model                    = false;
  bool output_stdout                   = false;
  bool output_compact                  = false;
  std::string output_path;
};

class ContactsSpec final : public CommandSpec {
public:
  ContactsSpec() {
    source_spec_.include_is_af2_model = true;
    source_spec_.default_extensions   = {".cif", ".cif.gz"};

    runtime_spec_.include_writer_threads = true;
    runtime_spec_.default_batch_size     = 200;
    runtime_spec_.default_writer_threads = P::get_default_backpressure_config().writer_threads;

    schema_.add({0,
                 "",
                 "",
                 validate::Unknown,
                 std::string("Usage: lahuta contacts [options]\n"
                             "Author: ")
                     .append(Author)
                     .append("\n\n")
                     .append(Summary)
                     .append("\n\n"
                             "<<< Detailed help will be added soon. >>>\n\n"
                             "Available reporters (--reporter <name>):\n"
                             "  summary     - Concise summary of computation results (negligible overhead).\n"
                             "  terse       - Fastest logging footprint.\n"
                             "  detail      - Detailed diagnostics for in-depth analysis (tiny overhead).")});

    schema_.add({0, "", "", option::Arg::None, "\nInput Options (choose one):"});
    add_source_options(schema_, source_spec_);
    schema_.add({contacts_opts::SourceMD,
                 "",
                 "md",
                 validate::Required,
                 "  --md <path>                  \tMD input file. Specify once for structure (PDB/GRO) and "
                 "once or more for trajectories (XTC)."});

    schema_.add({0, "", "", option::Arg::None, "\nOutput Options:"});
    schema_.add({contacts_opts::OutputPath,
                 "o",
                 "output",
                 validate::Required,
                 "  --output, -o <path>          \tWrite JSONL to file (default: "
                 "contacts.jsonl). Use '-' for stdout."});
    schema_.add({contacts_opts::OutputStdout,
                 "",
                 "stdout",
                 option::Arg::None,
                 "  --stdout                     \tWrite JSONL to stdout (same as --output -)."});
    schema_.add({contacts_opts::OutputCompact,
                 "",
                 "compact",
                 option::Arg::None,
                 "  --compact                    \tWrite compact JSONL (short keys, grouped contacts, 2dp "
                 "distances)."});

    schema_.add({0, "", "", option::Arg::None, "\nCompute Options:"});
    schema_.add({contacts_opts::Provider,
                 "p",
                 "provider",
                 contacts_validate::Provider,
                 "  --provider, -p <provider>    \tContact provider: 'molstar', 'arpeggio', or 'getcontacts' "
                 "(default: molstar)."});
    schema_.add({contacts_opts::InteractionType,
                 "i",
                 "interaction",
                 contacts_validate::ContactType,
                 "  --interaction, -i <type>     \tInteraction type(s): 'hbond', 'hydrophobic', 'ionic', "
                 "etc. Repeat or comma-separate to combine.\n"});
    schema_.add({contacts_opts::MolstarWaterWaterHbonds,
                 "",
                 "ms-water-water-hbonds",
                 option::Arg::None,
                 "  --ms-water-water-hbonds\n"
                 "                               \tInclude Mol* water-water hydrogen bonds (default: false).\n"});

    schema_.add({0, "", "", option::Arg::None, "\nReporting Options:"});
    add_report_options(schema_);

    schema_.add({0, "", "", option::Arg::None, "\nRuntime Options:"});
    add_runtime_options(schema_, runtime_spec_);

    schema_.add({0, "", "", option::Arg::None, "\nGlobal Options:"});
    add_global_options(schema_);
  }

  [[nodiscard]] std::string_view name() const override { return "contacts"; }
  [[nodiscard]] std::string_view summary() const override { return Summary; }
  [[nodiscard]] const OptionSchema &schema() const override { return schema_; }

  [[nodiscard]] std::any parse_config(const ParsedArgs &args) const override {
    ContactsConfig config;
    config.runtime      = parse_runtime_config(args, runtime_spec_);
    config.report       = parse_report_config(args);
    config.is_af2_model = args.get_flag(shared_opts::SourceIsAf2Model);

    const bool has_md              = args.has(contacts_opts::SourceMD);
    const bool has_standard_source = args.has(shared_opts::SourceDirectory) ||
                                     args.has(shared_opts::SourceVector) ||
                                     args.has(shared_opts::SourceFileList) ||
                                     args.has(shared_opts::SourceDatabase);

    if (!has_md && !has_standard_source) {
      throw CliUsageError(
          "Must specify exactly one source option: --directory, --files, --file-list, --database, or --md");
    }
    if (has_md && has_standard_source) {
      throw CliUsageError("Cannot specify multiple source options");
    }

    if (has_md) {
      const auto md_files = args.get_all_strings(contacts_opts::SourceMD);
      config.md_inputs    = contacts::parse_md_inputs(md_files);
      config.source_mode  = ContactsSourceMode::MD;
    } else {
      config.source       = parse_source_config(args, source_spec_);
      config.is_af2_model = config.source.is_af2_model;
      switch (config.source.mode) {
        case SourceConfig::Mode::Directory:
          config.source_mode = ContactsSourceMode::Directory;
          break;
        case SourceConfig::Mode::Vector:
          config.source_mode = ContactsSourceMode::Vector;
          break;
        case SourceConfig::Mode::FileList:
          config.source_mode = ContactsSourceMode::FileList;
          break;
        case SourceConfig::Mode::Database:
          config.source_mode = ContactsSourceMode::Database;
          break;
      }
    }

    if (args.has(contacts_opts::Provider)) {
      const std::string provider_arg = args.get_string(contacts_opts::Provider);
      if (auto provider = A::contact_provider_from_string(provider_arg)) {
        config.provider = *provider;
      } else {
        throw CliUsageError("Invalid provider '" + provider_arg +
                            "'. Must be 'molstar', 'arpeggio', or 'getcontacts'");
      }
    }

    if (args.has(contacts_opts::InteractionType)) {
      bool has_selection = false;
      InteractionTypeSet selected;
      const auto raw_values = args.get_all_strings(contacts_opts::InteractionType);
      for (const auto &raw : raw_values) {
        auto parsed = parse_interaction_type_sequence(raw, ',');
        if (!parsed) parsed = parse_interaction_type_sequence(raw, '|');
        if (!parsed) {
          throw CliUsageError("Invalid interaction type specification '" + raw + "'");
        }
        if (!has_selection) {
          selected      = *parsed;
          has_selection = true;
        } else {
          selected |= *parsed;
        }
      }
      if (has_selection) {
        config.interaction_types = selected;
      }
    }

    if (args.get_flag(contacts_opts::MolstarWaterWaterHbonds)) {
      config.molstar_include_water_water_hbonds = true;
    }

    config.output_stdout = args.get_flag(contacts_opts::OutputStdout);
    config.output_compact = args.get_flag(contacts_opts::OutputCompact);

    std::string output_arg;
    if (args.has(contacts_opts::OutputPath)) {
      output_arg = args.get_string(contacts_opts::OutputPath);
      if (output_arg.empty()) {
        throw CliUsageError("--output requires a value.");
      }
      if (output_arg == "-") {
        config.output_stdout = true;
      } else {
        config.output_path = output_arg;
      }
    }

    if (config.output_stdout && !output_arg.empty() && output_arg != "-") {
      throw CliUsageError("--stdout cannot be combined with --output <path>. Use --output - for stdout.");
    }

    if (!config.output_stdout && config.output_path.empty()) {
      config.output_path = "contacts.jsonl";
    }

    return std::make_any<ContactsConfig>(std::move(config));
  }

  [[nodiscard]] PipelinePlan build_plan(const std::any &config) const override {
    const auto &cfg = std::any_cast<const ContactsConfig &>(config);
    PipelinePlan plan;
    plan.report_label                       = "contacts";
    plan.reporter                           = cfg.report.reporter;
    plan.save_run_report                    = cfg.report.save_run_report;
    plan.run_report_prefix                  = "contacts";
    plan.threads                            = static_cast<std::size_t>(cfg.runtime.threads);
    plan.auto_builtins                      = true;
    plan.override_topology_params           = true;
    plan.topology_params.atom_typing_method = A::typing_for_provider(cfg.provider);
    plan.success_message                    = "Contact computation completed successfully!";

    const bool is_md              = cfg.source_mode == ContactsSourceMode::MD;
    const bool is_db              = cfg.source_mode == ContactsSourceMode::Database;
    const bool use_model_pipeline = cfg.is_af2_model || is_db;
    if (!is_md) {
      plan.override_system_params = true;
      plan.system_params.is_model = use_model_pipeline;
    }

    if (cfg.source_mode == ContactsSourceMode::Directory) {
      Logger::get_logger()->info("Source directory: {}", cfg.source.directory_path);
      Logger::get_logger()->info("Extensions: {}", describe_extensions(cfg.source.extensions));
      Logger::get_logger()->info("Recursive: {}", cfg.source.recursive ? "Yes" : "No");
    } else if (cfg.source_mode == ContactsSourceMode::MD) {
      Logger::get_logger()->info("MD Structure file: {}", cfg.md_inputs.structure_path);
      Logger::get_logger()->info("MD Trajectory files ({}):", cfg.md_inputs.trajectory_paths.size());
      for (const auto &xtc : cfg.md_inputs.trajectory_paths) {
        Logger::get_logger()->info("  - {}", xtc);
      }
    }

    plan.source_factory = [cfg]() -> PipelinePlan::SourcePtr {
      using Mode = ContactsSourceMode;
      switch (cfg.source_mode) {
        case Mode::MD: {
          return contacts::make_md_source_descriptor(cfg.md_inputs.structure_path,
                                                     cfg.md_inputs.trajectory_paths);
        }
        case Mode::Database: {
          auto source = P::from_lmdb(cfg.source.database_path,
                                     std::string{},
                                     cfg.runtime.batch_size,
                                     {static_cast<std::size_t>(cfg.runtime.threads) + 1});
          return PipelinePlan::SourcePtr(std::move(source));
        }
        case Mode::Directory: {
          auto source = P::from_directory(cfg.source.directory_path,
                                          cfg.source.extensions,
                                          cfg.source.recursive,
                                          cfg.runtime.batch_size);
          return PipelinePlan::SourcePtr(std::move(source));
        }
        case Mode::Vector: {
          auto source = P::from_vector(cfg.source.file_vector);
          return PipelinePlan::SourcePtr(std::move(source));
        }
        case Mode::FileList: {
          auto source = P::from_filelist(cfg.source.file_list_path);
          return PipelinePlan::SourcePtr(std::move(source));
        }
      }
      throw CliUsageError("contacts does not support this source mode");
    };

    P::ContactsParams params;
    params.provider = cfg.provider;
    params.type     = cfg.interaction_types;
    params.molstar_include_water_water_hbonds = cfg.molstar_include_water_water_hbonds;
    params.channel  = "contacts";
    params.format   = cfg.output_compact ? P::ContactsOutputFormat::JsonCompact : P::ContactsOutputFormat::Json;

    PipelineComputation computation;
    computation.name    = "contacts";
    computation.factory = [params]() { return std::make_unique<A::ContactsComputation>("contacts", params); };
    computation.thread_safe = true;
    plan.computations.push_back(std::move(computation));

    auto sink_cfg           = P::get_default_backpressure_config();
    sink_cfg.writer_threads = cfg.runtime.writer_threads;

    PipelineSink sink;
    sink.channel      = "contacts";
    sink.backpressure = sink_cfg;
    if (cfg.output_stdout) {
      sink.sink = std::make_shared<P::LoggingSink>();
      Logger::get_logger()->info("Writing to: stdout");
    } else {
      sink.sink = std::make_shared<P::NdjsonFileSink>(cfg.output_path);
      Logger::get_logger()->info("Writing to: {}", cfg.output_path);
      plan.output_files.push_back(cfg.output_path);
    }
    plan.sinks.push_back(std::move(sink));

    return plan;
  }

private:
  OptionSchema schema_;
  SourceOptionSpec source_spec_;
  RuntimeOptionSpec runtime_spec_;
};

} // namespace

const CommandSpec &get_contacts_spec() noexcept {
  static const ContactsSpec spec;
  return spec;
}

} // namespace lahuta::cli
