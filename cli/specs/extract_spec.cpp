#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

#include "analysis/extract/extract_tasks.hpp"
#include "analysis/system/model_parse_task.hpp"
#include "logging/logging.hpp"
#include "parsing/arg_validation.hpp"
#include "parsing/extension_utils.hpp"
#include "parsing/usage_error.hpp"
#include "pipeline/runtime/api.hpp"
#include "schemas/shared_options.hpp"
#include "sinks/logging.hpp"
#include "sinks/ndjson.hpp"
#include "specs/command_spec.hpp"

namespace lahuta::cli {
namespace A = lahuta::analysis;
namespace P = lahuta::pipeline;
namespace {

constexpr std::string_view Summary = "Extract data from AlphaFold2 model files or databases.";

namespace extract_opts {
constexpr unsigned BaseIndex = 200;
enum : unsigned { Fields = BaseIndex, Output, OutputStdout };
} // namespace extract_opts

constexpr std::string_view FIELD_SEQUENCE = "sequence";
constexpr std::string_view FIELD_PLDDT    = "plddt";
constexpr std::string_view FIELD_DSSP     = "dssp";
constexpr std::string_view FIELD_ORGANISM = "organism";

bool is_valid_field(std::string_view field) {
  return field == FIELD_SEQUENCE || field == FIELD_PLDDT || field == FIELD_DSSP || field == FIELD_ORGANISM;
}

void parse_fields_argument(std::string_view raw, std::vector<std::string> &out) {
  detail::split_argument_list(raw, false, out);
}

std::shared_ptr<P::ITask> build_extract_task(std::string_view field, std::string output_channel) {
  using namespace analysis;
  if (field == FIELD_SEQUENCE) {
    return std::make_shared<SequenceExtractTask>(std::move(output_channel));
  }
  if (field == FIELD_PLDDT) {
    return std::make_shared<PlddtExtractTask>(std::move(output_channel));
  }
  if (field == FIELD_DSSP) {
    return std::make_shared<DsspExtractTask>(std::move(output_channel));
  }
  if (field == FIELD_ORGANISM) {
    return std::make_shared<OrganismExtractTask>(std::move(output_channel));
  }
  throw std::logic_error("Invalid extract field");
}

struct ExtractConfig {
  SourceConfig source;
  RuntimeConfig runtime;
  ReportConfig report;
  std::vector<std::string> fields;
  bool output_stdout   = false;
  bool output_override = false;
  std::string output_path;
};

class ExtractSpec final : public CommandSpec {
public:
  ExtractSpec() {
    source_spec_.include_is_af2_model = true;
    source_spec_.default_extensions   = {".cif", ".cif.gz", ".pdb", ".pdb.gz"};

    runtime_spec_.include_writer_threads = true;
    runtime_spec_.default_batch_size     = 512;
    runtime_spec_.default_writer_threads = P::get_default_backpressure_config().writer_threads;

    schema_.add(
        {0,
         "",
         "",
         validate::Unknown,
         std::string("Usage: lahuta extract --fields <fields> [options]\n"
                     "Author: ")
             .append(Author)
             .append("\n\n")
             .append(Summary)
             .append("\n"
                     "Note: file-based inputs must be AF2 model files (mmCIF). Generic structures are not "
                     "supported.\n\n"
                     "Available reporters (--reporter <name>):\n"
                     "  summary     - Balanced totals and item counts (default; negligible overhead).\n"
                     "  terse       - Single-line throughput summary (fastest logging footprint).\n"
                     "  detail      - Summary plus concurrency and stage breakdown details.\n\n"
                     "Fields (repeat --fields or comma-separate):\n"
                     "  sequence   Extract amino acid sequences.\n"
                     "  plddt      Extract per-residue pLDDT confidence scores.\n"
                     "  dssp       Extract secondary structure assignments (DSSP).\n"
                     "  organism   Extract organism metadata.")});

    schema_.add({0, "", "", option::Arg::None, "\nInput Options (choose one):"});
    schema_.add({extract_opts::Fields,
                 "",
                 "fields",
                 validate::Required,
                 "  --fields <field1,field2>     \tFields to extract (repeat or comma-separate)."});
    schema_.add({shared_opts::SourceDatabase,
                 "",
                 "database",
                 validate::Required,
                 "  --database <path>            \tProcess structures from database."});
    schema_.add({shared_opts::SourceDirectory,
                 "d",
                 "directory",
                 validate::Required,
                 "  --directory, -d <path>       \tProcess all files in directory."});
    schema_.add({shared_opts::SourceVector,
                 "f",
                 "files",
                 validate::Required,
                 "  --files, -f <file1,file2>    \tProcess specific files (comma-separated or repeat -f)."});
    schema_.add({shared_opts::SourceFileList,
                 "l",
                 "file-list",
                 validate::Required,
                 "  --file-list, -l <path>       \tProcess files listed in text file (one per line)."});

    schema_.add({shared_opts::SourceExtension,
                 "e",
                 "extension",
                 validate::Required,
                 "  --extension, -e <ext>        \tFile extension(s) for directory mode. Repeat or "
                 "comma-separate values (default: .cif, .cif.gz, .pdb, .pdb.gz)."});
    schema_.add({shared_opts::SourceRecursive,
                 "r",
                 "recursive",
                 option::Arg::None,
                 "  --recursive, -r              \tRecursively search subdirectories."});

    schema_.add({shared_opts::SourceIsAf2Model,
                 "",
                 "is_af2_model",
                 option::Arg::None,
                 "  --is_af2_model               \tRequired for file-based sources; inputs are AlphaFold2 "
                 "models (AF2-like mmCIF)."});

    schema_.add({0, "", "", option::Arg::None, "\nOutput Options:"});
    schema_.add({extract_opts::Output,
                 "o",
                 "output",
                 validate::Required,
                 "  --output, -o <path>          \tWrite JSONL to file (default: <field>_data.jsonl). Use "
                 "'-' for stdout."});
    schema_.add({extract_opts::OutputStdout,
                 "",
                 "stdout",
                 option::Arg::None,
                 "  --stdout                     \tWrite JSONL to stdout (same as --output -)."});
    schema_.add({0, "", "", option::Arg::None, "\nReporting Options:"});
    add_report_options(schema_);

    schema_.add({0, "", "", option::Arg::None, "\nRuntime Options:"});
    add_runtime_options(schema_, runtime_spec_);

    schema_.add({0, "", "", option::Arg::None, "\nGlobal Options:"});
    add_global_options(schema_);
  }

  [[nodiscard]] std::string_view name() const override { return "extract"; }

  [[nodiscard]] std::string_view summary() const override { return Summary; }

  [[nodiscard]] const OptionSchema &schema() const override { return schema_; }

  [[nodiscard]] std::any parse_config(const ParsedArgs &args) const override {
    ExtractConfig config;
    config.source  = parse_source_config(args, source_spec_);
    config.runtime = parse_runtime_config(args, runtime_spec_);
    config.report  = parse_report_config(args);

    (void)require_arg(args,
                      extract_opts::Fields,
                      "--fields",
                      "Missing required --fields option. Must include one or more of: sequence, plddt, dssp, "
                      "organism");

    std::vector<std::string> field_tokens;
    const auto raw_values = args.get_all_strings(extract_opts::Fields);
    for (const auto &raw : raw_values) {
      parse_fields_argument(raw, field_tokens);
    }

    if (field_tokens.empty()) {
      throw CliUsageError(
          "Missing required --fields option. Must include one or more of: sequence, plddt, dssp, organism");
    }

    std::unordered_set<std::string> seen;
    for (const auto &token : field_tokens) {
      if (!is_valid_field(token)) {
        throw CliUsageError("Invalid field '" + token + "'. Must be one of: sequence, plddt, dssp, organism");
      }
      if (seen.insert(token).second) {
        config.fields.push_back(token);
      }
    }

    config.output_stdout = args.get_flag(extract_opts::OutputStdout);

    if (args.has(extract_opts::Output)) {
      const std::string output_arg = args.get_string(extract_opts::Output);
      if (output_arg.empty()) {
        throw CliUsageError("--output requires a value.");
      }
      if (output_arg == "-") {
        config.output_stdout = true;
      } else {
        config.output_path     = output_arg;
        config.output_override = true;
      }
    }

    if (config.output_stdout && config.output_override) {
      throw CliUsageError("--stdout cannot be combined with --output <path>. Use --output - for stdout.");
    }

    if (config.fields.size() > 1 && config.output_override) {
      throw CliUsageError("--output can only be used with a single field. Omit --output for per-field "
                          "outputs or use --output - for stdout.");
    }

    if (config.source.mode == SourceConfig::Mode::Database) {
      config.source.is_af2_model = true;
    }

    if (config.source.mode != SourceConfig::Mode::Database && !config.source.is_af2_model) {
      throw CliUsageError("extract expects AlphaFold2 model inputs. For file-based sources, pass "
                          "--is_af2_model (or use --database).");
    }

    return std::make_any<ExtractConfig>(std::move(config));
  }

  [[nodiscard]] PipelinePlan build_plan(const std::any &config) const override {
    const auto &cfg = std::any_cast<const ExtractConfig &>(config);
    PipelinePlan plan;
    plan.report_label      = "extract";
    plan.reporter          = cfg.report.reporter;
    plan.save_run_report   = cfg.report.save_run_report;
    plan.run_report_prefix = "extract";
    plan.threads           = static_cast<std::size_t>(cfg.runtime.threads);
    plan.success_message   = "Extraction completed successfully!";

    plan.source_factory = [cfg]() -> PipelinePlan::SourcePtr {
      using Mode = SourceConfig::Mode;
      switch (cfg.source.mode) {
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
      throw CliUsageError("extract does not support this source mode");
    };

    auto sink_cfg           = P::get_default_backpressure_config();
    sink_cfg.writer_threads = cfg.runtime.writer_threads;

    const bool needs_parse_task = cfg.source.mode != SourceConfig::Mode::Database;
    if (needs_parse_task) {
      PipelineTask parse_task;
      parse_task.name        = "parse_model";
      parse_task.task        = std::make_shared<A::ModelParseTask>();
      parse_task.thread_safe = true;
      plan.tasks.push_back(std::move(parse_task));
    }

    std::shared_ptr<P::LoggingSink> stdout_sink;
    if (cfg.output_stdout) {
      stdout_sink = std::make_shared<P::LoggingSink>();
      Logger::get_logger()->info("Extract output sink -> stdout (-)");
    }

    for (const auto &field : cfg.fields) {
      const std::string output_channel = field + "_data";
      PipelineTask task;
      task.name = output_channel;
      task.task = build_extract_task(field, output_channel);
      if (needs_parse_task) {
        task.deps.emplace_back("parse_model");
      }
      task.thread_safe = true;
      plan.tasks.push_back(std::move(task));

      PipelineSink sink;
      sink.channel      = output_channel;
      sink.backpressure = sink_cfg;
      if (cfg.output_stdout) {
        sink.sink = stdout_sink;
      } else {
        const std::string output_path = cfg.output_override ? cfg.output_path : output_channel + ".jsonl";
        sink.sink                     = std::make_shared<P::NdjsonFileSink>(output_path);
        Logger::get_logger()->info("Extract output sink -> file: {}", output_path);
      }
      plan.sinks.push_back(std::move(sink));
    }

    return plan;
  }

private:
  OptionSchema schema_;
  SourceOptionSpec source_spec_;
  RuntimeOptionSpec runtime_spec_;
};

} // namespace

const CommandSpec &get_extract_spec() noexcept {
  static const ExtractSpec spec;
  return spec;
}

} // namespace lahuta::cli
