#include <filesystem>
#include <string>
#include <string_view>
#include <utility>

#include "analysis/extract/extract_tasks.hpp"
#include "logging/logging.hpp"
#include "parsing/arg_validation.hpp"
#include "parsing/usage_error.hpp"
#include "pipeline/ingest/factory.hpp"
#include "runner/time_utils.hpp"
#include "schemas/shared_options.hpp"
#include "specs/command_spec.hpp"
#include "tasks/positions_task.hpp"

namespace lahuta::cli {
namespace A = lahuta::analysis;
namespace P = lahuta::pipeline;
namespace {

constexpr std::string_view Summary = "Extract 3D atomic coordinates from model files or databases.";

namespace positions_opts {
constexpr unsigned BaseIndex = 200;
enum : unsigned { Output = BaseIndex, TreeDepth };
} // namespace positions_opts

struct PositionsConfig {
  SourceConfig source;
  RuntimeConfig runtime;
  std::filesystem::path output_dir;
  int tree_depth = 1;
};

class PositionsSpec final : public CommandSpec {
public:
  PositionsSpec() {
    source_spec_.include_is_af2_model = true;
    source_spec_.default_extensions   = {".cif", ".cif.gz", ".pdb", ".pdb.gz"};

    runtime_spec_.include_writer_threads = false;
    runtime_spec_.default_batch_size     = 512;

    schema_.add(
        {0,
         "",
         "",
         validate::Unknown,
         std::string("Usage: lahuta positions [--output <dir>] [options]\n"
                     "Author: ")
             .append(Author)
             .append("\n\n")
             .append(Summary)
             .append(
                 "\n"
                 "If you are running this on millions of files, use a --tree-depth of 2 to not hit any OS "
                 "filesystem issues.\n"
                 "Binary format: float32 XYZ triplets.\n"
                 "Python: np.fromfile(path, dtype=np.float32).reshape(-1, 3)")});

    schema_.add({0, "", "", option::Arg::None, "\nInput Options (choose one):"});
    add_source_options(schema_, source_spec_);

    schema_.add({0, "", "", option::Arg::None, "\nOutput Options:"});
    schema_.add({positions_opts::Output,
                 "o",
                 "output",
                 validate::Required,
                 "  --output, -o <dir>           \tOutput directory for sharded binary files "
                 "(default: positions_<timestamp>)."});
    schema_.add({positions_opts::TreeDepth,
                 "",
                 "tree-depth",
                 validate::Required,
                 "  --tree-depth <0|1|2>         \tSharding depth (default: 1)."});

    schema_.add({0, "", "", option::Arg::None, "\nRuntime Options:"});
    add_runtime_options(schema_, runtime_spec_);

    schema_.add({0, "", "", option::Arg::None, "\nGlobal Options:"});
    add_global_options(schema_);
  }

  [[nodiscard]] std::string_view name() const override { return "positions"; }

  [[nodiscard]] std::string_view summary() const override { return Summary; }

  [[nodiscard]] const OptionSchema &schema() const override { return schema_; }

  [[nodiscard]] std::any parse_config(const ParsedArgs &args) const override {
    PositionsConfig config;
    config.source  = parse_source_config(args, source_spec_);
    config.runtime = parse_runtime_config(args, runtime_spec_);

    if (config.source.mode == SourceConfig::Mode::Database) {
      config.source.is_af2_model = true;
    }

    if (config.source.mode != SourceConfig::Mode::Database && !config.source.is_af2_model) {
      throw CliUsageError("positions expects AlphaFold2 model inputs. For file-based sources, pass "
                          "--is_af2_model (or use --database).");
    }

    std::string output_arg;
    if (args.has(positions_opts::Output)) {
      output_arg = args.get_string(positions_opts::Output);
      if (output_arg.empty()) {
        throw CliUsageError("--output requires a value.");
      }
    } else {
      output_arg = "positions_" + current_timestamp_string();
    }
    if (args.has(positions_opts::TreeDepth)) {
      config.tree_depth = std::stoi(args.get_string(positions_opts::TreeDepth));
      if (config.tree_depth < 0 || config.tree_depth > 2) {
        throw CliUsageError("--tree-depth must be 0, 1, or 2.");
      }
    }

    config.output_dir = validate_output_dir(output_arg);

    return std::make_any<PositionsConfig>(std::move(config));
  }

  [[nodiscard]] PipelinePlan build_plan(const std::any &config) const override {
    const auto &cfg = std::any_cast<const PositionsConfig &>(config);
    PipelinePlan plan;
    plan.report_label           = "positions";
    plan.threads                = static_cast<std::size_t>(cfg.runtime.threads);
    plan.override_system_params = true;
    plan.success_message        = "Positions extraction completed successfully!";
    plan.system_params.is_model = (cfg.source.mode == SourceConfig::Mode::Database) ? true
                                                                                    : cfg.source.is_af2_model;

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
      throw CliUsageError("positions does not support this source mode");
    };

    const bool needs_parse_task = cfg.source.mode != SourceConfig::Mode::Database;
    if (needs_parse_task) {
      PipelineTask parse_task;
      parse_task.name        = "parse_model";
      parse_task.task        = std::make_shared<A::ModelParseTask>();
      parse_task.thread_safe = true;
      plan.tasks.push_back(std::move(parse_task));
    }

    PipelineTask positions_task;
    positions_task.name = "positions";
    positions_task.task = positions::make_positions_task(cfg.output_dir, cfg.tree_depth);
    if (needs_parse_task) {
      positions_task.deps.emplace_back("parse_model");
    }
    positions_task.thread_safe = true;
    plan.tasks.push_back(std::move(positions_task));

    Logger::get_logger()->info("Positions output directory: {}", cfg.output_dir.string());
    Logger::get_logger()->info("Positions output tree depth: {}", cfg.tree_depth);

    return plan;
  }

private:
  OptionSchema schema_;
  SourceOptionSpec source_spec_;
  RuntimeOptionSpec runtime_spec_;
};

} // namespace

const CommandSpec &get_positions_spec() noexcept {
  static const PositionsSpec spec;
  return spec;
}

} // namespace lahuta::cli
