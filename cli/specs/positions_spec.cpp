#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include "analysis/extract/extract_tasks.hpp"
#include "logging/logging.hpp"
#include "parsing/arg_validation.hpp"
#include "pipeline/ingest/factory.hpp"
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
    runtime_spec_.default_threads        = 8;
    runtime_spec_.default_batch_size     = 512;

    schema_.add(
        {0,
         "",
         "",
         validate::Unknown,
         std::string("Usage: lahuta positions --output <dir> [options]\n\n")
             .append(Summary)
             .append(
                 "\n"
                 "If you are running this on millions of files, use a --tree-depth of 2 to not hit any OS "
                 "filesystem issues.\n"
                 "Binary format: float32 XYZ triplets.\n"
                 "Python: np.fromfile(path, dtype=np.float32).reshape(-1, 3)")});

    schema_.add({0, "", "", option::Arg::None, "\nInput Options (choose one):"});
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

    schema_.add({0, "", "", option::Arg::None, "\nDirectory Options:"});
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

    schema_.add({0, "", "", option::Arg::None, "\nRequired:"});
    schema_.add({positions_opts::Output,
                 "o",
                 "output",
                 validate::Required,
                 "  --output, -o <dir>           \tOutput directory for sharded binary files."});
    schema_.add({positions_opts::TreeDepth,
                 "",
                 "tree-depth",
                 validate::Required,
                 "  --tree-depth <0|1|2>         \tSharding depth (default: 1)."});

    schema_.add({0, "", "", option::Arg::None, "\nRuntime Options:"});
    add_runtime_options(schema_, runtime_spec_);
    schema_.add({shared_opts::SourceIsAf2Model,
                 "",
                 "is_af2_model",
                 option::Arg::None,
                 "  --is_af2_model               \tTreat file inputs as AlphaFold2 models (ignored for "
                 "--database)."});

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
      throw std::runtime_error("positions expects AlphaFold2 model inputs. For file-based sources, pass "
                               "--is_af2_model (or use --database).");
    }

    if (!args.has(positions_opts::Output)) {
      throw std::runtime_error("Missing required --output option.");
    }

    const std::string output_arg = args.get_string(positions_opts::Output);
    if (output_arg.empty()) {
      throw std::runtime_error("--output requires a value.");
    }
    config.output_dir = std::filesystem::path(output_arg);

    if (args.has(positions_opts::TreeDepth)) {
      config.tree_depth = std::stoi(args.get_string(positions_opts::TreeDepth));
      if (config.tree_depth < 0 || config.tree_depth > 2) {
        throw std::runtime_error("--tree-depth must be 0, 1, or 2.");
      }
    }

    if (config.output_dir.empty()) {
      throw std::runtime_error("Output directory cannot be empty.");
    }

    std::error_code ec;
    if (std::filesystem::exists(config.output_dir, ec)) {
      if (ec) {
        throw std::runtime_error("Unable to access output directory: " + output_arg);
      }
      if (!std::filesystem::is_directory(config.output_dir, ec)) {
        throw std::runtime_error("Output path is not a directory: " + output_arg);
      }
    } else {
      if (!std::filesystem::create_directories(config.output_dir, ec) || ec) {
        throw std::runtime_error("Unable to create output directory: " + output_arg);
      }
    }

    return std::make_any<PositionsConfig>(std::move(config));
  }

  [[nodiscard]] PipelinePlan build_plan(const std::any &config) const override {
    const auto &cfg = std::any_cast<const PositionsConfig &>(config);
    PipelinePlan plan;
    plan.report_label           = "positions";
    plan.threads                = static_cast<std::size_t>(cfg.runtime.threads);
    plan.override_system_params = true;
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
      throw std::runtime_error("positions does not support this source mode");
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
