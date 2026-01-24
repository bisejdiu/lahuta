#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

#include <Geometry/point.h>

#include "analysis/extract/extract_tasks.hpp"
#include "cli/arg_validation.hpp"
#include "cli/extension_utils.hpp"
#include "commands/positions.hpp"
#include "commands/reporting.hpp"
#include "hash/fnv1a.hpp"
#include "logging/logging.hpp"
#include "pipeline/dynamic/manager.hpp"
#include "pipeline/dynamic/sources.hpp"
#include "runtime.hpp"

// clang-format off
namespace lahuta::cli {

namespace dyn = pipeline::dynamic;
using pipeline::DataField;
using pipeline::DataFieldSet;

struct PositionsOptions {
  enum class SourceMode { Directory, Vector, FileList, Database };

  SourceMode source_mode = SourceMode::Directory;
  std::string directory_path;
  std::vector<std::string> extensions{".cif", ".cif.gz", ".pdb", ".pdb.gz"};
  bool recursive = false;
  std::vector<std::string> file_vector;
  std::string file_list_path;
  std::string database_path;

  bool is_af2_model = false;
  std::string output_dir;
  int tree_depth = 1;
  const PipelineReporter* reporter = nullptr;
  int threads = 8;
  size_t batch_size = 512;
};

namespace {

void initialize_runtime(int num_threads) {
  LahutaRuntime::ensure_initialized(static_cast<std::size_t>(num_threads));
}

std::string sanitize_model_name(std::string_view model_path) {
  std::filesystem::path path(model_path);
  auto name = path.filename().string();
  if (!name.empty()) return name;
  return std::string(model_path);
}

std::string hex_byte(uint8_t value) {
  constexpr char kHex[] = "0123456789abcdef";
  std::string out(2, '0');
  out[0] = kHex[(value >> 4) & 0x0F];
  out[1] = kHex[value & 0x0F];
  return out;
}

std::filesystem::path
output_path_for_model(const std::filesystem::path& base_dir, const std::string& model_name, int tree_depth) {
  if (tree_depth <= 0) return base_dir / (model_name + ".bin");

  const uint64_t hash = fnv1a_64(model_name);
  const auto level1   = hex_byte(static_cast<uint8_t>((hash >> 8) & 0xFF));
  if (tree_depth == 1) return base_dir / level1 / (model_name + ".bin");

  const auto level2 = hex_byte(static_cast<uint8_t>(hash & 0xFF));
  return base_dir / level1 / level2 / (model_name + ".bin");
}

void write_positions_as_float(const RDGeom::POINT3D_VECT& pts, const std::filesystem::path& path) {
  std::ofstream ofs(path, std::ios::binary);
  if (!ofs) throw std::runtime_error("Failed to open file for writing: " + path.string());

  std::vector<float> data;
  data.reserve(pts.size() * 3);
  for (const auto& pt : pts) {
    data.push_back(static_cast<float>(pt.x));
    data.push_back(static_cast<float>(pt.y));
    data.push_back(static_cast<float>(pt.z));
  }

  ofs.write(reinterpret_cast<const char*>(data.data()), static_cast<std::streamsize>(data.size() * sizeof(float)));
}

class PositionsTask final : public dyn::ITask {
public:
  explicit PositionsTask(std::filesystem::path output_dir, int tree_depth)
      : output_dir_(std::move(output_dir)), tree_depth_(tree_depth) {}

  dyn::TaskResult run(const std::string& item_path, dyn::TaskContext& ctx) override {
    auto payload = ctx.model_payload();
    const RDGeom::POINT3D_VECT* positions = nullptr;
    std::shared_ptr<const ModelParserResult> parsed;
    if (payload && payload->positions && !payload->positions->empty()) {
      positions = payload->positions.get();
    } else if (!payload) {
      parsed = analysis::extract::get_cached_model_parser_result(ctx);
      if (parsed && parsed->coords_size() > 0) {
        positions = &parsed->coords;
      }
    }

    if (!positions || positions->empty()) {
      Logger::get_logger()->warn("[positions] Missing position data for '{}'", item_path);
      return {};
    }

    return write_positions_from_vector(item_path, *positions);
  }

  DataFieldSet data_requirements() const override {
    return DataFieldSet::of({DataField::Positions});
  }

private:
  dyn::TaskResult write_positions_from_vector(const std::string& item_path, const RDGeom::POINT3D_VECT& positions) {
    const std::string model_name = sanitize_model_name(item_path);
    if (model_name.empty()) {
      Logger::get_logger()->error("[positions] Empty model name for '{}'", item_path);
      return dyn::TaskResult{false, {}};
    }

    const auto output_path = output_path_for_model(output_dir_, model_name, tree_depth_);

    dyn::TaskResult result;
    try {
      std::filesystem::create_directories(output_path.parent_path());
      write_positions_as_float(positions, output_path);
      Logger::get_logger()->debug("[positions] Wrote {} atoms to {}", positions.size(), output_path.string());
      result.ok = true;
    } catch (const std::exception& e) {
      Logger::get_logger()->error("[positions] Failed to write '{}': {}", output_path.string(), e.what());
      result.ok = false;
    }
    return result;
  }

  std::filesystem::path output_dir_;
  int tree_depth_ = 1;
};

} // namespace

namespace positions_opts {
const option::Descriptor usage[] = {
  {PositionsOptionIndex::Unknown, 0, "", "", validate::Unknown,
   "Usage: lahuta positions --output <dir> [options]\n\n"
   "Extract 3D atomic coordinates from model files or databases.\n"
   "If you are running this on millions of files, use a --tree-depth of 2 to not hit any OS filesystem issues.\n"
   "Binary format: float32 XYZ triplets.\n"
   "Python: np.fromfile(path, dtype=np.float32).reshape(-1, 3)\n\n"
   "Input Options (choose one):"},
  {PositionsOptionIndex::Help, 0, "h", "help", option::Arg::None,
   "  --help, -h                   \tPrint this help message and exit."},
  {PositionsOptionIndex::SourceDatabase, 0, "", "database", validate::Required,
   "  --database <path>            \tProcess structures from database."},
  {PositionsOptionIndex::SourceDirectory, 0, "d", "directory", validate::Required,
   "  --directory, -d <path>       \tProcess all files in directory."},
  {PositionsOptionIndex::SourceVector, 0, "f", "files", validate::Required,
   "  --files, -f <file1,file2>    \tProcess specific files (comma-separated or repeat -f)."},
  {PositionsOptionIndex::SourceFileList, 0, "l", "file-list", validate::Required,
   "  --file-list, -l <path>       \tProcess files listed in text file (one per line)."},
  {0, 0, "", "", option::Arg::None,
   "\nDirectory Options:"},
  {PositionsOptionIndex::Extension, 0, "e", "extension", validate::Required,
   "  --extension, -e <ext>        \tFile extension(s) for directory mode. Repeat or comma-separate values (default: .cif, .cif.gz, .pdb, .pdb.gz)."},
  {PositionsOptionIndex::Recursive, 0, "r", "recursive", option::Arg::None,
   "  --recursive, -r              \tRecursively search subdirectories."},
  {0, 0, "", "", option::Arg::None,
   "\nRequired:"},
  {PositionsOptionIndex::Output, 0, "o", "output", validate::Required,
   "  --output, -o <dir>           \tOutput directory for sharded binary files."},
  {PositionsOptionIndex::TreeDepth, 0, "", "tree-depth", validate::Required,
   "  --tree-depth <0|1|2>         \tSharding depth (default: 1)."},
  {0, 0, "", "", option::Arg::None,
   "\nRuntime Options:"},
  {PositionsOptionIndex::Threads, 0, "t", "threads", validate::Required,
   "  --threads, -t <num>          \tNumber of threads to use (default: 8)."},
  {PositionsOptionIndex::BatchSize, 0, "b", "batch-size", validate::Required,
   "  --batch-size, -b <size>      \tBatch size for processing (default: 512)."},
  {PositionsOptionIndex::IsAf2Model, 0, "", "is_af2_model", option::Arg::None,
   "  --is_af2_model               \tTreat file inputs as AlphaFold2 models (ignored for --database)."},
  {0, 0, 0, 0, 0, 0}
};
} // namespace positions_opts

[[nodiscard]] std::unique_ptr<CliCommand> PositionsCommand::create() {
  return std::unique_ptr<CliCommand>(new PositionsCommand());
}

int PositionsCommand::run(int argc, char* argv[]) {
  option::Stats stats(true, positions_opts::usage, argc, const_cast<const char**>(argv));
  std::vector<option::Option> options(stats.options_max);
  std::vector<option::Option> buffer (stats.buffer_max);
  option::Parser parse(true, positions_opts::usage, argc, const_cast<const char**>(argv), options.data(), buffer.data());

  if (parse.error()) return 1;

  if (options[positions_opts::PositionsOptionIndex::Help]) {
    option::printUsage(std::cout, positions_opts::usage);
    return 0;
  }

  try {
    if (parse.nonOptionsCount() > 0) {
      Logger::get_logger()->error("Unexpected positional argument '{}'. Use --output instead.", parse.nonOption(0));
      option::printUsage(std::cerr, positions_opts::usage);
      return 1;
    }

    PositionsOptions cli;
    cli.reporter = &default_pipeline_reporter();

    // Parse source options
    int source_count = 0;
    if (options[positions_opts::PositionsOptionIndex::SourceDatabase]) {
      cli.source_mode = PositionsOptions::SourceMode::Database;
      cli.database_path = options[positions_opts::PositionsOptionIndex::SourceDatabase].arg;
      source_count++;
    }
    if (options[positions_opts::PositionsOptionIndex::SourceDirectory]) {
      cli.source_mode = PositionsOptions::SourceMode::Directory;
      cli.directory_path = options[positions_opts::PositionsOptionIndex::SourceDirectory].arg;
      source_count++;
    }
    if (options[positions_opts::PositionsOptionIndex::SourceVector]) {
      cli.source_mode = PositionsOptions::SourceMode::Vector;
      cli.file_vector.clear();
      for (const option::Option* opt = &options[positions_opts::PositionsOptionIndex::SourceVector];
           opt != nullptr;
           opt = opt->next()) {
        if (opt->arg) parse_file_argument(opt->arg, cli.file_vector);
      }
      source_count++;
    }
    if (options[positions_opts::PositionsOptionIndex::SourceFileList]) {
      cli.source_mode = PositionsOptions::SourceMode::FileList;
      cli.file_list_path = options[positions_opts::PositionsOptionIndex::SourceFileList].arg;
      source_count++;
    }

    if (source_count == 0) {
      Logger::get_logger()->error("Must specify exactly one source option: --database, --directory, --files, or --file-list");
      return 1;
    }
    if (source_count > 1) {
      Logger::get_logger()->error("Cannot specify multiple source options");
      return 1;
    }

    if (options[positions_opts::PositionsOptionIndex::Extension]) {
      cli.extensions.clear();
      for (const option::Option* opt = &options[positions_opts::PositionsOptionIndex::Extension]; opt != nullptr; opt = opt->next()) {
        parse_extension_argument(opt->arg ? opt->arg : "", cli.extensions);
      }
    }
    if (cli.extensions.empty()) cli.extensions.emplace_back();

    cli.recursive    = options[positions_opts::PositionsOptionIndex::Recursive]  ? true : false;
    cli.is_af2_model = options[positions_opts::PositionsOptionIndex::IsAf2Model] ? true : false;

    if (cli.source_mode == PositionsOptions::SourceMode::Database) {
      cli.is_af2_model = true;
    }

    if (cli.source_mode != PositionsOptions::SourceMode::Database && !cli.is_af2_model) {
      Logger::get_logger()->error(
          "positions expects AlphaFold2 model inputs. For file-based sources, pass --is_af2_model (or use --database).");
      option::printUsage(std::cerr, positions_opts::usage);
      return 1;
    }

    if (!options[positions_opts::PositionsOptionIndex::Output]) {
      Logger::get_logger()->error("Missing required --output option.");
      option::printUsage(std::cerr, positions_opts::usage);
      return 1;
    }
    std::string_view output_arg = options[positions_opts::PositionsOptionIndex::Output].arg
                                    ? options[positions_opts::PositionsOptionIndex::Output].arg
                                    : std::string_view{};
    if (output_arg.empty()) {
      Logger::get_logger()->error("--output requires a value.");
      return 1;
    }
    cli.output_dir = std::string(output_arg);

    if (options[positions_opts::PositionsOptionIndex::TreeDepth]) {
      cli.tree_depth = std::stoi(options[positions_opts::PositionsOptionIndex::TreeDepth].arg);
      if (cli.tree_depth < 0 || cli.tree_depth > 2) {
        Logger::get_logger()->error("--tree-depth must be 0, 1, or 2.");
        return 1;
      }
    }

    if (options[positions_opts::PositionsOptionIndex::Threads]) {
      cli.threads = std::stoi(options[positions_opts::PositionsOptionIndex::Threads].arg);
      if (cli.threads <= 0) {
        Logger::get_logger()->error("Threads must be positive");
        return 1;
      }
    }

    if (options[positions_opts::PositionsOptionIndex::BatchSize]) {
      cli.batch_size = std::stoull(options[positions_opts::PositionsOptionIndex::BatchSize].arg);
      if (cli.batch_size == 0) {
        Logger::get_logger()->error("Batch size must be positive");
        return 1;
      }
    }

    std::filesystem::path output_dir(cli.output_dir);
    std::error_code ec;
    if (output_dir.empty()) {
      Logger::get_logger()->error("Output directory cannot be empty.");
      return 1;
    }
    if (std::filesystem::exists(output_dir, ec)) {
      if (ec) {
        Logger::get_logger()->error("Unable to access output directory: {}", cli.output_dir);
        return 1;
      }
      if (!std::filesystem::is_directory(output_dir, ec)) {
        Logger::get_logger()->error("Output path is not a directory: {}", cli.output_dir);
        return 1;
      }
    } else {
      if (!std::filesystem::create_directories(output_dir, ec) || ec) {
        Logger::get_logger()->error("Unable to create output directory: {}", cli.output_dir);
        return 1;
      }
    }

    initialize_runtime(cli.threads);

    auto source = std::unique_ptr<sources::IDescriptor>{};
    switch (cli.source_mode) {
      case PositionsOptions::SourceMode::Database:
        source = dyn::sources_factory::from_lmdb(
            cli.database_path,
            std::string{},
            cli.batch_size,
            {static_cast<std::size_t>(cli.threads) + 1});
        break;
      case PositionsOptions::SourceMode::Directory:
        source = dyn::sources_factory::from_directory(
            cli.directory_path,
            cli.extensions,
            cli.recursive,
            cli.batch_size);
        break;
      case PositionsOptions::SourceMode::Vector:
        source = dyn::sources_factory::from_vector(cli.file_vector);
        break;
      case PositionsOptions::SourceMode::FileList:
        source = dyn::sources_factory::from_filelist(cli.file_list_path);
        break;
    }

    dyn::StageManager mgr(std::move(source));
    mgr.get_system_params().is_model =
        (cli.source_mode == PositionsOptions::SourceMode::Database) ? true : cli.is_af2_model;

    const bool needs_parse_task = cli.source_mode != PositionsOptions::SourceMode::Database;
    if (needs_parse_task) {
      auto parse_task = std::make_shared<analysis::extract::ModelParseTask>();
      mgr.add_task("parse_model", {}, std::move(parse_task), /*thread_safe=*/true);
    }

    auto task = std::make_shared<PositionsTask>(output_dir, cli.tree_depth);
    std::vector<std::string> deps;
    if (needs_parse_task) deps.emplace_back("parse_model");
    mgr.add_task("positions", std::move(deps), std::move(task), /*thread_safe=*/true);

    Logger::get_logger()->info("Positions output directory: {}", output_dir.string());
    Logger::get_logger()->info("Positions output tree depth: {}", cli.tree_depth);

    mgr.compile();
    auto progress = attach_progress_observer(mgr);
    const auto report = mgr.run(static_cast<std::size_t>(cli.threads));
    if (progress) progress->finish();
    const auto* reporter = cli.reporter ? cli.reporter : &default_pipeline_reporter();
    reporter->emit("positions", report);
    return 0;
  } catch (const std::exception& e) {
    Logger::get_logger()->error("Error: {}", e.what());
    return 1;
  }
}

} // namespace lahuta::cli
