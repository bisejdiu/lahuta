#include <iostream>
#include <string>
#include <string_view>
#include <type_traits>
#include <variant>
#include <vector>

#include "analysis/system/model_pack_task.hpp"
#include "analysis/system/records.hpp"
#include "cli/arg_validation.hpp"
#include "cli/extension_utils.hpp"
#include "commands/createdb.hpp"
#include "commands/reporting.hpp"
#include "db/db.hpp"
#include "logging.hpp"
#include "pipeline/dynamic/manager.hpp"
#include "pipeline/dynamic/sources.hpp"
#include "runtime.hpp"
#include "sinks/lmdb.hpp"

// clang-format off
namespace lahuta::cli {

using namespace lahuta;
using namespace lahuta::pipeline;

struct CreateDbOptions {
  enum class SourceMode { Directory, Vector, FileList };

  SourceMode source_mode = SourceMode::Directory;
  std::string directory_path;
  std::vector<std::string> extensions{".cif", ".cif.gz"};
  bool recursive = true;
  std::vector<std::string> file_vector;
  std::string file_list_path;

  std::string database_path;
  size_t batch_size = 1000;
  int threads = 8;
  size_t max_size_gb = 500;
};

using WriterRes = analysis::system::ModelRecord;
using Source = std::variant<sources::Directory, std::vector<std::string>, sources::FileList>;


static Source pick_source(const CreateDbOptions& cli) {
  switch (cli.source_mode) {
    case CreateDbOptions::SourceMode::Directory:
      return sources::Directory{cli.directory_path, cli.extensions, cli.recursive, cli.batch_size};
    case CreateDbOptions::SourceMode::Vector:
      return cli.file_vector;
    case CreateDbOptions::SourceMode::FileList:
      return sources::FileList{cli.file_list_path};
  }
  throw std::logic_error("Invalid source mode");
}

static void initialize_runtime(int num_threads) {
  lahuta::LahutaRuntime::ensure_initialized(static_cast<std::size_t>(num_threads));
}

namespace createdb_opts {
const option::Descriptor usage[] = {
  {CreateDbOptionIndex::Unknown, 0, "", "", validate::Unknown,
   "Usage: lahuta createdb [options]\n\n"
   "Create a database from alphafold2 models.\n\n"
   "Input Options (choose one):"},
  {CreateDbOptionIndex::Help, 0, "h", "help", option::Arg::None,
   "  --help, -h                   \tPrint this help message and exit."},
  {CreateDbOptionIndex::SourceDirectory, 0, "d", "directory", validate::Required,
   "  --directory, -d <path>       \tProcess all files in directory."},
  {CreateDbOptionIndex::SourceVector, 0, "f", "files", validate::Required,
   "  --files, -f <file1,file2>    \tProcess specific files (comma-separated)."},
  {CreateDbOptionIndex::SourceFileList, 0, "l", "file-list", validate::Required,
   "  --file-list, -l <path>       \tProcess files listed in text file (one per line)."},
  {0, 0, "", "", option::Arg::None,
   "\nDirectory Source Options:"},
  {CreateDbOptionIndex::DatabasePath, 0, "o", "output", validate::Required,
   "  --output, -o <path>          \tOutput database path."},
  {CreateDbOptionIndex::Extension, 0, "e", "extension", validate::Required,
   "  --extension, -e <ext>        \tFile extension(s) for directory mode. Repeat or comma-separate values (default: .cif, .cif.gz)."},
  {CreateDbOptionIndex::Recursive, 0, "r", "recursive", option::Arg::None,
   "  --recursive, -r              \tRecursively search subdirectories."},
  {CreateDbOptionIndex::MaxSize, 0, "m", "max-size", validate::Required,
   "  --max-size, -m <size>        \tMaximum database size in GB (default: 500)."},
  {0, 0, "", "", option::Arg::None,
   "\nPerformance Options:"},
  {CreateDbOptionIndex::BatchSize, 0, "b", "batch-size", validate::Required,
   "  --batch-size, -b <size>      \tBatch size for database writes (default: 1000)."},
  {CreateDbOptionIndex::Threads, 0, "t", "threads", validate::Required,
   "  --threads, -t <num>          \tNumber of threads to use (default: 8)."},
  {0, 0, 0, 0, 0, 0}
};
} // namespace createdb_opts

[[nodiscard]] std::unique_ptr<CliCommand> CreateDbCommand::create() {
  return std::unique_ptr<CliCommand>(new CreateDbCommand());
}

int CreateDbCommand::run(int argc, char* argv[]) {
  option::Stats stats(true, createdb_opts::usage, argc, const_cast<const char**>(argv));
  std::vector<option::Option> options(stats.options_max);
  std::vector<option::Option> buffer (stats.buffer_max);
  option::Parser parse(true, createdb_opts::usage, argc, const_cast<const char**>(argv), options.data(), buffer.data());

  if (parse.error()) return 1;

  if (options[createdb_opts::CreateDbOptionIndex::Help]) {
    option::printUsage(std::cout, createdb_opts::usage);
    return 0;
  }

  try {
    CreateDbOptions cli;

    // parse source options
    int source_count = 0;
    if (options[createdb_opts::CreateDbOptionIndex::SourceDirectory]) {
      cli.source_mode = CreateDbOptions::SourceMode::Directory;
      cli.directory_path = options[createdb_opts::CreateDbOptionIndex::SourceDirectory].arg;
      source_count++;
    }
    if (options[createdb_opts::CreateDbOptionIndex::SourceVector]) {
      cli.source_mode = CreateDbOptions::SourceMode::Vector;
      cli.file_vector.clear();
      for (const option::Option* opt = &options[createdb_opts::CreateDbOptionIndex::SourceVector];
           opt != nullptr;
           opt = opt->next()) {
        if (opt->arg) parse_file_argument(opt->arg, cli.file_vector);
      }
      source_count++;
    }
    if (options[createdb_opts::CreateDbOptionIndex::SourceFileList]) {
      cli.source_mode = CreateDbOptions::SourceMode::FileList;
      cli.file_list_path = options[createdb_opts::CreateDbOptionIndex::SourceFileList].arg;
      source_count++;
    }

    if (source_count == 0) {
      Logger::get_logger()->error("Must specify exactly one source option: --directory, --files, or --file-list");
      return 1;
    }
    if (source_count > 1) {
      Logger::get_logger()->error("Cannot specify multiple source options");
      return 1;
    }

    if (!options[createdb_opts::CreateDbOptionIndex::DatabasePath]) {
      Logger::get_logger()->error("Database output path is required (--output)");
      return 1;
    }
    cli.database_path = options[createdb_opts::CreateDbOptionIndex::DatabasePath].arg;

    // optional options
    if (options[createdb_opts::CreateDbOptionIndex::Extension]) {
      cli.extensions.clear();
      for (const option::Option* opt = &options[createdb_opts::CreateDbOptionIndex::Extension];
           opt != nullptr;
           opt = opt->next()) {
        parse_extension_argument(opt->arg ? opt->arg : "", cli.extensions);
      }
    }
    if (cli.extensions.empty()) cli.extensions.emplace_back();

    cli.recursive = options[createdb_opts::CreateDbOptionIndex::Recursive] ? true : false;

    if (options[createdb_opts::CreateDbOptionIndex::BatchSize]) {
      cli.batch_size = std::stoull(options[createdb_opts::CreateDbOptionIndex::BatchSize].arg);
      if (cli.batch_size == 0) {
        Logger::get_logger()->error("Batch size must be positive");
        return 1;
      }
    }

    if (options[createdb_opts::CreateDbOptionIndex::Threads]) {
      cli.threads = std::stoi(options[createdb_opts::CreateDbOptionIndex::Threads].arg);
      if (cli.threads <= 0) {
        Logger::get_logger()->error("Threads must be positive");
        return 1;
      }
    }

    if (options[createdb_opts::CreateDbOptionIndex::MaxSize]) {
      cli.max_size_gb = std::stoull(options[createdb_opts::CreateDbOptionIndex::MaxSize].arg);
      if (cli.max_size_gb == 0) {
        Logger::get_logger()->error("Max size must be positive");
        return 1;
      }
    }

    Logger::get_logger()->info("Creating database...");
    switch (cli.source_mode) {
      case CreateDbOptions::SourceMode::Directory:
        Logger::get_logger()->info("Source directory: {}", cli.directory_path);
        Logger::get_logger()->info("Extensions: {}", describe_extensions(cli.extensions));
        Logger::get_logger()->info("Recursive: {}", cli.recursive ? "Yes" : "No");
        break;
      case CreateDbOptions::SourceMode::Vector:
        Logger::get_logger()->info("Source files: {} file(s)", cli.file_vector.size());
        break;
      case CreateDbOptions::SourceMode::FileList:
        Logger::get_logger()->info("Source file list: {}", cli.file_list_path);
        break;
    }
    Logger::get_logger()->info("Database path: {}", cli.database_path);
    Logger::get_logger()->info("Batch size: {}", cli.batch_size);
    Logger::get_logger()->info("Threads: {}", cli.threads);
    Logger::get_logger()->info("Max size: {} GB", cli.max_size_gb);

    initialize_runtime(cli.threads);

    auto db = std::make_shared<LMDBDatabase>(cli.database_path, cli.max_size_gb);

    Logger::get_logger()->debug("Processing files ...");

    Source source_variant = pick_source(cli);
    std::visit([&](auto&& src) {
      using SrcT = std::decay_t<decltype(src)>;
      auto source_ptr = std::unique_ptr<sources::IDescriptor>{};
      if constexpr (std::is_same_v<SrcT, sources::Directory>) {
        source_ptr = dynamic::sources_factory::from_directory(std::move(src));
      } else if constexpr (std::is_same_v<SrcT, std::vector<std::string>>) {
        source_ptr = dynamic::sources_factory::from_vector(std::move(src));
      } else if constexpr (std::is_same_v<SrcT, sources::FileList>) {
        source_ptr = dynamic::sources_factory::from_filelist(std::move(src));
      } else {
        static_assert(sizeof(SrcT) == 0, "Unsupported source type");
      }
      dynamic::StageManager mgr(std::move(source_ptr));

      auto task = std::make_shared<analysis::system::ModelPackTask>("db");
      mgr.add_task("createdb", /*deps*/{}, task, /*thread_safe=*/true);

      mgr.connect_sink("db", std::make_shared<dynamic::LmdbSink>(db, cli.batch_size));

      mgr.compile();
      const auto report = mgr.run(static_cast<std::size_t>(cli.threads));
      log_pipeline_report("createdb", report);
    }, source_variant);

    Logger::get_logger()->info("Database creation completed successfully!");
    return 0;

  } catch (const std::exception& e) {
    Logger::get_logger()->error("Error: {}", e.what());
    return 1;
  }
}

} // namespace lahuta::cli
