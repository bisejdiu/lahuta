#include <chrono>
#include <iostream>
#include <sstream>
#include <string>
#include <variant>
#include <vector>

#include "analysis/system/records.hpp"
#include "cli/arg_validation.hpp"
#include "commands/createdb.hpp"
#include "db/db.hpp"
#include "gemmi/gz.hpp"
#include "io/collector.hpp"
#include "io/db_spill_policy.hpp"
#include "io/lmdb_backend.hpp"
#include "logging.hpp"
#include "mmap/MemoryMapped.h"
#include "runtime.hpp"
#include "models/topology.hpp"
#include "pipeline/dsl.hpp"
#include "serialization/formats.hpp"

// clang-format off
namespace lahuta::cli {

using namespace lahuta;
using namespace lahuta::pipeline;

struct CreateDbOptions {
  enum class SourceMode { Directory, Vector, FileList };

  SourceMode source_mode = SourceMode::Directory;
  std::string directory_path;
  std::string extension = ".cif.gz";
  bool recursive = true;
  std::vector<std::string> file_vector;
  std::string file_list_path;

  std::string database_path;
  size_t batch_size = 1000;
  int threads = 8;
};

using WriterRes = analysis::system::ModelRecord;
using Payload = std::shared_ptr<const WriterRes>;
using Source = std::variant<sources::DirectorySource, sources::VectorSource, sources::FileListSource>;

using LMDBPolicy = DBSpillPolicy<fmt::binary, Payload, LMDBBackend>;

static Source pick_source(const CreateDbOptions& cli) {
  switch (cli.source_mode) {
    case CreateDbOptions::SourceMode::Directory:
      return sources::DirectorySource{cli.directory_path, cli.extension, cli.recursive, cli.batch_size};
    case CreateDbOptions::SourceMode::Vector:
      return sources::VectorSource{cli.file_vector};
    case CreateDbOptions::SourceMode::FileList:
      return sources::FileListSource{cli.file_list_path};
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
   "  --extension, -e <ext>        \tFile extension for directory mode (default: .cif.gz)."},
  {CreateDbOptionIndex::Recursive, 0, "r", "recursive", option::Arg::None,
   "  --recursive, -r              \tRecursively search subdirectories."},
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
      std::string files_str = options[createdb_opts::CreateDbOptionIndex::SourceVector].arg;
      std::stringstream ss(files_str);
      std::string file;
      while (std::getline(ss, file, ',')) {
        if (!file.empty()) cli.file_vector.push_back(file);
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
      cli.extension = options[createdb_opts::CreateDbOptionIndex::Extension].arg;
    }

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

    Logger::get_logger()->info("Creating database...");
    switch (cli.source_mode) {
      case CreateDbOptions::SourceMode::Directory:
        Logger::get_logger()->info("Source directory: {}", cli.directory_path);
        Logger::get_logger()->info("Extension: {}", cli.extension);
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

    initialize_runtime(cli.threads);

    LMDBDatabase db(cli.database_path);
    auto writer = db.get_writer();

    Logger::get_logger()->debug("Processing files...");
    const auto t0 = std::chrono::high_resolution_clock::now();

    // build pipeline
    Collector<Payload, LMDBPolicy> db_collector{cli.batch_size, writer};
    auto compute_stage = stage(dsl::thread_safe,
        [](std::string path, IEmitter<Payload>& out) {
          WriterRes result;
          result.file_path = std::move(path);
          result.success = false;
          try {
            if (gemmi::iends_with(result.file_path, ".gz")) {
              gemmi::CharArray buffer = gemmi::MaybeGzipped(result.file_path).uncompress_into_buffer();
              result.data = parse_model(buffer.data(), buffer.size());
              if (!mock_build_model_topology(result.data)) {
                Logger::get_logger()->error("Failed to build topology for file: {}", result.file_path);
              } else {
                result.success = true;
              }
            } else {
              MemoryMapped mm(result.file_path);
              if (!mm.isValid()) {
                Logger::get_logger()->critical("Error opening file: {}", result.file_path);
              } else {
                const char *data = reinterpret_cast<const char *>(mm.getData());
                size_t size = static_cast<size_t>(mm.size());
                result.data = parse_model(data, size);
                if (!mock_build_model_topology(result.data)) {
                  Logger::get_logger()->error("Failed to build topology for file: {}", result.file_path);
                } else {
                  result.success = true;
                }
              }
            }
          } catch (const std::exception &e) {
            Logger::get_logger()->error("Error processing file {}: {}", result.file_path, e.what());
          }
          out.emit(std::make_shared<const WriterRes>(std::move(result)));
        });

    Source source_variant = pick_source(cli);
    std::visit([&](auto&& src) {
      auto pipe = dsl::source_t{std::move(src)} | compute_stage | db_collector;
      dsl::run(pipe, cli.threads);
    }, source_variant);
    db_collector.finish();

    const auto t1 = std::chrono::high_resolution_clock::now();
    const auto duration = std::chrono::duration<double>(t1 - t0).count();

    Logger::get_logger()->info("Database creation completed in {:.2f} seconds!", duration);
    return 0;

  } catch (const std::exception& e) {
    Logger::get_logger()->error("Error: {}", e.what());
    return 1;
  }
}

} // namespace lahuta::cli
