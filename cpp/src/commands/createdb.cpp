#include "commands/createdb.hpp"
#include "cli/arg_validation.hpp"
#include "logging.hpp"

#include "db/db.hpp"
#include "io/collector.hpp"
#include "io/db_spill_policy.hpp"
#include "io/lmdb_backend.hpp"
#include "models/factory.hpp"
#include "pipeline/dsl.hpp"
#include "serialization/formats.hpp"
#include "tasks/model_db_writer.hpp"

#include <chrono>
#include <iostream>
#include <string>

// clang-format off
namespace lahuta::cli {

using namespace lahuta;
using namespace lahuta::pipeline;

struct CreateDbOptions {
  std::string directory_path;
  std::string extension = ".cif.gz";
  bool recursive = true;
  std::string database_path;
  size_t batch_size = 1000;
  int threads = 8;
};

static void initialize_factories(int num_threads) {
  lahuta::InfoPoolFactory::initialize(num_threads);
  lahuta::BondPoolFactory::initialize(num_threads);
  lahuta::AtomPoolFactory::initialize(num_threads);
}

namespace createdb_opts {
const option::Descriptor usage[] = {
  {CreateDbOptionIndex::Unknown, 0, "", "", validate::Unknown,
   "Usage: lahuta createdb [options]\n\n"
   "Create a database from alphafold2 models.\n\n"
   "Required Options:"},
  {CreateDbOptionIndex::Help, 0, "h", "help", option::Arg::None,
   "  --help, -h                   \tPrint this help message and exit."},
  {CreateDbOptionIndex::SourceDirectory, 0, "d", "directory", validate::Required,
   "  --directory, -d <path>       \tSource directory containing model files."},
  {CreateDbOptionIndex::DatabasePath, 0, "o", "output", validate::Required,
   "  --output, -o <path>          \tOutput database path."},
  {0, 0, "", "", option::Arg::None,
   "\nSource Options:"},
  {CreateDbOptionIndex::Extension, 0, "e", "extension", validate::Required,
   "  --extension, -e <ext>        \tFile extension (default: .cif.gz)."},
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

    // required options
    if (!options[createdb_opts::CreateDbOptionIndex::SourceDirectory]) {
      Logger::get_logger()->error("Directory path is required (--directory)");
      return 1;
    }
    cli.directory_path = options[createdb_opts::CreateDbOptionIndex::SourceDirectory].arg;

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
    Logger::get_logger()->info("Source directory: {}", cli.directory_path);
    Logger::get_logger()->info("Database path: {}", cli.database_path);
    Logger::get_logger()->info("Extension: {}", cli.extension);
    Logger::get_logger()->info("Recursive: {}", cli.recursive ? "Yes" : "No");
    Logger::get_logger()->info("Batch size: {}", cli.batch_size);
    Logger::get_logger()->info("Threads: {}", cli.threads);

    initialize_factories(cli.threads);

    LMDBDatabase db(cli.database_path);
    auto writer = db.get_writer();

    Logger::get_logger()->info("Processing files...");
    const auto t0 = std::chrono::high_resolution_clock::now();

    using WriterRes = tasks::ModelWriteTask::result_type;
    tasks::ModelWriteTask db_task;

    using Payload = std::shared_ptr<const WriterRes>;
    using LMDBPolicy = DBSpillPolicy<fmt::binary, Payload, LMDBBackend>;

    // build pipeline
    Collector<Payload, LMDBPolicy> db_collector{cli.batch_size, writer};
    auto pipe = dsl::directory(cli.directory_path, cli.extension, cli.recursive)
        | stage(dsl::thread_safe,
            [db_task](std::string path, IEmitter<Payload>& out) mutable {
                 auto res = db_task(std::move(path));
                 out.emit(std::make_shared<const WriterRes>(std::move(res)));
            })
        | db_collector;

    dsl::run(pipe, cli.threads);
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
