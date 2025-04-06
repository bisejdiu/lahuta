#include <string>
#include <vector>

#include "file_system.hpp"
#include "gemmi/gz.hpp"

#include "ctpl/ctpl.h"
#include "logging.hpp"

#include "mmap/MemoryMapped.h"
#include "models/factory.hpp"
#include "db/db.hpp"
#include "models/topology.hpp"

// clang-format off
using namespace lahuta;
using gemmi::CharArray, gemmi::MaybeGzipped, gemmi::iends_with;

class CreateLahutaDB {
public:
  explicit CreateLahutaDB(LMDBWriter &manager, int concurrency)
      : db_writer_(manager), processed_files_(0), is_running_(false) {

    int concurr = (concurrency <= 0) ? static_cast<int>(std::thread::hardware_concurrency()) : concurrency;
    initialize(concurr);
  }

  void process_files(DirectoryHandler &dir_handler, size_t batch_size = 10'000) {
    if (is_running_) {
      Logger::get_logger()->error("Processing already in progress");
      return;
    }
    is_running_ = true;

    try {
      while (is_running_) {
        auto maybe_chunk = dir_handler.next_chunk(batch_size);
        if (!maybe_chunk) {
          Logger::get_logger()->info("No more files to process");
          break;
        }

        const FileChunk &chunk = *maybe_chunk;

        db_writer_.begin_txn();

        futures_.clear();
        futures_.reserve(chunk.size());

        for (const auto &file_path : chunk) {

          futures_.push_back(pool_.push([this, file_path](int /*thread_id*/) {
            try {
              ModelParserResult result;
              if (!process_file(file_path, result)) {
                Logger::get_logger()->error("Failed to process file: {}", file_path);
                return;
              }

              {
                std::lock_guard<std::mutex> lock(db_mutex_);
                db_writer_.store(file_path, result, false);
              }

            } catch (const std::exception &e) {
              Logger::get_logger()->error("Error processing file {}: {}", file_path, e.what());
            }

            increment_progress();
          }));
        }

        wait_for_completion(); // wait for all futures for this batch to finish

        db_writer_.commit_txn(); // commit the transaction
        Logger::get_logger()->info("Committed batch of {} items", chunk.size());
      }
    } catch (const std::exception &e) {
      Logger::get_logger()->critical("Error in process_files: {}", e.what());

      try { // maybe use RAII for this
        db_writer_.abort_txn();
      } catch (...) {
        Logger::get_logger()->error("Failed to abort transaction after error");
      }
    }

    Logger::get_logger()->info("All file processing complete");
    is_running_ = false;
  }

  void stop() {
    is_running_ = false;
    wait_for_completion();
  }

  void wait_for_completion() {
    for (auto &future : futures_) {
      future.get();
    }
    futures_.clear();
    pool_.clear_queue();
  }

  void shutdown() {
    stop();
    pool_.stop(true);
  }

private:
  bool process_file(const std::string &file_path, ModelParserResult &result) {

    if (iends_with(file_path, ".gz")) {
      CharArray buffer = MaybeGzipped(file_path).uncompress_into_buffer();
      result = parse_model(buffer.data(), buffer.size());
      if (!mock_build_model_topology(result)) {
        Logger::get_logger()->error("Failed to build topology for file: {}", file_path);
        increment_progress();
        return false;
      }
      return true;
    }

    MemoryMapped mm(file_path);
    if (!mm.isValid()) {
      Logger::get_logger()->critical("Error opening file: {}", file_path);
      increment_progress();
      return false;
    }
    const char *data = reinterpret_cast<const char *>(mm.getData());
    size_t size      = static_cast<size_t>(mm.size());
    result = parse_model(data, size);
    if (!mock_build_model_topology(result)) {
      Logger::get_logger()->error("Failed to build topology for file: {}", file_path);
      increment_progress();
      return false;
    }
    return true;
  }

  void increment_progress() {
    static constexpr size_t progress_interval = 500;
    size_t count = ++processed_files_;
    if (count % progress_interval == 0) {
      Logger::get_logger()->info("Processed {} files.", count);
    }
  }

  void initialize(int c = -1) {
    pool_.resize(c);
    InfoPoolFactory::initialize(pool_.size());
    BondPoolFactory::initialize(pool_.size());
    AtomPoolFactory::initialize(pool_.size());
  }

  LMDBWriter &db_writer_;
  ctpl::thread_pool pool_;
  std::vector<std::future<void>> futures_;
  std::atomic<size_t> processed_files_;
  std::atomic<bool> is_running_;
  std::mutex db_mutex_; // to protect concurrent writes to the database
};

// clang-format off
int main(int argc, char  *argv[]) {

  std::string data_path = argv[1];
  std::string db_path = argv[2];

  LMDBDatabase db(db_path);
  auto writer = db.get_writer();

  auto logger = Logger::get_instance();
  logger.set_log_level(Logger::LogLevel::Info);
  logger.set_format(Logger::FormatStyle::Simple);

  logger.log(Logger::LogLevel::Critical, "Starting the program");

  DirectoryHandler dir_handler(data_path, ".cif.gz", true);
  auto processor = CreateLahutaDB(writer, -1);
  processor.process_files(dir_handler);
  processor.shutdown();
}
