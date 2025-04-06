#include "ctpl/ctpl.h"
#include "file_system.hpp"
#include "logging.hpp"

#include "models/factory.hpp"
#include "models/topology.hpp"
#include "db/db.hpp"

// clang-format off
using namespace lahuta;

class SimpleFileProcessor {
public:
  explicit SimpleFileProcessor(LMDBReader &manager, int concurrency)
      : db_reader_(manager), processed_files_(0), is_running_(false), samples_(0), total_time_(0) {

    int concurr = (concurrency <= 0) ? static_cast<int>(std::thread::hardware_concurrency()) : concurrency;

    pool_.resize(concurr);
    initialize_memory_pools();
  }

  void process_files(DirectoryHandler &dir_handler, size_t batch_size = 10'000) {
    if (is_running_) {
      Logger::get_logger()->error("Processing already in progress");
      return;
    }
    is_running_ = true;

    while (is_running_) {
      auto maybe_chunk = dir_handler.next_chunk(batch_size);
      if (!maybe_chunk) break;

      const FileChunk &chunk = *maybe_chunk;

      futures_.clear();
      futures_.reserve(chunk.size());

      for (const auto &file_path : chunk) {
        futures_.push_back(
          pool_.push([this, file_path](int /*thread_id*/) {
            this->process_file(file_path);
          })
        );
      }
    }

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

  void print_benchmark() {
    if (samples_ == 0) return;
    auto avg_time = total_time_ / samples_;
    std::cout << "Average time per file: " << avg_time << " microseconds" << std::endl;
  }

private:
  void initialize_memory_pools() {
    InfoPoolFactory::initialize(pool_.size());
    BondPoolFactory::initialize(pool_.size());
    AtomPoolFactory::initialize(pool_.size());
  }

  void process_file(const std::string &file_path) {
    auto start = std::chrono::high_resolution_clock::now();
    ModelParserResult result;
    try {
      if (!db_reader_.retrieve(file_path, result)) {
        throw std::runtime_error("Failed to retrieve model");
      }
      std::shared_ptr<RDKit::RWMol> mol = std::make_shared<RDKit::RWMol>();
      build_model_topology(mol, result, ModelTopologyMethod::CSR);
      if (!mol) {
        throw std::runtime_error("Failed to build model topology");
      }
    } catch (const std::exception &e) {
      Logger::get_logger()->error("Error parsing file {}: {}", file_path, e.what());
      record_timing(start);
      increment_progress();
      return;
    }
    record_timing(start);
    increment_progress();
  }

  void record_timing(const std::chrono::high_resolution_clock::time_point &start) {
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    total_time_.fetch_add(duration, std::memory_order_relaxed);
    samples_.fetch_add(1, std::memory_order_relaxed);
  }

  void increment_progress() {
    static constexpr size_t progress_interval = 500;
    size_t count = ++processed_files_;
    if (count % progress_interval == 0) {
      std::cout << "Processed " << count << " files." << std::endl;
    }
  }

  LMDBReader &db_reader_;
  ctpl::thread_pool pool_;
  std::vector<std::future<void>> futures_;
  std::atomic<size_t> processed_files_;
  std::atomic<bool> is_running_;

  // benchmark stuff
  std::atomic<size_t> samples_;
  std::atomic<size_t> total_time_;
};

int main(int argc, char  *argv[]) {

  std::string data_path = argv[1];
  std::string db_path   = argv[2];

  LMDBDatabase db(db_path);
  auto reader = db.get_reader();

  auto logger = Logger::get_instance();
  logger.set_log_level(Logger::LogLevel::Info);
  logger.set_format(Logger::FormatStyle::Simple);

  DirectoryHandler dir_handler(data_path, ".cif.gz", true);
  auto processor = SimpleFileProcessor(reader, 1);
  processor.process_files(dir_handler);
  processor.shutdown();
  processor.print_benchmark();
}

