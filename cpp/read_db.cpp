#include "ctpl/ctpl.h"
#include "logging.hpp"

#include "db/db.hpp"
#include "models/factory.hpp"
#include "models/topology.hpp"

// clang-format off
using namespace lahuta;

class SimpleKeyProcessor {
public:
  explicit SimpleKeyProcessor(LMDBDatabase &database, int concurrency)
      : db_(database), processed_files_(0), is_running_(false) {
    int concurr = (concurrency <= 0) ? static_cast<int>(std::thread::hardware_concurrency()) : concurrency;
    pool_.resize(concurr);
    initialize_memory_pools();
  }

  void process_keys() {
    if (is_running_) {
      Logger::get_logger()->error("Processing already in progress");
      return;
    }
    is_running_ = true;

    futures_.clear();

    db_.for_each_key([this](const std::string &key) {
      futures_.push_back(
        pool_.push([this, key](int /*thread_id*/) {
          this->process_key(key);
        })
      );
    });

    wait_for_completion();
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
  void initialize_memory_pools() {
    InfoPoolFactory::initialize(pool_.size());
    BondPoolFactory::initialize(pool_.size());
    AtomPoolFactory::initialize(pool_.size());
  }

  void process_key(const std::string &key) {
    ModelParserResult result;
    try {
      if (!db_.get_reader().retrieve(key, result)) {
        throw std::runtime_error("Failed to retrieve model");
      }
      std::shared_ptr<RDKit::RWMol> mol = std::make_shared<RDKit::RWMol>();
      build_model_topology(mol, result, ModelTopologyMethod::CSR);
      if (!mol) {
        throw std::runtime_error("Failed to build model topology");
      }
      Logger::get_logger()->info("Processed key: {}, mol has {} atoms", key, mol->getNumAtoms());
    } catch (const std::exception &e) {
      Logger::get_logger()->error("Error processing key {}: {}", key, e.what());
      increment_progress();
      return;
    }
    increment_progress();
  }

  void increment_progress() {
    static constexpr size_t progress_interval = 500;
    size_t count = ++processed_files_;
    if (count % progress_interval == 0) {
      std::cout << "Processed " << count << " keys." << std::endl;
    }
  }

  LMDBDatabase &db_;
  ctpl::thread_pool pool_;
  std::vector<std::future<void>> futures_;
  std::atomic<size_t> processed_files_;
  std::atomic<bool> is_running_;
};

int main(int argc, char *argv[]) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " <db_path>" << "n_jobs" << std::endl;
    return 1;
  }
  std::string db_path = argv[1];
  int n_jobs = std::stoi(argv[2]);

  LMDBDatabase db(db_path);
  auto logger = Logger::get_instance();
  logger.set_log_level(Logger::LogLevel::Info);
  logger.set_format(Logger::FormatStyle::Simple);

  SimpleKeyProcessor processor(db, n_jobs);
  processor.process_keys();
  processor.shutdown();

  return 0;
}
