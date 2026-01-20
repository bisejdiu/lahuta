#include <cstdint>
#include <iostream>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

#include "logging.hpp"
#include "pipeline/data_requirements.hpp"
#include "pipeline/dynamic/manager.hpp"
#include "pipeline/dynamic/sources.hpp"
#include "pipeline/dynamic/types.hpp"
#include "serialization/json.hpp"
#include "sinks/ndjson.hpp"

// clang-format off
namespace {

using namespace lahuta;
using namespace lahuta::pipeline::dynamic;
using pipeline::DataField;
using pipeline::DataFieldSet;

constexpr std::size_t DefaultBatchSize = 512;
constexpr char OutputFile[] = "organism_data.jsonl";

void log_pipeline_summary(const char* label, const StageManager::RunReport& report) {
  auto logger = Logger::get_logger();
  const double throughput = (report.total_seconds > 0.0 && report.items_processed > 0)
      ? static_cast<double>(report.items_processed) / report.total_seconds : 0.0;

  if (report.metrics_enabled) {
    logger->info("{} pipeline summary: total={:.3f}s cpu={:.3f}s io={:.3f}s "
                 "(ingest={:.3f}s prepare={:.3f}s flush={:.3f}s) "
                 "setup={:.3f}s compute={:.3f}s",
                 label, report.total_seconds, report.cpu_seconds,
                 report.io_seconds, report.ingest_seconds, report.prepare_seconds,
                 report.flush_seconds, report.setup_seconds, report.compute_seconds);
    logger->info("{} throughput: {:.2f} items/sec ({} items)",
                 label, throughput, report.items_processed);
  } else {
    logger->info("{} pipeline summary: total={:.3f}s ({} items, {:.2f} items/sec)",
                 label, report.total_seconds, report.items_processed, throughput);
  }
}

class OrganismStatsTask final : public ITask {
public:
  struct Snapshot {
    std::uint64_t total = 0;
    std::uint64_t missing = 0;
    std::unordered_map<std::string, std::uint64_t> organism_counts;
  };

  explicit OrganismStatsTask(const std::string &output_channel) : output_channel_(output_channel) {}

  TaskResult run(const std::string &item_path, TaskContext &ctx) override {
    auto payload = ctx.model_payload();
    if (!payload || !payload->metadata) {
      record_observation({}, {});
      return {};
    }

    const auto &organism = payload->metadata->organism_scientific;
    const auto &taxonomy_id = payload->metadata->ncbi_taxonomy_id;
    record_observation(organism, taxonomy_id);

    JsonBuilder json(256);
    json.key("model").value(item_path)
        .key("organism").value(organism.empty() ? "Unknown" : organism)
        .key("taxonomy_id").value(taxonomy_id.empty() ? "N/A" : taxonomy_id);

    TaskResult result;
    result.ok = true;
    result.emits.push_back(Emission{output_channel_, json.str()});
    return result;
  }

  DataFieldSet data_requirements() const override { return DataFieldSet::of({DataField::Metadata}); }

  Snapshot snapshot() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return Snapshot{total_, missing_, organism_counts_};
  }

private:
  void record_observation(const std::string &organism, const std::string &taxonomy_id) {
    std::lock_guard<std::mutex> lock(mutex_);
    ++total_;
    if (organism.empty()) {
      ++missing_;
      return;
    }

    // Key with organism name and taxonomy ID
    std::string organism_key = organism;
    if (!taxonomy_id.empty()) {
      organism_key += " (TaxID: " + taxonomy_id + ")";
    }

    ++organism_counts_[organism_key];
  }

  std::string output_channel_;
  mutable std::mutex mutex_;
  std::uint64_t total_ = 0;
  std::uint64_t missing_ = 0;
  std::unordered_map<std::string, std::uint64_t> organism_counts_;
};

int run_organism_data(const std::string &db_path, std::size_t threads, std::size_t batch_size) {
  constexpr char OutputChannel[] = "organism_data";

  StageManager manager(sources_factory::from_lmdb(db_path, std::string{}, batch_size, {threads + 1}));
  manager.set_auto_builtins(true);
  manager.get_system_params().is_model = true;

  auto task = std::make_shared<OrganismStatsTask>(OutputChannel);
  manager.add_task("organism_data", {}, task, /*thread_safe=*/true);
  manager.connect_sink(OutputChannel, std::make_shared<NdjsonFileSink>(OutputFile));

  const auto report = manager.run(threads);
  log_pipeline_summary("organism_data", report);

  const auto stats = task->snapshot();
  const auto logger = Logger::get_logger();
  const auto with_species = stats.total >= stats.missing ? stats.total - stats.missing : 0;
  logger->info("[organism_stats] models={} with_species={} missing={} unique_species={}",
               stats.total, with_species, stats.missing, stats.organism_counts.size());
  logger->info("[organism_stats] Results written to {}", OutputFile);
  return 0;
}

void print_usage(const char *prog) {
  std::cerr << "Usage: " << prog << " <database_path> <threads> [batch_size]\n";
  std::cerr << "\nDescription:\n";
  std::cerr << "  Analyzes organism/species distribution in a Lahuta LMDB database.\n";
  std::cerr << "  Reports total models, those with species annotations, and top 500 species.\n";
  std::cerr << "\nArguments:\n";
  std::cerr << "  <database_path>  - Path to LMDB database directory\n";
  std::cerr << "  <threads>        - Number of worker threads\n";
  std::cerr << "  [batch_size]     - Batch size for LMDB source (default: 512)\n";
  std::cerr << "\nExample:\n";
  std::cerr << "  " << prog << " swissprot_db 16\n";
  std::cerr << "  " << prog << " swissprot_db 16 1024\n";
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 3) {
    print_usage(argv[0]);
    return 1;
  }

  const std::string db_path = argv[1];
  const std::size_t threads = std::stoul(argv[2]);
  const std::size_t batch_size = (argc > 3) ? std::stoul(argv[3]) : DefaultBatchSize;

  return run_organism_data(db_path, threads, batch_size);
}
