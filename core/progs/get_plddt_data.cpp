#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "logging.hpp"
#include "models/plddt.hpp"
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
constexpr char OutputFile[] = "plddt_data.jsonl";

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

char plddt_to_char(pLDDTCategory plddt) {
  switch (plddt) {
    case pLDDTCategory::VeryHigh: return 'E';  // Excellent (>90)
    case pLDDTCategory::High:     return 'H';  // High (70-90)
    case pLDDTCategory::Low:      return 'L';  // Low (50-70)
    case pLDDTCategory::VeryLow:  return 'P';  // Poor (<50)
  }
  return '?';
}

class PlddtDataTask final : public ITask {
public:
  explicit PlddtDataTask(const std::string &output_channel) : output_channel_(output_channel) {}

  TaskResult run(const std::string &item_path, TaskContext &ctx) override {
    auto payload = ctx.model_payload();
    if (!payload || !payload->plddts || payload->plddts->empty()) {
      Logger::get_logger()->warn("[plddt_data] Missing pLDDT data for '{}'", item_path);
      return {};
    }

    const auto &plddt_vec = *payload->plddts;

    std::string plddt_sequence;
    plddt_sequence.reserve(plddt_vec.size());
    for (const auto &plddt : plddt_vec) {
      plddt_sequence.push_back(plddt_to_char(plddt));
    }

    JsonBuilder json(512);
    json.key("model").value(item_path)
        .key("length").value(static_cast<unsigned int>(plddt_sequence.size()))
        .key("plddt_sequence").value(plddt_sequence);

    TaskResult result;
    result.ok = true;
    result.emits.push_back(Emission{output_channel_, json.str()});
    return result;
  }

  DataFieldSet data_requirements() const override {
    return DataFieldSet::of({DataField::Plddt});
  }

private:
  std::string output_channel_;
};

int run_plddt_data(const std::string &db_path, std::size_t threads, std::size_t batch_size) {
  constexpr char OutputChannel[] = "plddt_data";

  StageManager manager(sources_factory::from_lmdb(db_path, std::string{}, batch_size, {threads + 1}));
  manager.set_auto_builtins(true);
  manager.get_system_params().is_model = true;

  auto task = std::make_shared<PlddtDataTask>(OutputChannel);
  manager.add_task("plddt_data", {}, task, /*thread_safe=*/true);
  manager.connect_sink(OutputChannel, std::make_shared<NdjsonFileSink>(OutputFile));

  const auto report = manager.run(threads);
  log_pipeline_summary("plddt_data", report);

  Logger::get_logger()->info("[plddt_data] Results written to {}", OutputFile);
  return 0;
}

void print_usage(const char *prog) {
  std::cerr << "Usage: " << prog << " <database_path> <threads> [batch_size]\n";
  std::cerr << "\nDescription:\n";
  std::cerr << "  Extracts pLDDT (confidence) sequences from a Lahuta LMDB database.\n";
  std::cerr << "  Each residue is encoded with a single character representing its confidence level.\n";
  std::cerr << "\n  pLDDT character codes:\n";
  std::cerr << "    E - Excellent (VeryHigh: pLDDT > 90)\n";
  std::cerr << "    H - High (pLDDT 70-90)\n";
  std::cerr << "    L - Low (pLDDT 50-70)\n";
  std::cerr << "    P - Poor (VeryLow: pLDDT < 50)\n";
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
  Logger::get_instance().set_log_level(Logger::LogLevel::Info);

  if (argc < 3) {
    print_usage(argv[0]);
    return 1;
  }

  const std::string db_path = argv[1];
  const std::size_t threads = std::stoul(argv[2]);
  const std::size_t batch_size = (argc > 3) ? std::stoul(argv[3]) : DefaultBatchSize;

  try {
    return run_plddt_data(db_path, threads, batch_size);
  } catch (const std::exception &ex) {
    Logger::get_logger()->error("[plddt_data] Fatal error: {}", ex.what());
    return 1;
  }
}
