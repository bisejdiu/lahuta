#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "logging.hpp"
#include "models/dssp.hpp"
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
constexpr char OutputFile[] = "dssp_data.jsonl";

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
    logger->info("{} throughput: {:.2f} items/sec ({} items)", label, throughput, report.items_processed);
  } else {
    logger->info("{} pipeline summary: total={:.3f}s ({} items, {:.2f} items/sec)",
                 label, report.total_seconds, report.items_processed, throughput);
  }
}

char dssp_to_char(DSSPAssignment dssp) {
  switch (dssp) {
    case DSSPAssignment::Coil:             return 'C';
    case DSSPAssignment::AlphaHelix:       return 'H';
    case DSSPAssignment::Helix3_10:        return 'G';
    case DSSPAssignment::HelixPi:          return 'I';
    case DSSPAssignment::PolyProlineHelix: return 'P';
    case DSSPAssignment::Strand:           return 'E';
    case DSSPAssignment::Turn:             return 'T';
    case DSSPAssignment::Bend:             return 'S';
  }
  return '?';
}

class DsspStatsTask final : public ITask {
public:
  explicit DsspStatsTask(const std::string &output_channel) : output_channel_(output_channel) {}

  TaskResult run(const std::string &item_path, TaskContext &ctx) override {
    auto payload = ctx.model_payload();
    if (!payload) {
      throw std::runtime_error("Missing model payload for '" + item_path + "'");
    }

    if (!payload->dssp || payload->dssp->empty()) {
      throw std::runtime_error("Missing DSSP data for '" + item_path + "'");
    }

    if (!payload->plddts || payload->plddts->empty()) {
      throw std::runtime_error("Missing pLDDT data for '" + item_path + "'");
    }

    const auto &dssp_vec = *payload->dssp;
    const auto &plddt_vec = *payload->plddts;

    // Validate lengths match
    if (dssp_vec.size() != plddt_vec.size()) {
      throw std::runtime_error(
          "DSSP/pLDDT length mismatch for '" + item_path + "': " +
          "DSSP=" + std::to_string(dssp_vec.size()) + " vs " +
          "pLDDT=" + std::to_string(plddt_vec.size()));
    }

    std::string dssp_sequence;
    dssp_sequence.reserve(dssp_vec.size());
    for (const auto &dssp : dssp_vec) {
      dssp_sequence.push_back(dssp_to_char(dssp));
    }

    JsonBuilder json(512);
    json.key("model").value(item_path)
        .key("length").value(static_cast<unsigned int>(dssp_sequence.size()))
        .key("dssp_sequence").value(dssp_sequence);

    TaskResult result;
    result.ok = true;
    result.emits.push_back(Emission{output_channel_, json.str()});
    return result;
  }

  DataFieldSet data_requirements() const override {
    return DataFieldSet::of({DataField::Dssp, DataField::Plddt});
  }

private:
  std::string output_channel_;
};

int run_dssp_data(const std::string &db_path, std::size_t threads, std::size_t batch_size) {
  constexpr char OutputChannel[] = "dssp_data";

  StageManager manager(sources_factory::from_lmdb(db_path, std::string{}, batch_size, {threads + 1}));
  manager.set_auto_builtins(true);
  manager.get_system_params().is_model = true;

  auto task = std::make_shared<DsspStatsTask>(OutputChannel);
  manager.add_task("dssp_data", {}, task, /*thread_safe=*/true);
  manager.connect_sink(OutputChannel, std::make_shared<NdjsonFileSink>(OutputFile));

  const auto report = manager.run(threads);
  log_pipeline_summary("dssp_data", report);

  Logger::get_logger()->info("[dssp_data] Results written to {}", OutputFile);
  return 0;
}

void print_usage(const char *prog) {
  std::cerr << "Usage: " << prog << " <database_path> <threads> [batch_size]\n";
  std::cerr << "\nDescription:\n";
  std::cerr << "  Extracts DSSP (secondary structure) sequences from a Lahuta LMDB database.\n";
  std::cerr << "  Validates that DSSP length matches pLDDT vector length for each model.\n";
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
    return run_dssp_data(db_path, threads, batch_size);
  } catch (const std::exception &ex) {
    std::cerr << "[dssp_data] Error: " << ex.what() << '\n';
    return 1;
  }
}
