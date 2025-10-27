#include <algorithm>
#include <array>
#include <atomic>
#include <cstring>
#include <iostream>
#include <memory>
#include <numeric>
#include <string>

#include "analysis/contacts/computation.hpp"
#include "analysis/contacts/provider.hpp"
#include "entities/interaction_types.hpp"
#include "logging.hpp"
#include "models/plddt.hpp"
#include "pipeline/compute/parameters.hpp"
#include "pipeline/data_requirements.hpp"
#include "pipeline/dynamic/keys.hpp"
#include "pipeline/dynamic/manager.hpp"
#include "pipeline/dynamic/sources.hpp"
#include "pipeline/dynamic/types.hpp"
#include "pipeline/model_payload.hpp"
#include "serialization/json.hpp"
#include "sinks/ndjson.hpp"
#include "spatial/fastns.hpp"
#include "topology.hpp"

// clang-format off
namespace {

using namespace lahuta;
using namespace lahuta::pipeline::dynamic;
using pipeline::DataField;
using pipeline::DataFieldSet;

constexpr std::size_t DefaultBatchSize = 512;
constexpr char OutputFile[] = "output.jsonl";

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

    logger->info("{} pipeline items: total={} processed={} skipped={} "
                 "throughput={:.2f} items/s",
                 label, report.items_total, report.items_processed,
                 report.items_skipped, throughput);
  } else {
    logger->info("{} pipeline summary: total={:.3f}s (metrics disabled)", label, report.total_seconds);
    logger->info("{} pipeline items: metrics disabled; totals unavailable", label);
  }

  logger->info("{} pipeline resources: stages={} threads_requested={} "
               "threads_used={} all_thread_safe={} run_token={}",
               label, report.stage_count, report.threads_requested,
               report.threads_used, report.all_thread_safe ? "yes" : "no",
               report.run_token);
}

class NoOpTask final : public ITask {
public:
  TaskResult run(const std::string &, TaskContext &) override { return TaskResult{true, {}}; }
};

class TopologyOnlyTask final : public ITask {
public:
  TaskResult run(const std::string &item_path, TaskContext &ctx) override {
    auto topology = ctx.get_object<const Topology>(pipeline::CTX_TOPOLOGY_KEY);
    if (!topology) {
      std::cerr << "[topology] Missing topology for '" << item_path << "'\n";
      return TaskResult{false, {}};
    }
    return TaskResult{true, {}};
  }
};

class NeighborCountTask final : public ITask {
public:
  NeighborCountTask(double cutoff, std::shared_ptr<std::atomic<std::uint64_t>> total_pairs)
      : cutoff_(cutoff), total_pairs_(std::move(total_pairs)) {}

  TaskResult run(const std::string &item_path, TaskContext &ctx) override {
    auto conformer = ctx.get_object<const RDKit::Conformer>(pipeline::CTX_CONFORMER_KEY);
    if (!conformer) {
      std::cerr << "[neighbors] Missing conformer for '" << item_path << "'\n";
      return TaskResult{false, {}};
    }

    FastNS grid(conformer->getPositions());
    if (!grid.build(cutoff_)) {
      std::cerr << "[neighbors] FastNS build failed for '" << item_path << "'\n";
      return TaskResult{false, {}};
    }

    const auto pairs = grid.self_search();
    total_pairs_->fetch_add(static_cast<std::uint64_t>(pairs.size()), std::memory_order_relaxed);
    return TaskResult{true, {}};
  }

private:
  double cutoff_;
  std::shared_ptr<std::atomic<std::uint64_t>> total_pairs_;
};

class PlddtSummaryTask final : public ITask {
public:
  explicit PlddtSummaryTask(const std::string &output_channel) : output_channel_(output_channel) {}

  TaskResult run(const std::string &item_path, TaskContext &ctx) override {
    auto payload = ctx.get_object<const pipeline::ModelPayloadSlices>(pipeline::CTX_MODEL_PAYLOAD_KEY);
    if (!payload || !payload->plddts) {
      std::cerr << "[plddt_summary] Missing pLDDT scores for '" << item_path << "'\n";
      return TaskResult{false, {}};
    }
    const auto &categories = *payload->plddts;

    constexpr std::size_t kBucketCount = static_cast<std::size_t>(pLDDTCategory::VeryLow) + 1;
    std::array<std::size_t, kBucketCount> counts{};
    for (auto cat : categories) {
      const auto idx = static_cast<std::size_t>(cat);
      if (idx < counts.size()) ++counts[idx];
    }

    const auto total = std::accumulate(counts.begin(), counts.end(), std::size_t{0});
    JsonBuilder json(256);
    json.begin_object().key("model").value(item_path).key("total").value(static_cast<unsigned int>(total));

    auto emit_pct = [&](const char *key, pLDDTCategory cat) {
      const auto idx = static_cast<std::size_t>(cat);
      const double pct = (total == 0) ? 0.0 : (static_cast<double>(counts[idx]) / static_cast<double>(total));
      json.key(key).value(pct);
    };

    emit_pct("very_low",  pLDDTCategory::VeryLow);
    emit_pct("low",       pLDDTCategory::Low);
    emit_pct("high",      pLDDTCategory::High);
    emit_pct("very_high", pLDDTCategory::VeryHigh);
    json.end_object();

    TaskResult result;
    result.ok = true;
    result.emits.push_back(Emission{output_channel_, json.str()});
    return result;
  }

  DataFieldSet data_requirements() const override { return DataFieldSet::of({DataField::Plddt}); }

private:
  std::string output_channel_;
};

int run_noop(const std::string &db_path, std::size_t threads, std::size_t batch_size) {
  StageManager manager(sources_factory::from_lmdb(db_path, std::string{}, batch_size));
  manager.set_auto_builtins(false);
  manager.add_task("noop", {}, std::make_shared<NoOpTask>(), /*thread_safe=*/true);

  const auto report = manager.run(threads);
  log_pipeline_summary("noop", report);
  return 0;
}

int run_topology(const std::string &db_path, std::size_t threads, std::size_t batch_size) {
  StageManager manager(sources_factory::from_lmdb(db_path, std::string{}, batch_size));
  manager.set_auto_builtins(true);

  auto &sys_params = manager.get_system_params();
  sys_params.is_model = true;

  auto &topo_params = manager.get_topology_params();
  topo_params.flags = TopologyComputation::All;

  manager.add_task("topology_guard", {"topology"}, std::make_shared<TopologyOnlyTask>(), true);

  const auto report = manager.run(threads);
  log_pipeline_summary("topology", report);
  return 0;
}

int run_ionic(const std::string &db_path, std::size_t threads, std::size_t batch_size) {
  StageManager manager(sources_factory::from_lmdb(db_path, std::string{}, batch_size));
  manager.set_auto_builtins(true);

  auto &sys_params = manager.get_system_params();
  sys_params.is_model = true;

  auto &topo_params = manager.get_topology_params();
  topo_params.flags = TopologyComputation::All;
  topo_params.atom_typing_method = analysis::contacts::typing_for_provider(analysis::contacts::ContactProvider::MolStar);

  pipeline::compute::ContactsParams params{};
  params.provider = analysis::contacts::ContactProvider::MolStar;
  params.type = InteractionTypeSet{InteractionType::Ionic};
  params.channel = "ionic";
  params.format = pipeline::compute::ContactsOutputFormat::Json;

  manager.add_computation(
      "ionic",
      {},
      [params]() {
        return std::make_unique<analysis::contacts::ContactsComputation>("ionic", params);
      },
      /*thread_safe=*/true);

  manager.connect_sink(params.channel, std::make_shared<NdjsonFileSink>(OutputFile));

  const auto report = manager.run(threads);
  log_pipeline_summary("ionic", report);
  Logger::get_logger()->info("[ionic] Results -> {}", OutputFile);
  return 0;
}

int run_contacts(const std::string &db_path, std::size_t threads, std::size_t batch_size) {
  StageManager manager(sources_factory::from_lmdb(db_path, std::string{}, batch_size));
  manager.set_auto_builtins(true);

  auto &sys_params = manager.get_system_params();
  sys_params.is_model = true;

  auto &topo_params = manager.get_topology_params();
  topo_params.flags = TopologyComputation::All;
  topo_params.atom_typing_method = analysis::contacts::typing_for_provider(analysis::contacts::ContactProvider::MolStar);

  pipeline::compute::ContactsParams params{};
  params.provider = analysis::contacts::ContactProvider::MolStar;
  params.type = InteractionTypeSet::all();
  params.channel = "contacts";
  params.format = pipeline::compute::ContactsOutputFormat::Json;

  manager.add_computation(
      "contacts",
      {},
      [params]() {
        return std::make_unique<analysis::contacts::ContactsComputation>("contacts", params);
      },
      /*thread_safe=*/true);

  manager.connect_sink(params.channel, std::make_shared<NdjsonFileSink>(OutputFile));

  const auto report = manager.run(threads);
  log_pipeline_summary("contacts", report);
  Logger::get_logger()->info("[contacts] Results -> {}", OutputFile);
  return 0;
}

int run_neighbors(const std::string &db_path, std::size_t threads, std::size_t batch_size) {
  constexpr double NeighborCutoff = 5.0;
  auto total_pairs = std::make_shared<std::atomic<std::uint64_t>>(0);

  StageManager manager(sources_factory::from_lmdb(db_path, std::string{}, batch_size));
  manager.set_auto_builtins(true);
  manager.get_system_params().is_model = true;

  auto task = std::make_shared<NeighborCountTask>(NeighborCutoff, total_pairs);
  manager.add_task("neighbor_counts", {"system"}, std::move(task), /*thread_safe=*/true);

  const auto report = manager.run(threads);
  log_pipeline_summary("neighbors", report);
  Logger::get_logger()->info("[neighbors] total neighbor pairs={} (cutoff={} A)", total_pairs->load(), NeighborCutoff);
  return 0;
}

int run_plddt_stats(const std::string &db_path, std::size_t threads, std::size_t batch_size) {
  constexpr char OutputChannel[] = "plddt_stats";

  StageManager manager(sources_factory::from_lmdb(db_path, std::string{}, batch_size));
  manager.set_auto_builtins(true);
  manager.get_system_params().is_model = true;

  manager.add_task("plddt_stats", {}, std::make_shared<PlddtSummaryTask>(OutputChannel), /*thread_safe=*/true);
  manager.connect_sink(OutputChannel, std::make_shared<NdjsonFileSink>(OutputFile));

  manager.compile();
  const auto report = manager.run(threads);
  log_pipeline_summary("plddt_stats", report);
  Logger::get_logger()->info("[plddt_stats] Results -> {}", OutputFile);
  return 0;
}

void print_usage(const char *prog) {
  std::cerr << "Usage: " << prog << " [--example=<name>] <database_path> <threads> [batch_size]\n";
  std::cerr << "\nAvailable examples:\n";
  std::cerr << "  noop               - No-op task (default)\n";
  std::cerr << "  topology           - Load topology only\n";
  std::cerr << "  ionic              - Compute ionic contacts (output -> output.jsonl)\n";
  std::cerr << "  contacts           - Compute all MolStar contacts (output -> output.jsonl)\n";
  std::cerr << "  neighbors          - Count neighbor pairs\n";
  std::cerr << "  plddt_stats        - Summarize pLDDT scores (output -> output.jsonl)\n";
  std::cerr << "\nOptions:\n";
  std::cerr << "  --example=<name>   - Select example (default: noop)\n";
  std::cerr << "  [batch_size]       - Batch size for LMDB source (default: 512)\n";
  std::cerr << "\nExamples:\n";
  std::cerr << "  " << prog << " --example=noop db.lmdb 4\n";
  std::cerr << "  " << prog << " --example=plddt_stats db.lmdb 4\n";
  std::cerr << "  " << prog << " db.lmdb 4              # defaults to noop\n";
}

} // namespace

int main(int argc, char **argv) {
  Logger::get_instance().set_log_level(Logger::LogLevel::Info);

  if (argc < 3) {
    print_usage(argv[0]);
    return 1;
  }

  std::string example = "noop";
  int arg_offset = 1;
  if (std::strncmp(argv[1], "--example=", 10) == 0) {
    example = argv[1] + 10;
    arg_offset = 2;
  }

  if (argc < arg_offset + 2) {
    print_usage(argv[0]);
    return 1;
  }

  const std::string db_path = argv[arg_offset];
  const std::size_t threads = std::max<std::size_t>(1, std::stoul(argv[arg_offset + 1]));

  try {
    if (example == "noop") {
      const std::size_t batch_size = (argc > arg_offset + 2) ? std::stoul(argv[arg_offset + 2]) : DefaultBatchSize;
      return run_noop(db_path, threads, batch_size);
    } else if (example == "topology") {
      const std::size_t batch_size = (argc > arg_offset + 2) ? std::stoul(argv[arg_offset + 2]) : DefaultBatchSize;
      return run_topology(db_path, threads, batch_size);
    } else if (example == "ionic") {
      const std::size_t batch_size = (argc > arg_offset + 2) ? std::stoul(argv[arg_offset + 2]) : DefaultBatchSize;
      return run_ionic(db_path, threads, batch_size);
    } else if (example == "contacts") {
      const std::size_t batch_size = (argc > arg_offset + 2) ? std::stoul(argv[arg_offset + 2]) : DefaultBatchSize;
      return run_contacts(db_path, threads, batch_size);
    } else if (example == "neighbors") {
      const std::size_t batch_size = (argc > arg_offset + 2) ? std::stoul(argv[arg_offset + 2]) : DefaultBatchSize;
      return run_neighbors(db_path, threads, batch_size);
    } else if (example == "plddt_stats") {
      const std::size_t batch_size = (argc > arg_offset + 2) ? std::stoul(argv[arg_offset + 2]) : DefaultBatchSize;
      return run_plddt_stats(db_path, threads, batch_size);
    }

    std::cerr << "Error: Unknown example '" << example << "'\n";
    print_usage(argv[0]);
    return 1;

  } catch (const std::exception &ex) {
    std::cerr << "[lmdb_examples] Error: " << ex.what() << '\n';
    return 1;
  }
}
