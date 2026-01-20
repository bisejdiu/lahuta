#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <Geometry/point.h>

#include "logging.hpp"
#include "pipeline/data_requirements.hpp"
#include "pipeline/dynamic/manager.hpp"
#include "pipeline/dynamic/sources.hpp"
#include "pipeline/dynamic/types.hpp"

// clang-format off
namespace {

using namespace lahuta;
using namespace lahuta::pipeline::dynamic;
using pipeline::DataField;
using pipeline::DataFieldSet;

constexpr std::size_t DefaultBatchSize = 512;
const std::string OutputDir = "/Users/bsejdiu/projects/lahuta_dev/lahuta/position_data";

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

void write_positions_as_float(const std::vector<RDGeom::Point3D>& pts, const std::string& filename) {
    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs) throw std::runtime_error("Failed to open file for writing: " + filename);

    // Write as flat array of floats to ensure proper memory layout
    std::vector<float> data;
    data.reserve(pts.size() * 3);
    for (const auto& pt : pts) {
        data.push_back(static_cast<float>(pt.x));
        data.push_back(static_cast<float>(pt.y));
        data.push_back(static_cast<float>(pt.z));
    }

    ofs.write(reinterpret_cast<const char*>(data.data()),
              static_cast<std::streamsize>(data.size() * sizeof(float)));
}

std::string sanitize_filename(const std::string& model_path) {
  auto last_slash = model_path.find_last_of("/\\");
  if (last_slash != std::string::npos) {
    return model_path.substr(last_slash + 1);
  }
  return model_path;
}

class PositionDataTask final : public ITask {
public:
  explicit PositionDataTask(const std::string& output_dir) : output_dir_(output_dir) {}

  TaskResult run(const std::string &item_path, TaskContext &ctx) override {
    auto payload = ctx.model_payload();
    if (!payload) {
      throw std::runtime_error("Missing model payload for '" + item_path + "'");
    }

    if (!payload->positions || payload->positions->empty()) {
      throw std::runtime_error("Missing position data for '" + item_path + "'");
    }

    const auto& positions = *payload->positions;

    std::string safe_name = sanitize_filename(item_path);
    std::string output_path = output_dir_ + "/" + safe_name + ".bin";

    try {
      write_positions_as_float(positions, output_path);

      auto logger = Logger::get_logger();
      logger->debug("[position_data] Wrote {} atoms to {}", positions.size(), output_path);
    } catch (const std::exception& e) {
      throw std::runtime_error("Failed to write position data for '" + item_path + "': " + e.what());
    }

    TaskResult result;
    result.ok = true;
    return result;
  }

  DataFieldSet data_requirements() const override {
    return DataFieldSet::of({DataField::Positions});
  }

private:
  std::string output_dir_;
};

int run_position_data(const std::string &db_path, std::size_t threads, std::size_t batch_size) {
  std::filesystem::create_directories(OutputDir);
  auto logger = Logger::get_logger();
  logger->info("[position_data] Output directory: {}", OutputDir);
  logger->info("[position_data] Data type: float32 (converted from Point3D)");

  StageManager manager(sources_factory::from_lmdb(db_path, std::string{}, batch_size, {threads + 1}));
  manager.set_auto_builtins(true);
  manager.get_system_params().is_model = true;

  auto task = std::make_shared<PositionDataTask>(OutputDir);
  manager.add_task("position_data", {}, task, /*thread_safe=*/true);

  const auto report = manager.run(threads);
  log_pipeline_summary("position_data", report);

  logger->info("[position_data] Binary files written to {}/", OutputDir);
  return 0;
}

void print_usage(const char *prog) {
  std::cerr << "Usage: " << prog << " <database_path> <threads> [batch_size]\n";
  std::cerr << "\nDescription:\n";
  std::cerr << "  Extracts 3D atomic coordinates from a Lahuta LMDB database.\n";
  std::cerr << "  Each model's coordinates are saved as a NumPy-loadable binary file.\n";
  std::cerr << "\nArguments:\n";
  std::cerr << "  <database_path>  - Path to LMDB database directory\n";
  std::cerr << "  <threads>        - Number of worker threads\n";
  std::cerr << "  [batch_size]     - Batch size for LMDB source (default: 512)\n";
  std::cerr << "\nOutput:\n";
  std::cerr << "  Files are written to: " << OutputDir << "/\n";
  std::cerr << "  Format: One binary file per model, structured as Nx3 arrays (float32)\n";
  std::cerr << "  Can be loaded in Python with: np.fromfile(path, dtype=np.float32).reshape(-1, 3)\n";
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
    return run_position_data(db_path, threads, batch_size);
  } catch (const std::exception &ex) {
    std::cerr << "[position_data] Error: " << ex.what() << '\n';
    return 1;
  }
}
