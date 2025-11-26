#include <filesystem>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <GraphMol/Conformer.h>

#include "analysis/contacts/computation.hpp"
#include "pipeline/dynamic/manager.hpp"
#include "pipeline/ingestion.hpp"
#include "sinks/ndjson.hpp"

// clang-format off
namespace {

using namespace lahuta;
using lahuta::sources::IDescriptor;
using namespace lahuta::pipeline::dynamic;
using namespace lahuta::pipeline::compute;

class SingleTrajectoryDescriptor final : public IDescriptor {
public:
  SingleTrajectoryDescriptor(std::string structure, std::vector<std::string> xtcs)
      : structure_(std::move(structure)), xtcs_(std::move(xtcs)) {}

  std::optional<IngestDescriptor> next() override {
    if (done_) return std::nullopt;
    done_ = true;
    IngestDescriptor desc;
    desc.id = structure_;
    desc.origin = MDRef{structure_, xtcs_};
    return desc;
  }

  void reset() override { done_ = false; }

private:
  std::string structure_;
  std::vector<std::string> xtcs_;
  bool done_ = false;
};

bool has_extension(const std::string& path, const std::vector<std::string>& extensions) {
  std::filesystem::path p(path);
  std::string ext = p.extension().string();
  std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
  for (const auto& valid_ext : extensions) {
    if (ext == valid_ext) return true;
  }
  return false;
}

} // namespace

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " <structure.pdb|.gro> <traj1.xtc> [traj2.xtc ...] [output.jsonl]\n";
    std::cerr << "  structure: PDB or GRO file (required)\n";
    std::cerr << "  trajectories: one or more XTC files (required)\n";
    std::cerr << "  output: optional output JSONL file (default: md_contact_data.jsonl)\n";
    return 1;
  }

  std::string structure_path = argv[1];

  // Validate structure file has .pdb or .gro extension
  if (!has_extension(structure_path, {".pdb", ".gro"})) {
    std::cerr << "Error: Structure file must have .pdb or .gro extension\n";
    return 1;
  }

  if (!std::filesystem::exists(structure_path)) {
    std::cerr << "Error: Structure file not found: " << structure_path << "\n";
    return 1;
  }

  // Collect XTC files and determine output file
  std::vector<std::string> xtc_files;
  std::string output_file = "md_contact_data.jsonl";

  for (int i = 2; i < argc; ++i) {
    std::string arg = argv[i];

    // Check if this is the last argument and looks like an output file
    if (i == argc - 1 && has_extension(arg, {".jsonl", ".json"})) {
      output_file = arg;
    } else {
      // Validate XTC file
      if (!has_extension(arg, {".xtc"})) {
        std::cerr << "Error: Trajectory file must have .xtc extension: " << arg << "\n";
        return 1;
      }
      if (!std::filesystem::exists(arg)) {
        std::cerr << "Error: Trajectory file not found: " << arg << "\n";
        return 1;
      }
      xtc_files.push_back(arg);
    }
  }

  if (xtc_files.empty()) {
    std::cerr << "Error: At least one XTC trajectory file is required\n";
    return 1;
  }

  std::cout << "Structure file: " << structure_path << "\n";
  std::cout << "Trajectory files (" << xtc_files.size() << "):\n";
  for (const auto& xtc : xtc_files) {
    std::cout << "  - " << xtc << "\n";
  }
  std::cout << "Output file: " << output_file << "\n\n";

  Logger::get_instance().set_log_level(Logger::LogLevel::Warn);
  auto source = std::make_unique<SingleTrajectoryDescriptor>(structure_path, xtc_files);
  StageManager manager(std::move(source));
  manager.set_auto_builtins(true);

  ContactsParams p{};
  p.provider = analysis::contacts::ContactProvider::MolStar;
  manager.add_computation(
      "contacts",
      {},
      [label = std::string("contacts"), p]() {
        return std::make_unique<analysis::contacts::ContactsComputation>(label, p);
      },
      /*thread_safe=*/true);

  auto contacts_data_sink = std::make_shared<lahuta::pipeline::dynamic::NdjsonFileSink>(output_file);
  manager.connect_sink("contacts", contacts_data_sink);

  manager.run(16);

  std::cout << "\nProcessing complete. Results written to: " << output_file << "\n";
  return 0;
}
