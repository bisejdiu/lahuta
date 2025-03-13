#include <optional>
#include <string>
#include <vector>

#include "contacts/interactions.hpp"
#include "file_system.hpp"
#include "lahuta.hpp"

#include "ctpl/ctpl.h"
#include "logging.hpp"
#include "parallel.hpp"

using namespace lahuta;

std::vector<std::string> get_file_paths(const std::string &file_name) {
  std::ifstream file_list(file_name);
  if (!file_list.is_open()) {
    throw std::runtime_error("Failed to open file: " + file_name);
  }

  std::vector<std::string> file_paths;
  std::string line;
  while (std::getline(file_list, line)) {
    if (!line.empty()) {
      file_paths.push_back(line);
    }
  }
  file_list.close();
  return file_paths;
}

int main(int argc, char const *argv[]) {

  Logger::get_instance().set_log_level(Logger::LogLevel::Info);
  Logger::get_instance().set_format(Logger::FormatStyle::Detailed);

  auto analyzer = [](const Luni &luni) -> LuniResult {
    LuniResult r;
    r.values = luni.indices();
    r.labels = luni.names();
    return r;
  };

  LuniFileProcessor<decltype(analyzer), LuniResult> processor(
      /* concurrency */ 4,
      analyzer);

  std::string file_name = argv[1];
  std::vector<std::string> file_paths = get_file_paths(file_name);
  spdlog::info("Number of files: {}", file_paths.size());

  processor.process_files(file_paths);
  processor.wait_for_completion();

  for (const auto &file : file_paths) {
    auto res_opt = processor.get_result(file);
    if (res_opt) {
      spdlog::info("Test:: File {} => #labels={}, #values={}", file, res_opt->labels.size(), res_opt->values.size());
    } else {
      spdlog::warn("No result found for file {}", file);
    }
  }

  return 0;
}

int dir_based_main() {
  try {
    DirectoryHandler dir_handler("/Users/bsejdiu/data/mini_pdb", ".cif.gz", true);

    int num_threads = std::thread::hardware_concurrency();
    ctpl::thread_pool pool(num_threads);

    std::vector<std::future<void>> futures;

    for (const auto &file : dir_handler) {
      auto file_name = file.get_path();
      futures.emplace_back(pool.push([file_name](int id) {
        try {
          Luni luni(file_name);
          spdlog::info("Molecule: {}", luni.get_molecule().getNumAtoms());
        } catch (const std::exception &e) {
          std::cerr << "Error processing file " << file_name << ": " << e.what() << std::endl;
        }
      }));
    }

    for (auto &fut : futures) {
      fut.get();
    }

  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }
  return 0;
}
