#include <optional>
#include <random>
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

std::vector<std::string> random_shuffle(const std::vector<std::string> &file_paths) {
  std::vector<std::string> shuffled_paths = file_paths;
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(shuffled_paths.begin(), shuffled_paths.end(), g);
  return shuffled_paths;
}

int main(int argc, char const *argv[]) {

  Logger::get_instance().set_log_level(Logger::LogLevel::Warn);
  Logger::get_instance().set_format(Logger::FormatStyle::Detailed);

  auto analyzer = [](const Luni &luni) -> LuniResult {
    LuniResult r;
    r.values = luni.indices();
    r.labels = luni.names();
    return r;
  };

  std::string file_name2 = argv[1];
  std::vector<std::string> file_paths2 = get_file_paths(file_name2);
  spdlog::info("Number of files: {}", file_paths2.size());

  LuniFileProcessorMinimal<decltype(analyzer), LuniResult> processor(/* concurrency */-1, analyzer);

  processor.process_files(file_paths2);
  processor.wait_for_completion();

  // Print final stats
  for (const auto &file : file_paths2) {
    auto res_opt = processor.get_result(file);
    if (res_opt) {
      spdlog::info("Test:: File {} => #labels={}, #values={}",
        file, res_opt->labels.size(), res_opt->values.size());
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
