#include <random>
#include <string>
#include <vector>

#include "contacts/interactions.hpp"
#include "file_system.hpp"
#include "lahuta.hpp"

#include "ctpl/ctpl.h"
#include "logging.hpp"

#include "luni_props.hpp"
#include "properties/query.hpp"
#include "properties/analyzer.hpp"
#include "luni_processor.hpp"

// clang-format off

using namespace lahuta;

std::vector<std::string> get_file_paths( std::string &file_name) {
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


std::vector<std::string> random_shuffle( std::vector<std::string> &file_paths) {
  std::vector<std::string> shuffled_paths = file_paths;
  std::random_device rd;
  std::mt19937 g(rd());
  std::shuffle(shuffled_paths.begin(), shuffled_paths.end(), g);
  return shuffled_paths;
}


int main(int argc, char  *argv[]) {
  LuniProperties::initialize();

  auto logger = Logger::get_instance();
  logger.set_log_level(Logger::LogLevel::Info);
  logger.set_format(Logger::FormatStyle::Simple);

  logger.log(Logger::LogLevel::Critical, "Starting the program");

  std::string file_name = argv[1];
  std::vector<std::string> file_paths = get_file_paths(file_name);
  Logger::get_logger()->info("Number of files: {}", file_paths.size());

  auto query     = PropertyQuery<Luni>().select(PropertyKey::Names, PropertyKey::Indices);
  auto analyzer  = PropertyAnalyzer<Luni>(query);
  auto processor = FileProcessor(1, analyzer, false);


  /*Luni luni(file_paths2[0]);*/
  /*auto v = luniAnalyzer(luni);*/
  /**/
  /*const auto& names = v.get<PropertyKey::Indices>();*/
  /*const auto& indices = v.get<PropertyKey::Names>();*/
  /*Logger::get_logger()->info("Test:: File {} => #labels={}, #values={}", file_paths2[0], names.size(), indices.size());*/

  /*auto identityAnalyzer = IdentityAnalyzer<Luni>{};*/
  /*FileProcessor processorL(-1, identityAnalyzer, true);*/
  /**/
  /*processorL.process_files(file_paths);*/
  /*processorL.wait_for_completion();*/

  //! auto res_opt = processorL.get_result(file_paths2[0]);
  //! if (res_opt) {
  //!   const auto& names = res_opt->names();
  //!   const auto& indices = res_opt->indices();
  //!   Logger::get_logger()->info("Test:: File {} => #labels={}, #values={}", file_paths2[0], names.size(), indices.size());
  //! } else {
  //!   Logger::get_logger()->warn("No result found for file {}", file_paths2[0]);
  //! }

  /*auto start = std::chrono::high_resolution_clock::now();*/
  processor.process_files(file_paths);
  processor.wait_for_completion();
  /*auto end = std::chrono::high_resolution_clock::now();*/
  /*std::chrono::duration<double> elapsed_seconds = end - start;*/
  /*Logger::get_logger()->info("Elapsed time: {}", elapsed_seconds.count());*/


  return 0;

  /*start = std::chrono::high_resolution_clock::now();*/
  /*for (auto &file : file_paths) {*/
  /*  auto res_opt = processor.get_result(file);*/
  /*  if (res_opt) {*/
  /**/
  /*    const auto& names = res_opt->get<PropertyKey::Indices>();*/
  /*    const auto& indices = res_opt->get<PropertyKey::Names>();*/
  /**/
  /*    Logger::get_logger()->info("Test:: File {} => #labels={}, #values={}", file, names.size(), indices.size());*/
  /*  } else {*/
  /*    Logger::get_logger()->warn("No result found for file {}", file);*/
  /*  }*/
  /*}*/
  /*end = std::chrono::high_resolution_clock::now();*/
  /*elapsed_seconds = end - start;*/
  /*Logger::get_logger()->info("Data access time: {}", elapsed_seconds.count());*/

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
          Logger::get_logger()->info("Molecule: {}", luni.get_molecule().getNumAtoms());
        } catch (const std::exception &e) {
          Logger::get_logger()->error("Error processing file {}: {}", file_name.string(), e.what());
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
