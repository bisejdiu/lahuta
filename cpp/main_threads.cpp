#include "contacts/interactions.hpp"
#include "file_system.hpp"
#include "lahuta.hpp"

#include "ctpl/ctpl.h"
#include "spdlog/common.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"

using namespace lahuta;

void set_logger_pattern(spdlog::level::level_enum level) {
  if (level == spdlog::level::debug) {
    spdlog::set_pattern("[%T] [%^%l%$] [thread %t] %v");
  } else {
    spdlog::set_pattern("[%^%l%$] %v");
  }
}

int process_file(const std::string &file_path) {
  try {
    Luni luni(file_path);
  } catch (const std::exception &e) {
    std::cerr << "Error processing file " << file_path << ": " << e.what() << std::endl;
    return 1;
  }
  return 0;
}

std::vector<std::string> get_file_paths(const std::string &file_name) {
  std::ifstream file_list(file_name);
  if (!file_list.is_open()) {
    throw std::runtime_error("Failed to open FILES.txt");
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

int _main(int argc, char const *argv[]) {
  auto logger = spdlog::stdout_color_mt("console");

  spdlog::level::level_enum log_level = spdlog::level::info;
  spdlog::set_level(log_level);
  set_logger_pattern(log_level);

  try {
    std::string file_name = argv[1];
    /*std::string file_name = "/Users/bsejdiu/projects/lahuta/cpp/build/F1000.txt";*/
    std::vector<std::string> file_paths = get_file_paths(file_name);

    int num_threads = std::thread::hardware_concurrency();
    if (num_threads == 0) num_threads = 1;
    ctpl::thread_pool pool(num_threads);

    for (const auto &file_name : file_paths) {
      pool.push([file_name](int id) {
        try {
          Luni luni(file_name);
        } catch (const std::exception &e) {
          spdlog::error("Error processing file {}: {}", file_name, e.what());
        }
      });
    }

  } catch (const std::exception &e) {
    spdlog::critical("Error: {}", e.what());
  }
  return 0;
}

int main(int argc, char const *argv[]) {

  auto logger = spdlog::stdout_color_mt("console");

  spdlog::level::level_enum log_level = spdlog::level::info;
  spdlog::set_level(log_level);
  set_logger_pattern(log_level);

  try {

    std::string file_name = argv[1];
    /*std::string file_name = "/Users/bsejdiu/projects/lahuta/cpp/build/F1000.txt";*/
    std::vector<std::string> file_paths = get_file_paths(file_name);

    int num_threads = std::thread::hardware_concurrency();
    if (num_threads == 0) num_threads = 1;
    std::cout << "num_threads: " << num_threads << std::endl;
    ctpl::thread_pool pool(num_threads);

    std::vector<std::future<void>> futures;
    for (const auto &file_name : file_paths) {
      futures.emplace_back(pool.push([file_name](int id) {
        try {
          Luni luni(file_name);
          /*Luni* luni = new Luni(file_name);*/
        } catch (const std::exception &e) {
          spdlog::error("Error processing file {}: {}", file_name, e.what());
        }
      }));
    }

    for (auto &fut : futures) {
      fut.get();
    }

  } catch (const std::exception &e) {
    spdlog::critical("Error: {}", e.what());
  }
  std::_Exit(EXIT_SUCCESS);
}

int file_based_main() {
  try {
    DirectoryHandler dir_handler("/Users/bsejdiu/data/PDB_ARCHIVE", ".cif.gz", true);

    int num_threads = std::thread::hardware_concurrency();
    ctpl::thread_pool pool(num_threads);

    std::vector<std::future<void>> futures;

    for (const auto &file : dir_handler) {
      std::string file_name = file.get_absolute_path();
      futures.emplace_back(pool.push([file_name](int id) {
        try {
          Luni luni(file_name);
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
