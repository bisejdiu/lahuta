#ifndef LAHUTA_PARALLEL_HPP
#define LAHUTA_PARALLEL_HPP

#include <mutex>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "contacts/interactions.hpp"
#include "lahuta.hpp"

#include "ctpl/ctpl.h"
#include "spdlog/common.h"
#include "spdlog/spdlog.h"

// clang-format off

namespace lahuta {

template <typename T>
class ThreadSafeStore {
public:
  void add_result(const std::string &file_name, const T &value) {
    std::lock_guard<std::mutex> lock(mtx_);
    results_[file_name] = value;
  }

  std::optional<T> get_result(const std::string &file_name) const {
    std::lock_guard<std::mutex> lock(mtx_);
    auto it = results_.find(file_name);
    if (it != results_.end()) {
      return it->second;
    }
    return std::nullopt;
  }

  std::unordered_map<std::string, T> get_all_results() const {
    std::lock_guard<std::mutex> lock(mtx_);
    return results_;
  }

private:
  mutable std::mutex mtx_;
  std::unordered_map<std::string, T> results_;
};

// prototype for a parallel file processor taking a callable analyzer returning type T
template <typename Analyzer, typename T>
class LuniFileProcessor {
public:
  explicit LuniFileProcessor(int concurrency, Analyzer analyzer) : analyzer_(analyzer) {
    concurrency = (concurrency <= 0) ? static_cast<int>(std::thread::hardware_concurrency()) : concurrency;
    pool_.resize(concurrency);
  }

  void process_files(const std::vector<std::string> &file_paths) {

    futures_.reserve(file_paths.size());
    for (const auto &path : file_paths) {
      futures_.emplace_back(
        // ctpl::thread_pool push returns an std::future<RetType>
        pool_.push([this, path](int /*thread_id*/) -> void {
          try {
            Luni luni(path);
            luni.build_topology();

            // invoke the analyzer on the Luni giving us the T result
            T result = analyzer_(luni);

            // keep the result keyed by file path
            // FIX: this won't work if:
            // 1. the file path is not unique
            // 2. (more seriously) number of files is very large
            threadpool_store_.add_result(path, result);

            spdlog::info("Successfully processed file: {}", path);
          } catch (const std::exception &e) {
            spdlog::error("Error processing file {}: {}", path, e.what());
          }
        })
      );
    }
  }

  /// blocks until all files have been processed
  void wait_for_completion() {
    for (auto &future : futures_) {
      future.get();
    }
    futures_.clear();
  }

  std::optional<T> get_result(const std::string &file_name) const { return threadpool_store_.get_result(file_name); }
  std::unordered_map<std::string, T> get_all_results()      const { return threadpool_store_.get_all_results(); }

private:
  ctpl::thread_pool pool_;
  Analyzer analyzer_;
  ThreadSafeStore<T> threadpool_store_;
  std::vector<std::future<void>> futures_;
};

struct LuniResult {
  std::vector<int> values;
  std::vector<std::string> labels;
};


} // namespace lahuta

#endif // LAHUTA_PARALLEL_HPP
