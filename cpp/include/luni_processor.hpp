#ifndef LAHUTA_PROCESSOR_HPP
#define LAHUTA_PROCESSOR_HPP

#include "ctpl/ctpl.h"
#include "logging.hpp"
#include <memory>
#include <string>
#include <unordered_map>

// clang-format off

namespace lahuta {

struct NoOpCallback {
  inline void operator()(const std::string&) const noexcept { }
};

/// Functor that wraps spinner updates
struct SpinnerCallback {
  indicators::MinimalProgressSpinner* spinner;

  inline void operator()(const std::string& file_path) const {
    spinner->set_postfix_text("Processed: " + file_path.substr(file_path.find_last_of('/') + 1));
    spinner->tick();
  }
};


template <typename SourceType, typename Analyzer>
class FileProcessor {
public:
  using ResultType = typename Analyzer::ResultType;

  explicit FileProcessor(int concurrency, Analyzer analyzer, bool use_spinner = false) : analyzer_(std::move(analyzer)), use_spinner_(use_spinner) {
    int concurr = (concurrency <= 0) ? static_cast<int>(std::thread::hardware_concurrency()) : concurrency;
    pool_.resize(concurr);
  }

  void process_files(const std::vector<std::string>& file_paths) {
    futures_.reserve(file_paths.size());

    if (use_spinner_) {
      hide_cursor();

      auto existing_logger = spdlog::default_logger();

      spinner_.set_max_progress(file_paths.size());
      Logger::get_instance().configure_for_spinner(&spinner_, spinner_.get_mutex());
      auto callback = SpinnerCallback{ &spinner_ };

      for (const auto& path : file_paths) {
        futures_.push_back(
          pool_.push([this, path, callback](int /*thread_id*/) {
            this->process_file(path, callback);
          })
        );
      }
      spdlog::set_default_logger(existing_logger);

    } else {

      auto callback = NoOpCallback{};

      for (const auto& path : file_paths) {
        futures_.push_back(
          pool_.push([this, path, callback](int /*thread_id*/) {
            this->process_file(path, callback);
          })
        );
      }
    }
  }

  void wait_for_completion() {
    for (auto& future : futures_) {
        future.get(); 
    }

    futures_.clear();
    if (use_spinner_) spinner_cleanup();

  }

  /// Return a pointer to the result. If there's no entry, returns nullptr.
  ResultType* get_result(const std::string& file_name) {
    std::lock_guard<std::mutex> lock(mutex_);
    auto it = results_.find(file_name);
    if (it == results_.end()) {
        return nullptr;
    }
    return it->second.get();
  }

  /// Return all results as a map of string -> pointer.
  std::unordered_map<std::string, ResultType*> get_all_results() {
    std::lock_guard<std::mutex> lock(mutex_);

    std::unordered_map<std::string, ResultType*> out;
    out.reserve(results_.size());
    for (auto & [key, value_ptr] : results_) {
        out[key] = value_ptr.get();
    }
    return out;
  }

private:
  template <typename OnTickCallback>
  void process_file(const std::string& file_path, OnTickCallback on_tick_callback) {
    try {
      SourceType source(file_path);
      source.build_topology();

      ResultType result = analyzer_(source);

      {
        std::lock_guard<std::mutex> lock(mutex_);
        results_.emplace(file_path, std::make_unique<ResultType>(std::move(result)));
      }

      on_tick_callback(file_path);

      Logger::get_logger()->info("Successfully processed file: {}", file_path);
    } catch (const std::exception& e) {
      Logger::get_logger()->error("Error processing file {}: {}", file_path, e.what());
    }
  }

  void spinner_cleanup() { // it's a bit annoying that FileProcessor has to know about spinner internals (oh well)
    spinner_.set_state("none");
    spinner_.set_postfix_text("Done! Run results:");
    spinner_.print_progress();
    show_cursor();
    Logger::get_instance().configure_for_spinner(nullptr, spinner_.get_mutex());
  }

  void hide_cursor() { std::cout << "\033[?25l"; }
  void show_cursor() { std::cout << "\033[?25h"; }

  // whether to use a spinner for progress reporting
  bool use_spinner_;
  indicators::MinimalProgressSpinner spinner_ = {"Processing files ...", 0};

  // Thread pool and futures for async processing
  ctpl::thread_pool pool_;
  std::vector<std::future<void>> futures_;

  Analyzer analyzer_;

  // Thread-safe store of file -> ResultType
  std::mutex mutex_;
  std::unordered_map<std::string, std::unique_ptr<ResultType>> results_;
};


//
// yo dawg, I heard you like templates
// Thhis automatically deduces FileProcessor<T, AnalyzerTemplate<T>> from AnalyzerTemplate<T>
// This is called a "template template parameter":
// See:
//    https://stackoverflow.com/questions/213761
//    https://en.cppreference.com/w/cpp/language/template_parameters
//                                                                          - Besian, March 2025
template <typename T, template <typename> typename AnalyzerTemplate>
FileProcessor(int, AnalyzerTemplate<T>, bool = false) -> FileProcessor<T, AnalyzerTemplate<T>>;



} // namespace lahuta

#endif // LAHUTA_PROCESSOR_HPP
