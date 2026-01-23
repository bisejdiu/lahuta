#ifndef LAHUTA_PROCESSOR_HPP
#define LAHUTA_PROCESSOR_HPP

#include <memory>
#include <string>
#include <unordered_map>

#include <ctpl/ctpl.h>

#include "logging/logging.hpp"
#include "models/factory.hpp"

// clang-format off
namespace lahuta {

constexpr size_t MAX_MEM_THRESHOLD = 1024ull * 1024 * 10000; // 10 GB

template <typename T>
class ResultStore {
public:
  // FIX: in Python, alias int to Bytes
  explicit ResultStore(size_t mem_threshold = MAX_MEM_THRESHOLD) { (void)mem_threshold; }

  // 'value' is constructed in-place as a std::unique_ptr<T>
  void add_result(const std::string &file_name, T &&value) {
      std::lock_guard<std::mutex> lock(mtx_);

      auto value_ = std::make_unique<T>(std::move(value));

      results_.emplace(
          std::piecewise_construct,
          std::forward_as_tuple(file_name),
          std::forward_as_tuple(std::move(value_)));
  }

  // Returns a pointer to the stored T, or nullptr if not present.
  // The caller does *not* own the returned pointer.
  T *get_result(const std::string &file_name) {
      std::lock_guard<std::mutex> lock(mtx_);
      auto it = results_.find(file_name);
      if (it == results_.end()) {
          return nullptr;
      }
      return it->second.get(); // raw pointer to T
  }

  // Return references to *all* stored results
  std::unordered_map<std::string, T *> get_all_results() {
      std::lock_guard<std::mutex> lock(mtx_);
      std::unordered_map<std::string, T *> out;
      out.reserve(results_.size());
      for (auto &kv : results_) {
          out.emplace(kv.first, kv.second.get());
      }
      return out;
  }

  void clear() {
      std::lock_guard<std::mutex> lock(mtx_);
      results_.clear();
  }

private:
  mutable std::mutex mtx_;
  std::unordered_map<std::string, std::unique_ptr<T>> results_;
};

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
class FileProcessor1 {
public:
  using ResultType = typename Analyzer::ResultType;

  explicit FileProcessor1(int concurrency, Analyzer analyzer, bool use_spinner = false) : analyzer_(std::move(analyzer)), use_spinner_(use_spinner) {
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

    pool_.clear_queue();
    pool_.stop(true);
  }

  // Return a pointer to the result. If there's no entry, returns nullptr.
  ResultType* get_result(const std::string& file_name) {
    return threadpool_store_.get_result(file_name);
  }

  // Return all results as a map of string -> pointer.
  std::unordered_map<std::string, ResultType*> get_all_results() {
    return threadpool_store_.get_all_results();
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
        threadpool_store_.add_result(file_path, std::move(result));
      }

      on_tick_callback(file_path);

      Logger::get_logger()->info("Successfully processed file: {}", file_path);
    } catch (const std::exception& e) {
      Logger::get_logger()->error("Error processing file {}: {}", file_path, e.what());
    }
  }

  void spinner_cleanup() { // it's a bit annoying that FileProcessor1 has to know about spinner internals (oh well)
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
  ResultStore<ResultType> threadpool_store_;
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

    InfoPoolFactory::initialize(pool_.size());
    BondPoolFactory::initialize(pool_.size());
    AtomPoolFactory::initialize(pool_.size());

    if (use_spinner_) {
      hide_cursor();

      auto current_format = Logger::get_instance().get_format_style();

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
      Logger::get_instance().set_format(current_format);

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

    pool_.clear_queue();
    pool_.stop(true);
  }

  // Return a pointer to the result. If there's no entry, returns nullptr.
  ResultType* get_result(const std::string& file_name) {
    return threadpool_store_.get_result(file_name);
  }

  // Return all results as a map of string -> pointer.
  std::unordered_map<std::string, ResultType*> get_all_results() {
    return threadpool_store_.get_all_results();
  }

private:
  template <typename OnTickCallback>
  void process_file(const std::string& file_path, OnTickCallback on_tick_callback) {
    try {
      auto source = std::make_unique<SourceType>(file_path, typename SourceType::ModelFileTag{});
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
  ResultStore<ResultType> threadpool_store_;
};

//
// Thhis automatically deduces FileProcessor<T, AnalyzerTemplate<T>> from AnalyzerTemplate<T>
// This is called a "template template parameter":
// See:
//    https://stackoverflow.com/questions/213761
//    https://en.cppreference.com/w/cpp/language/template_parameters
template <typename T, template <typename> typename AnalyzerTemplate>
FileProcessor(int, AnalyzerTemplate<T>, bool = false) -> FileProcessor<T, AnalyzerTemplate<T>>;

template <typename T, template <typename> typename AnalyzerTemplate>
FileProcessor1(int, AnalyzerTemplate<T>, bool = false) -> FileProcessor1<T, AnalyzerTemplate<T>>;

} // namespace lahuta

#endif // LAHUTA_PROCESSOR_HPP
