#ifndef LAHUTA_PIPELINE_SOURCES_DIRECTORY_SOURCE_HPP
#define LAHUTA_PIPELINE_SOURCES_DIRECTORY_SOURCE_HPP

#include "logging.hpp"
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

// clang-format off
namespace lahuta::sources {
namespace fs = std::filesystem;

// batch-oriented directory iterator.
class DirectorySource {
public:
  using value_type = std::string;

  explicit DirectorySource(
      std::string_view directory_path, std::string_view extension = "", bool recursive = true, std::size_t batch_size = 1024)
      : directory_path_(directory_path), extension_(extension), recursive_(recursive), batch_size_(batch_size) {

    if (!fs::exists(directory_path_) || !fs::is_directory(directory_path_)) {
      throw fs::filesystem_error(
          "Directory does not exist",
          std::make_error_code(std::errc::no_such_file_or_directory));
    }

    // correct iterator
    if (recursive_) recursive_it_.emplace(directory_path_, fs::directory_options::skip_permission_denied);
    else                  dir_it_.emplace(directory_path_, fs::directory_options::skip_permission_denied);

    Logger::get_logger()->info(
        "DirectorySource: scanning directory '{}'{}; batch_size={} (recursive={})",
        directory_path_,
        extension_.empty() ? "" : fmt::format(", filter='*{}'", extension_),
        batch_size_,
        recursive_ ? "yes" : "no");

    load_next_batch();
  }

  [[nodiscard]] std::optional<std::string> next() {
    if (current_index_ == current_batch_.size() && !load_next_batch())
      return std::nullopt;

    return current_batch_[current_index_++];
  }

  [[nodiscard]] std::size_t size() const noexcept { return total_files_; }

  void reset() {
    current_index_ = 0;
    total_files_ = 0;
    exhausted_ = false;
    current_batch_.clear();

    // Recreate iterators from the beginning
    if (recursive_) recursive_it_.emplace(directory_path_, fs::directory_options::skip_permission_denied);
    else                  dir_it_.emplace(directory_path_, fs::directory_options::skip_permission_denied);
    load_next_batch();
  }

private:
  bool load_next_batch() {
    current_batch_.clear();
    current_index_ = 0;
    if (exhausted_) return false;

    const auto push_if_match = [&](const fs::path &p) {
      if (!fs::is_regular_file(p)) return;
      if (!suffix_matches(p, extension_)) return;
      current_batch_.push_back(p.string());
      ++total_files_;
    };

    std::size_t loaded = 0;
    try {
      if (recursive_) {
        while (recursive_it_.has_value() && *recursive_it_ != fs::recursive_directory_iterator{} && loaded < batch_size_) {
          push_if_match((*recursive_it_)->path());
          ++(*recursive_it_);
          ++loaded;
        }

        if (*recursive_it_ == fs::recursive_directory_iterator{}) exhausted_ = true;
      } else {
        while (dir_it_.has_value() && *dir_it_ != fs::directory_iterator{} && loaded < batch_size_) {
          push_if_match((*dir_it_)->path());
          ++(*dir_it_);
          ++loaded;
        }

        if (*dir_it_ == fs::directory_iterator{}) exhausted_ = true;
      }
    } catch (const fs::filesystem_error &e) {
      Logger::get_logger()->warn("Error loading batch: {}", e.what());
      exhausted_ = true;
    }

    return !current_batch_.empty();
  }

  static bool suffix_matches(const fs::path &p, std::string_view ext) noexcept {
      if (ext.empty()) return true;

      fs::path fn_path = p.filename();

      auto* cstr = fn_path.c_str();
      auto  len  = std::char_traits<char>::length(cstr);

      if (len < ext.size()) return false;
      return std::memcmp(cstr + len - ext.size(), ext.data(), ext.size()) == 0;
  }


  // config stuff
  bool        recursive_;
  std::string directory_path_;
  std::string extension_;
  std::size_t batch_size_;

  //
  // iterator state, only one is engaged at a time (recursive or non-recursive)
  // technically, we could use a single iterator packed in some tagged union, or std::variant,
  // but it would be more complex and less readable. - Besian, May 16, 2025
  //
  std::optional<fs::recursive_directory_iterator> recursive_it_;
  std::optional<fs::directory_iterator> dir_it_;

  // batching
  std::vector<std::string> current_batch_;
  std::size_t current_index_ = 0;
  std::size_t total_files_ = 0;
  bool exhausted_ = false;
};

} // namespace lahuta::sources

#endif // LAHUTA_PIPELINE_SOURCES_DIRECTORY_SOURCE_HPP
