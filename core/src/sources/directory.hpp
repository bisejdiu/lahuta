#ifndef LAHUTA_SOURCES_DIRECTORY_HPP
#define LAHUTA_SOURCES_DIRECTORY_HPP

#include <algorithm>
#include <cstring>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "logging.hpp"

// clang-format off
namespace lahuta::sources {
namespace fs = std::filesystem;

// batch-oriented directory iterator.
class Directory {
public:
  using value_type = std::string;

  explicit Directory(
      std::string_view directory_path,
      std::string_view extension = "",
      bool recursive = true,
      std::size_t batch_size = 1024)
      : Directory(directory_path,
                  extension.empty() ? std::vector<std::string>{}
                                    : std::vector<std::string>{std::string(extension)},
                  recursive,
                  batch_size) {}

  Directory(std::string_view directory_path,
            std::vector<std::string> extensions,
            bool recursive = true,
            std::size_t batch_size = 1024)
      : directory_path_(directory_path),
        extensions_(normalize_extensions(std::move(extensions))),
        recursive_(recursive),
        batch_size_(batch_size) {

    if (!fs::exists(directory_path_) || !fs::is_directory(directory_path_)) {
      throw fs::filesystem_error( "Directory does not exist", std::make_error_code(std::errc::no_such_file_or_directory));
    }

    if (recursive_) recursive_it_.emplace(directory_path_, fs::directory_options::skip_permission_denied);
    else                  dir_it_.emplace(directory_path_, fs::directory_options::skip_permission_denied);

    const auto filters = describe_filters();
    Logger::get_logger()->info(
        "Directory: scanning directory '{}'{}; batch_size={} (recursive={})",
        directory_path_,
        filters.empty() ? std::string{} : fmt::format(", filter='{}'", filters),
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
    total_files_   = 0;
    exhausted_     = false;
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
      if (!suffix_matches(p, extensions_)) return;
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

  static std::vector<std::string> normalize_extensions(std::vector<std::string> exts) {
    std::vector<std::string> filtered;
    filtered.reserve(exts.size());
    for (auto &ext : exts) {
      if (ext.empty()) continue;
      filtered.emplace_back(std::move(ext));
    }
    std::sort(filtered.begin(), filtered.end());
    filtered.erase(std::unique(filtered.begin(), filtered.end()), filtered.end());
    return filtered;
  }

  [[nodiscard]] std::string describe_filters() const {
    if (extensions_.empty()) return {};
    std::string buf;
    for (std::size_t i = 0; i < extensions_.size(); ++i) {
      if (i > 0) buf.append(", ");
      buf.push_back('*');
      buf.append(extensions_[i]);
    }
    return buf;
  }

  static bool suffix_matches(const fs::path &p, const std::vector<std::string>& exts) noexcept {
      if (exts.empty()) return true;

      const auto filename = p.filename().string();
      for (const auto &ext : exts) {
        if (filename.size() < ext.size()) continue;
        if (std::memcmp(filename.data() + filename.size() - ext.size(), ext.data(), ext.size()) == 0) {
          return true;
        }
      }
      return false;
  }


  // config stuff
  std::string directory_path_;
  std::vector<std::string> extensions_;
  bool        recursive_;
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

#endif // LAHUTA_SOURCES_DIRECTORY_HPP
