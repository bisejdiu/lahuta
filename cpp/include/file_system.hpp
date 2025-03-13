#ifndef LAHUTA_FILE_SYSTEM_HPP
#define LAHUTA_FILE_SYSTEM_HPP

#include <filesystem>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

class FileHandler {
public:
  // takes in a path to a file or directory
  explicit FileHandler(const fs::path &path) : path_(path) {
    if (!fs::exists(path)) {
      throw std::runtime_error("Path does not exist: " + path.string());
    }
  }

  static std::string get_extension(const fs::path &path) {
    fs::path p = path;
    std::string ext;

    while (!p.extension().empty()) {
      ext = p.extension().string() + ext;
      p = p.stem();
    }
    return ext;
  }

  bool is_directory() const { return fs::is_directory(path_); }
  bool is_file() const { return fs::is_regular_file(path_); }
  const fs::path &get_path() const noexcept { return path_; }
  std::uintmax_t get_size() const { return fs::is_regular_file(path_) ? fs::file_size(path_) : 0; }

private:
  fs::path path_;
};

class FileChunk : public std::vector<std::string> {
public:
  using std::vector<std::string>::vector;
  operator bool() const { return !this->empty(); }
};

class DirectoryHandler {
public:
  /// takes in optional extension and recursive flag
  explicit DirectoryHandler(
      const fs::path &dir_path, const std::string extension = {}, bool recursive = false)
      : dir_path_(dir_path), extension_(std::move(extension)), recursive_(recursive) {
    if (!fs::is_directory(dir_path)) {
      throw std::runtime_error("Not a directory: " + dir_path.string());
    }
  }

  // could be useful for the Python binding
  void set_extension(const std::string &extension) { extension_ = extension; }
  void set_recursive(bool recursive) { recursive_ = recursive; }

  class Iterator {
  public:
    using iterator_category = std::input_iterator_tag;
    using value_type = FileHandler;
    using difference_type = std::ptrdiff_t;
    using pointer = FileHandler *;
    using reference = FileHandler &;

    Iterator() = default;

    Iterator(const fs::path &dir_path, const std::string_view extension, bool recursive)
        : extension_(extension), recursive_(recursive) {
      if (recursive_) {
        iter_.emplace<fs::recursive_directory_iterator>(dir_path);
      } else {
        iter_.emplace<fs::directory_iterator>(dir_path);
      }
      advance();
    }

    const FileHandler &operator*() const noexcept { return *current_file_; }
    const FileHandler *operator->() const noexcept { return current_file_.get(); }

    Iterator &operator++() {
      advance();
      return *this;
    }

    bool operator==(const Iterator &other) const {
      if (!iter_.index() && !other.iter_.index()) return true;
      if (iter_.index() != other.iter_.index()) return false;

      if (recursive_) {
        return std::get<fs::recursive_directory_iterator>(iter_)
               == std::get<fs::recursive_directory_iterator>(other.iter_);
      }
      return std::get<fs::directory_iterator>(iter_) == std::get<fs::directory_iterator>(other.iter_);
    }

    bool operator!=(const Iterator &other) const { return !(*this == other); }

  private:
    void advance() {
      while (true) {
        if (recursive_) {
          auto &it = std::get<fs::recursive_directory_iterator>(iter_);
          if (it == fs::recursive_directory_iterator{}) {
            iter_.emplace<std::monostate>();
            return;
          }
          if (fs::is_regular_file(it->path())
              && (extension_.empty() || FileHandler::get_extension(it->path()) == extension_)) {
            current_file_ = std::make_unique<FileHandler>(it->path());
            ++it;
            return;
          }
          ++it;
        } else {
          auto &it = std::get<fs::directory_iterator>(iter_);
          if (it == fs::directory_iterator{}) {
            iter_.emplace<std::monostate>();
            return;
          }
          if (fs::is_regular_file(it->path())
              && (extension_.empty() || FileHandler::get_extension(it->path()) == extension_)) {
            current_file_ = std::make_unique<FileHandler>(it->path());
            ++it;
            return;
          }
          ++it;
        }
      }
    }

    std::variant<std::monostate, fs::directory_iterator, fs::recursive_directory_iterator> iter_;
    std::string_view extension_;
    bool recursive_;
    std::unique_ptr<FileHandler> current_file_;
  };

  FileChunk next_chunk(size_t chunk_size) {
    FileChunk chunk;
    chunk.reserve(chunk_size);

    if (!current_iter_) {
      current_iter_ = std::make_unique<Iterator>(dir_path_, extension_, recursive_);
    }

    while (chunk.size() < chunk_size && *current_iter_ != end()) {
      chunk.push_back((*current_iter_)->get_path().string());
      ++(*current_iter_);
    }

    // comparing *current_iter_ == end() does not work
    if (chunk.empty()) {
      current_iter_.reset();
    }

    return chunk;
  }

  std::vector<std::string> get_all_files() {
    std::vector<std::string> files;
    for (const auto &file : *this) {
      files.push_back(file.get_path().string());
    }
    return files;
  }

  Iterator begin() const { return Iterator(dir_path_, extension_, recursive_); }
  Iterator end() const { return Iterator(); }

  const fs::path &get_path() const noexcept { return dir_path_; }

private:
  fs::path dir_path_;
  std::string extension_;
  bool recursive_;
  std::unique_ptr<Iterator> current_iter_;
  Iterator current_iterator_ = Iterator(dir_path_, extension_, recursive_);
};

inline int test_impl() {
  try {
    DirectoryHandler dir_handler("/Users/bsejdiu/data/PDB_ARCHIVE", ".cif.gz", true);
    for (const auto &file : dir_handler) {
      std::cout << "PDB file: " << file.get_path() << std::endl;
    }

  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }
  return 0;
}

inline std::size_t count_files(const DirectoryHandler &dir_handler) {
  std::size_t count = 0;
  for (const auto &file : dir_handler) {
    ++count;
  }
  return count;
}

#endif // LAHUTA_FILE_SYSTEM_HPP
