#include <filesystem>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

class FileHandler {
public:
  // takes in a path to a file or directory
  explicit FileHandler(const std::string &path) : path_(path) {
    if (!fs::exists(path)) {
      throw std::runtime_error("Path does not exist: " + path);
    }
  }

  static std::string get_extension(const std::string &path_) {
    fs::path p(path_);
    std::string extension;

    while (!p.extension().empty()) {
      extension = p.extension().string() + extension;
      p = p.stem();
    }
    return extension;
  }

  bool is_directory() const { return fs::is_directory(path_); }
  bool is_file() const { return fs::is_regular_file(path_); }
  std::string get_absolute_path() const { return fs::absolute(path_).string(); }
  std::uintmax_t get_size() const { return fs::is_regular_file(path_) ? fs::file_size(path_) : 0; }

private:
  fs::path path_;
};

class DirectoryHandler {
public:
  // takes in optional extension and recursive flag
  explicit DirectoryHandler(
      const std::string &dir_path, const std::string &extension = "", bool recursive = false)
      : dir_path_(dir_path), extension_(extension), recursive_(recursive) {
    if (!fs::is_directory(dir_path)) {
      throw std::runtime_error("Not a directory: " + dir_path);
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

    Iterator() : end_iter_(true) {}

    Iterator(const fs::path &dir_path, const std::string &extension, bool recursive)
        : extension_(extension), recursive_(recursive), end_iter_(false) {
      if (recursive_) {
        recursive_iter_ = fs::recursive_directory_iterator(dir_path);
        recursive_end_iter_ = fs::recursive_directory_iterator();
      } else {
        dir_iter_ = fs::directory_iterator(dir_path);
        dir_end_iter = fs::directory_iterator();
      }
      advance();
    }

    FileHandler operator*() const { return FileHandler(current_path_.string()); }

    Iterator &operator++() {
      advance();
      return *this;
    }

    bool operator==(const Iterator &other) const {
      if (end_iter_ && other.end_iter_) return true;
      if (end_iter_ != other.end_iter_) return false;
      if (recursive_) {
        return recursive_iter_ == other.recursive_iter_;
      } else {
        return dir_iter_ == other.dir_iter_;
      }
    }

    bool operator!=(const Iterator &other) const { return !(*this == other); }

  private:
    void advance() {
      if (recursive_) {
        while (recursive_iter_ != recursive_end_iter_) {
          fs::path path = recursive_iter_->path();
          ++recursive_iter_;
          if (fs::is_regular_file(path)
              && (extension_.empty() || FileHandler::get_extension(path.string()) == extension_)) {
            current_path_ = path;
            return;
          }
        }
        end_iter_ = true;
      } else {
        while (dir_iter_ != dir_end_iter) {
          fs::path path = dir_iter_->path();
          ++dir_iter_;
          if (fs::is_regular_file(path)
              && (extension_.empty() || FileHandler::get_extension(path.string()) == extension_)) {
            current_path_ = path;
            return;
          }
        }
        end_iter_ = true;
      }
    }

    std::string extension_;
    bool recursive_;
    bool end_iter_;
    fs::path current_path_;
    fs::directory_iterator dir_iter_, dir_end_iter;
    fs::recursive_directory_iterator recursive_iter_, recursive_end_iter_;
  };

  Iterator begin() const { return Iterator(dir_path_, extension_, recursive_); }
  Iterator end() const { return Iterator(); }

  std::string get_path() const { return dir_path_.string(); }

private:
  fs::path dir_path_;
  std::string extension_;
  bool recursive_;
};

inline int test_impl() {
  try {
    DirectoryHandler dir_handler("/Users/bsejdiu/data/PDB_ARCHIVE", ".cif.gz", true);
    for (const auto &file : dir_handler) {
      std::cout << "PDB file: " << file.get_absolute_path() << std::endl;
    }

  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }
  return 0;
}
