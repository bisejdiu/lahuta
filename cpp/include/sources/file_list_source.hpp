#ifndef LAHUTA_PIPELINE_SOURCES_FILE_LIST_SOURCE_HPP
#define LAHUTA_PIPELINE_SOURCES_FILE_LIST_SOURCE_HPP

#include <fstream>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace lahuta::sources {

// reads a newline-delimited list of file paths
class FileListSource {
public:
  struct FileData {
    std::string path;     // full path
    std::string contents; // whole file in memory
  };

  using value_type = FileData;

  explicit FileListSource(std::string_view list_file, std::size_t io_threads = 0) {
    std::ifstream in(std::string(list_file).c_str());
    if (!in) {
      throw std::runtime_error("FileListSource: cannot open list file " + std::string(list_file));
    }

    std::string line;
    while (std::getline(in, line)) {
      if (!line.empty()) {
        paths_.push_back(std::move(line));
      }
    }
  }

  // TODO: I'm not sure if we still need these ctors
  FileListSource(FileListSource &&other) noexcept : paths_(std::move(other.paths_)), idx_(other.idx_) {}
  FileListSource &operator=(FileListSource &&other) noexcept {
    if (this != &other) {
      paths_ = std::move(other.paths_);
      idx_ = other.idx_;
    }
    return *this;
  }

  FileListSource(const FileListSource &) = delete;
  FileListSource &operator=(const FileListSource &) = delete;

  std::optional<value_type> next() {
    const std::size_t i = idx_++;
    if (i >= paths_.size()) return std::nullopt;

    return load_file(paths_[i]);
  }

private:
  FileData load_file(std::string_view path) {
    std::ifstream in(std::string(path).c_str(), std::ios::binary | std::ios::ate);
    if (!in) {
      throw std::runtime_error("FileListSource: cannot open " + std::string(path));
    }

    const auto size = static_cast<std::size_t>(in.tellg());
    std::string contents(size, '\0');
    in.seekg(0);
    in.read(contents.data(), size);

    return {std::string(path), std::move(contents)};
  }

  std::size_t idx_{0};
  std::vector<std::string> paths_;
};

} // namespace lahuta::sources

#endif // LAHUTA_PIPELINE_SOURCES_FILE_LIST_SOURCE_HPP
