#ifndef LAHUTA_PIPELINE_SOURCES_FILE_READER_SOURCE_HPP
#define LAHUTA_PIPELINE_SOURCES_FILE_READER_SOURCE_HPP

#include "file_reader.hpp"
#include "logging.hpp"
#include <memory>
#include <optional>
#include <string>

// clang-format off
namespace lahuta::sources {

// reading serialized objects from a file.
template <typename FormatTag, typename T>
class FileReaderSource {
public:
  using value_type = std::shared_ptr<const T>;

  explicit FileReaderSource(const std::string &file_path, std::size_t batch_size = 1024)
      : file_path_(file_path), batch_size_(batch_size), reader_(std::make_shared<FileReader<FormatTag, T>>(file_path)), current_index_(0) {
    Logger::get_logger()->info("FileReaderSource: reading from {}; batch_size={}", file_path_, batch_size_);
    load_next_batch();
  }

  FileReaderSource(const FileReaderSource &) = default;
  FileReaderSource &operator=(const FileReaderSource &) = default;

  [[nodiscard]] std::optional<value_type> next() {
    if (!reader_ || (current_index_ >= current_batch_.size() && !load_next_batch())) return std::nullopt;
    return current_batch_[current_index_++];
  }

  void reset() {
    current_index_ = 0;
    current_batch_.clear();
    if (reader_) reader_->reset();
    load_next_batch();
  }

private:
  bool load_next_batch() {
    current_batch_.clear();
    current_index_ = 0;

    if (!reader_ || reader_->eof()) return false;

    std::size_t count = 0;
    while (count < batch_size_) {
      T item;
      if (!reader_->read_next(item)) break;
      current_batch_.push_back(std::make_shared<const T>(std::move(item)));
      ++count;
    }
    return !current_batch_.empty();
  }

  // config stuff
  std::string file_path_;
  std::size_t batch_size_;
  std::shared_ptr<FileReader<FormatTag, T>> reader_;

  // batching
  std::vector<value_type> current_batch_;
  std::size_t current_index_;
};

} // namespace lahuta::sources

#endif // LAHUTA_PIPELINE_SOURCES_FILE_READER_SOURCE_HPP
