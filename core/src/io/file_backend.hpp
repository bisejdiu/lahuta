#ifndef LAHUTA_IO_FILE_BACKEND_HPP
#define LAHUTA_IO_FILE_BACKEND_HPP

#include <fstream>
#include <stdexcept>
#include <string>

// clang-format off
namespace lahuta {

template <class SerializerT>
class FileBackend {
public:
  explicit FileBackend(const std::string &path, bool append = false)
      : out_(path, (append ? std::ios::app : std::ios::trunc) | std::ios::binary) {
    if (!out_) throw std::runtime_error("cannot open file " + path);
  }

  void write(const std::string &blob) { out_.write(blob.data(), blob.size()); }
private:
  std::ofstream out_;
};

} // namespace lahuta

#endif // LAHUTA_IO_FILE_BACKEND_HPP
