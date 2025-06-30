#ifndef LAHUTA_IO_GZIP_FILE_SPILL_POLICY_HPP
#define LAHUTA_IO_GZIP_FILE_SPILL_POLICY_HPP

#include "serialization/serializer.hpp"
#include <cstring>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
#include <zlib.h>

// clang-format off
namespace lahuta {

template<
  typename  FormatTag,
  typename  Payload
>
class GzipFileSpillPolicy {
  static_assert(
    std::is_same<Payload, std::shared_ptr<const typename Payload::element_type>>::value,
    "GzipFileSpillPolicy expects Payload = std::shared_ptr<const U>");

  using value_type      = typename Payload::element_type;
  using serializer_type = serialization::Serializer<FormatTag, value_type>;

public:
  using payload_type = Payload;

  explicit GzipFileSpillPolicy(std::size_t batch, const std::string& path, int compression_level = Z_DEFAULT_COMPRESSION)
    : batch_(batch), data_size_(0), out_(gzopen(path.c_str(), "wb")) {
    if (!out_) throw std::runtime_error("Cannot open gzip file: " + path);

    if (compression_level != Z_DEFAULT_COMPRESSION) {
      gzsetparams(out_, compression_level, Z_DEFAULT_STRATEGY);
    }
  }

  template<size_t N>
  explicit GzipFileSpillPolicy(std::size_t batch, const char (&path)[N], int compression_level = Z_DEFAULT_COMPRESSION)
    : batch_(batch), data_size_(0), out_(gzopen(path, "wb")) {
    if (!out_) {
      std::string path_str(path);
      throw std::runtime_error("Cannot open gzip file: " + path_str);
    }

    if (compression_level != Z_DEFAULT_COMPRESSION) {
      gzsetparams(out_, compression_level, Z_DEFAULT_STRATEGY);
    }
  }

  ~GzipFileSpillPolicy() {
    if (out_) {
      finish();
      gzclose(out_);
    }
  }

  bool needs_flush() const {
    return buf_.size() >= batch_ || data_size_ > 1 * 1024 * 1024;
  }

  std::size_t append_size(const payload_type& p) const {
    return serializer_type::serialize(*p).size();
  }

  void append(const payload_type& p) { push(p); }
  void append(payload_type&& p)      { push(p); }

  std::size_t flush() {
    if (buf_.empty()) return 0;

    std::size_t freed = data_size_;
    for (const auto& raw : buf_) {
      gzwrite(out_, raw.data(), raw.size());
      gzputc(out_, '\n');
    }
    gzflush(out_, Z_SYNC_FLUSH);
    buf_.clear();
    data_size_ = 0;
    return freed;
  }

  std::size_t finish() {
      return flush();
  }

private:
  void push(const payload_type& p) {
    auto raw = serializer_type::serialize(*p);
    data_size_ += raw.size();
    buf_.emplace_back(std::move(raw));
  }

  std::size_t              batch_;
  std::size_t              data_size_;
  gzFile                   out_;
  std::vector<std::string> buf_;
};

} // namespace lahuta 

#endif // LAHUTA_IO_GZIP_FILE_SPILL_POLICY_HPP
