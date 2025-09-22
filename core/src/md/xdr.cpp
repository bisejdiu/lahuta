#include <cstring>
#include <limits>
#include <stdexcept>

#include "md/xdr.hpp"

namespace lahuta::md::internal {

XDRReader::XDRReader(std::istream *s) : stream_(s), need_to_swab_(false) {
  int test = 0x1234;
  if (*(reinterpret_cast<char *>(&test)) == 0x34) {
    need_to_swab_ = true;
  }
}

template <typename T> bool XDRReader::read(T *p) {
  if (sizeof(T) > sizeof(block_type)) {
    throw std::runtime_error("XDR data size error");
  }

  T result;
  stream_->read(reinterpret_cast<char *>(&result), sizeof(block_type));
  if (sizeof(T) > 1 && need_to_swab_) {
    result = swab<T>(result);
  }

  *p = result;
  return !stream_->fail();
}

bool XDRReader::read(double *p) {
  double result;
  stream_->read(reinterpret_cast<char *>(&result), sizeof(double));
  if (need_to_swab_) {
    result = swab<double>(result);
  }

  *p = result;
  return !stream_->fail();
}

bool XDRReader::read(std::int64_t *p) {
  std::uint64_t raw = 0;
  stream_->read(reinterpret_cast<char *>(&raw), sizeof(raw));
  if (need_to_swab_) {
    raw = swab<std::uint64_t>(raw);
  }
  *p = static_cast<std::int64_t>(raw);
  return !stream_->fail();
}

bool XDRReader::read(std::uint64_t *p) {
  std::uint64_t raw = 0;
  stream_->read(reinterpret_cast<char *>(&raw), sizeof(raw));
  if (need_to_swab_) {
    raw = swab<std::uint64_t>(raw);
  }
  *p = raw;
  return !stream_->fail();
}

template <typename T> bool XDRReader::read(T *ary, std::size_t n) {
  for (std::size_t i = 0; i < n; ++i) {
    if (!read(ary + i)) {
      return false;
    }
  }
  return true;
}

bool XDRReader::read(char *p, std::size_t n) {
  std::size_t rndup = n % sizeof(block_type);
  if (rndup > 0) {
    rndup = sizeof(block_type) - rndup;
  }

  stream_->read(p, static_cast<std::streamsize>(n));
  if (stream_->fail()) {
    return false;
  }

  if (rndup) {
    char buf[sizeof(block_type)];
    stream_->read(buf, static_cast<std::streamsize>(rndup));
  }

  return true;
}

bool XDRReader::read(std::string &s) {
  unsigned int len = 0;
  if (!read(&len)) {
    return false;
  }
  if (len > s.max_size()) {
    throw std::runtime_error("XDR string length exceeds max size");
  }
  s.resize(len);
  return read(&s[0], static_cast<std::size_t>(len));
}

template <typename T> constexpr T XDRReader::swab(T value) {
  T result;
  auto *src = reinterpret_cast<char *>(&value);
  auto *dst = reinterpret_cast<char *>(&result);

  for (std::size_t i = 0; i < sizeof(T); ++i) {
    dst[i] = src[sizeof(T) - 1 - i];
  }

  return result;
}

XDRWriter::XDRWriter() : stream_(nullptr), need_to_swab_(false) {
  int test = 0x1234;
  if (*(reinterpret_cast<char *>(&test)) == 0x34) {
    need_to_swab_ = true;
  }
}

XDRWriter::XDRWriter(std::ostream *s) : stream_(s), need_to_swab_(false) {
  int test = 0x1234;
  if (*(reinterpret_cast<char *>(&test)) == 0x34) {
    need_to_swab_ = true;
  }
}

template <typename T> bool XDRWriter::write(const T &p) {
  if (sizeof(T) > sizeof(block_type)) {
    throw std::runtime_error("XDR data size error");
  }

  block_type u;
  auto *up = reinterpret_cast<T *>(&u);
  *up = p;

  if (sizeof(T) > 1 && need_to_swab_) {
    u = swab<block_type>(u);
  }

  stream_->write(reinterpret_cast<char *>(&u), sizeof(block_type));
  return !stream_->bad();
}

bool XDRWriter::write(const double &p) {
  unsigned long u;
  auto *up = reinterpret_cast<double *>(&u);
  *up = p;

  if (need_to_swab_) {
    u = swab<unsigned long>(u);
  }

  stream_->write(reinterpret_cast<char *>(&u), sizeof(double));
  return !stream_->fail();
}

bool XDRWriter::write(const std::uint64_t &p) {
  std::uint64_t raw = p;
  if (need_to_swab_) {
    raw = swab<std::uint64_t>(raw);
  }
  stream_->write(reinterpret_cast<char *>(&raw), sizeof(raw));
  return !stream_->fail();
}

template <typename T> bool XDRWriter::write(const T *ary, std::size_t n) {
  for (std::size_t i = 0; i < n; ++i) {
    if (!write(ary[i])) {
      return false;
    }
  }
  return true;
}

bool XDRWriter::write(const char *p, std::size_t n) {
  std::size_t rndup = n % sizeof(block_type);
  if (rndup > 0) {
    rndup = sizeof(block_type) - rndup;
  }

  char buf[sizeof(block_type)] = {0};

  stream_->write(p, static_cast<std::streamsize>(n));
  if (!stream_->fail()) {
    stream_->write(buf, static_cast<std::streamsize>(rndup));
  }

  return !stream_->fail();
}

bool XDRWriter::write(const char *p) {
  if (p == nullptr) {
    unsigned int len = 0;
    return write(len);
  }
  std::size_t n = std::strlen(p);
  if (n > std::numeric_limits<unsigned int>::max()) {
    throw std::runtime_error("XDR string too long to encode");
  }
  unsigned int len = static_cast<unsigned int>(n);
  if (!write(len)) {
    return false;
  }
  return write(p, len);
}

bool XDRWriter::write(std::string_view s) {
  if (s.size() > std::numeric_limits<unsigned int>::max()) {
    throw std::runtime_error("XDR string too long to encode");
  }
  unsigned int len = static_cast<unsigned int>(s.size());
  if (!write(len)) {
    return false;
  }
  return write(s.data(), len);
}

template <typename T> constexpr T XDRWriter::swab(T value) {
  T result;
  auto *src = reinterpret_cast<char *>(&value);
  auto *dst = reinterpret_cast<char *>(&result);

  for (std::size_t i = 0; i < sizeof(T); ++i) {
    dst[i] = src[sizeof(T) - 1 - i];
  }

  return result;
}

template bool XDRReader::read<int>(int *p);
template bool XDRReader::read<float>(float *p);
template bool XDRReader::read<unsigned int>(unsigned int *p);

template bool XDRReader::read<int>(int *ary, std::size_t n);
template bool XDRReader::read<float>(float *ary, std::size_t n);

template bool XDRWriter::write<int>(const int &p);
template bool XDRWriter::write<float>(const float &p);
template bool XDRWriter::write<unsigned int>(const unsigned int &p);

template bool XDRWriter::write<int>(const int *ary, std::size_t n);
template bool XDRWriter::write<float>(const float *ary, std::size_t n);

template bool XDRWriter::write<unsigned int>(const unsigned int *ary, std::size_t n);

} // namespace lahuta::md::internal
