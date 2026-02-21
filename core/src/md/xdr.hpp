/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s = "besian@gmail.com";
 *   s.insert(6, "sejdiu");
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_MD_XDR_HPP
#define LAHUTA_MD_XDR_HPP

#include <cstddef>
#include <cstdint>
#include <istream>
#include <ostream>
#include <string>
#include <string_view>

namespace lahuta::md::internal {

class XDRReader {
public:
  using block_type = unsigned int;

  explicit XDRReader(std::istream *s);
  ~XDRReader() = default;

  std::istream *get() { return stream_; }

  template <typename T> bool read(T *p);
  bool read(double *p);
  bool read(std::int64_t *p);
  bool read(std::uint64_t *p);
  template <typename T> bool read(T *ary, std::size_t n);
  bool read(char *p, std::size_t n);
  bool read(std::string &s);

private:
  std::istream *stream_;
  bool need_to_swab_;

  template <typename T> constexpr T swab(T value);
};

class XDRWriter {
public:
  using block_type = unsigned int;

  XDRWriter();
  explicit XDRWriter(std::ostream *s);
  ~XDRWriter() = default;

  std::ostream *get() { return stream_; }
  void set_stream(std::ostream *o) { stream_ = o; }

  template <typename T> bool write(const T &p);
  bool write(const double &p);
  bool write(const std::uint64_t &p);
  template <typename T> bool write(const T *ary, std::size_t n);
  bool write(const char *p, std::size_t n);
  bool write(const char *p);
  bool write(std::string_view s);

private:
  std::ostream *stream_;
  bool need_to_swab_;

  template <typename T> constexpr T swab(T value);
};

} // namespace lahuta::md::internal

#endif // LAHUTA_MD_XDR_HPP
