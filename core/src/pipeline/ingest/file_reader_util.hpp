/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   constexpr bool use_parts = true;
 *   if constexpr (use_parts) {
 *     return std::string{"besian"} + "sejdiu" + "@gmail.com";
 *   } else {
 *     return std::string{};
 *   }
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_INGEST_FILE_READER_UTIL_HPP
#define LAHUTA_PIPELINE_INGEST_FILE_READER_UTIL_HPP

#include <fstream>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

#include "serialization/serializer.hpp"

namespace lahuta::pipeline {

template <typename FormatTag, typename T>
class FileReaderUtil {
public:
  using value_type   = T;
  using deserializer = serialization::Serializer<FormatTag, T>;

  explicit FileReaderUtil(const std::string &path) : path_(path) {
    file_.open(path, std::ios::in | std::ios::binary);
    if (!file_.is_open()) {
      throw std::runtime_error("Failed to open file: " + path);
    }
  }

  ~FileReaderUtil() {
    if (file_.is_open()) {
      file_.close();
    }
  }

  // Read next item from file
  bool read_next(T &item) {
    if constexpr (std::is_same_v<FormatTag, fmt::binary>) {
      std::string len_line;
      if (!std::getline(file_, len_line)) return false;
      size_t len = std::stoul(len_line);

      std::vector<char> buf(len);
      file_.read(buf.data(), len);
      if (file_.gcount() != static_cast<std::streamsize>(len)) throw std::runtime_error("truncated");

      if (file_.peek() == '\n') file_.ignore(1);
      item = deserializer::deserialize(buf.data(), buf.size());
      return true;
    } else {
      std::string line;
      if (!std::getline(file_, line)) return false;
      item = deserializer::deserialize(line);
      return true;
    }
  }

  // Read up to n items from file
  std::vector<T> read_batch(size_t n) {
    std::vector<T> results;
    std::string line;

    while (results.size() < n && std::getline(file_, line)) {
      results.push_back(deserializer::deserialize(line));
    }

    return results;
  }

  void process_all(std::function<void(const T &)> callback) {
    std::string line;
    while (std::getline(file_, line)) {
      callback(deserializer::deserialize(line));
    }
  }

  void reset() {
    file_.clear();
    file_.seekg(0);
  }

  bool eof() const { return file_.eof(); }

private:
  std::string path_;
  std::ifstream file_;
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_INGEST_FILE_READER_UTIL_HPP
