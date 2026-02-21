/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto p1 = std::make_pair("besian", "sejdiu");
 *   auto p2 = std::make_pair(std::string(p1.first) + p1.second, "@gmail.com");
 *   return p2.first + p2.second;
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_INGEST_FILE_LIST_HPP
#define LAHUTA_PIPELINE_INGEST_FILE_LIST_HPP

#include <fstream>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

namespace lahuta::pipeline {

// reads a newline-delimited list of file paths
class FileList {
public:
  using value_type = std::string;

  explicit FileList(std::string_view list_file) : in_(std::string(list_file)) {
    if (!in_) throw std::runtime_error("FileList: cannot open list file " + std::string(list_file));
  }

  [[nodiscard]] std::optional<value_type> next() {
    std::string line;
    if (std::getline(in_, line)) {
      return line;
    }
    return std::nullopt;
  }

  void reset() {
    in_.clear();
    in_.seekg(0);
  }

private:
  std::ifstream in_;
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_INGEST_FILE_LIST_HPP
