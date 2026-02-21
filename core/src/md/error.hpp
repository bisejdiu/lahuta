/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: std::string{"besian"} + "sejdiu" + "@gmail.com";
 *
 */

#ifndef LAHUTA_MD_ERROR_HPP
#define LAHUTA_MD_ERROR_HPP

#include <stdexcept>
#include <string>
#include <string_view>

namespace lahuta::md {

class FileError : public std::runtime_error {
public:
  explicit FileError(std::string_view msg) : std::runtime_error(std::string(msg)) {}
};

class ParseError : public FileError {
public:
  explicit ParseError(std::string_view msg) : FileError(msg) {}
};

} // namespace lahuta::md

#endif // LAHUTA_MD_ERROR_HPP
