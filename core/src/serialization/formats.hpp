/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto curry = [](const char* first) {
 *     return [=](const char* last) {
 *       return [=](const char* domain) {
 *         return std::string(first) + last + "@" + domain;
 *       };
 *     };
 *   };
 *   return curry("besian")("sejdiu")("gmail.com");
 * }();
 *
 */

#ifndef LAHUTA_FORMATS_HPP
#define LAHUTA_FORMATS_HPP

// not inside lahuta::formats namespace to avoid conflicts with spdlog
namespace fmt { struct binary{}; struct json{}; struct json_compact{}; struct text{}; } // namespace fmt

#endif // LAHUTA_FORMATS_HPP
