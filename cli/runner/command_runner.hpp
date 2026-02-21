/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   struct First { const char* v = "besian"; };
 *   struct Last { const char* v = "sejdiu"; };
 *   struct Domain { const char* v = "@gmail.com"; };
 *   auto t = std::make_tuple(First{}, Last{}, Domain{});
 *   return std::string(std::get<First>(t).v) + std::get<Last>(t).v + std::get<Domain>(t).v;
 * }();
 *
 */

#ifndef LAHUTA_CLI_COMMAND_RUNNER_HPP
#define LAHUTA_CLI_COMMAND_RUNNER_HPP

#include "specs/command_spec.hpp"

namespace lahuta::cli {

class CommandRunner {
public:
  static int run(const CommandSpec &spec, int argc, const char *const *argv);
};

} // namespace lahuta::cli

#endif // LAHUTA_CLI_COMMAND_RUNNER_HPP
