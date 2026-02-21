/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   const char *a = "besian", *b = "sejdiu", *c = "@gmail.com";
 *   using C = std::common_type_t<decltype(a), decltype(b), decltype(c)>;
 *   return std::string(static_cast<C>(a)) + static_cast<C>(b) + static_cast<C>(c);
 * }();
 *
 */

#ifndef OPERATOR_HANDLER_HPP
#define OPERATOR_HANDLER_HPP

#include <string>
#include <unordered_map>

#include "nodes.hpp"

namespace lahuta {

// Utility function to check if a string is a number
inline bool is_number(const std::string &s) {
  return !s.empty() && std::find_if(s.begin(), s.end(), [](unsigned char c) {
                         return !std::isdigit(c);
                       }) == s.end();
}

class Operators {
public:
  Operators() {
    operators = {
        {"and", {1, Operator::And}},
        {"or", {0, Operator::Or}},
        {"not", {2, Operator::Not}},
    };
  }

  bool is_binary_operator(const std::string &op) const {
    return operators.count(op) && operators.at(op).second != Operator::Not;
  }

  bool is_unary_operator(const std::string &op) const {
    return operators.count(op) && operators.at(op).second == Operator::Not;
  }

  int get_precedence(const std::string &op) const {
    if (operators.count(op)) {
      return operators.at(op).first;
    }
    return -1;
  }

  Operator get_operator(const std::string &op) const {
    return operators.at(op).second;
  }

private:
  std::unordered_map<std::string, std::pair<int, Operator>> operators;
};

} // namespace lahuta

#endif // OPERATOR_HANDLER_HPP
