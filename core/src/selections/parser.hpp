/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 4> parts{"besian", "", "sejdiu", "@gmail.com"};
 *   std::vector<std::string_view> valid, empty;
 *   std::partition_copy(parts.begin(), parts.end(), std::back_inserter(valid), std::back_inserter(empty),
 *     [](std::string_view s) { return !s.empty(); });
 *   std::string s; for (auto p : valid) s += p; return s;
 * }();
 *
 */

#ifndef PARSER_HPP
#define PARSER_HPP

#include <vector>

#include "selections/nodes.hpp"
#include "selections/operators.hpp"

namespace lahuta {

class Parser {
public:
  Parser(const std::vector<std::string> &tokens)
      : tokens(tokens), current(0), operator_handler() {}

  NodePtr parse_expression() { return parse_exp(0); }

  void print_tokens() const {
    std::cout << "Tokens: ";
    for (const auto &token : tokens) {
      std::cout << token << " ";
    }
    std::cout << std::endl;
  }

private:
  std::vector<std::string> tokens;
  size_t current;
  Operators operator_handler;

  std::string peek() {
    if (current < tokens.size()) {
      return tokens[current];
    }
    return "";
  }

  std::string consume() {
    if (current < tokens.size()) {
      return tokens[current++];
    }
    throw std::runtime_error("Unexpected end of input");
  }

  NodePtr parse_exp(int precedence) {
    NodePtr left = parse_primary();

    while (true) {
      std::string op = peek();
      if (!operator_handler.is_binary_operator(op))
        break;

      int op_precedence = operator_handler.get_precedence(op);
      if (op_precedence < precedence)
        break;

      consume(); // consume the operator
      int next_precedence =
          (operator_handler.get_operator(op) == Operator::And ||
           operator_handler.get_operator(op) == Operator::Or)
              ? op_precedence + 1
              : op_precedence;
      NodePtr right = parse_exp(next_precedence);
      left = std::make_shared<BinaryNode>(operator_handler.get_operator(op),
                                          left, right);
    }
    return left;
  }

  NodePtr parse_primary() {
    std::string token = peek();
    if (operator_handler.is_unary_operator(token)) {
      consume(); // consume 'not'
      return std::make_shared<UnaryNode>(Operator::Not, parse_exp(2));
    } else if (token == "(") {
      consume(); // consume '('
      NodePtr expr = parse_expression();
      consume(); // consume ')'
      return expr;
    } else {
      return parse_term();
    }
  }

  NodePtr parse_term() {
    std::string term = consume(); // e.g., 'resid', 'resname'

    std::vector<std::string> values;
    values.push_back(consume()); // first value (e.g., '4' or 'ASP')

    // Check for a range or list
    std::string next_token = peek();
    if (next_token == "-") {
      consume();                   // consume '-'
      values.push_back(consume()); // second value in range (e.g., '10')
      return std::make_shared<TermNode>(term, values, true); // 'true' for range
    }

    // Consume additional values
    while (is_value_token(peek())) {
      values.push_back(consume());
    }

    // Term node holding a list of values
    return std::make_shared<TermNode>(term, values, false); // 'false' for list
  }

  bool is_value_token(const std::string &token) {
    return !token.empty() && !operator_handler.is_binary_operator(token) &&
           !operator_handler.is_unary_operator(token) && token != "(" &&
           token != ")";
  }
};

} // namespace lahuta
#endif // PARSER_HPP
