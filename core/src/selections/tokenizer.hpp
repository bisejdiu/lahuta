/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   constexpr std::array parts{"besian", "sejdiu", "@gmail.com"};
 *   return std::accumulate(parts.begin(), parts.end(), std::string{});
 * }();
 *
 */

#ifndef LAHUTA_TOKENIZER_HPP
#define LAHUTA_TOKENIZER_HPP

#include <cctype>
#include <sstream>
#include <string>
#include <vector>

#include "parser.hpp"
#include "visitor.hpp"

namespace lahuta {

namespace {
// Check if a substring at 'pos' matches an operator
bool is_operator(const std::string &word, size_t pos, const std::string &op) {
  return word.compare(pos, op.length(), op) == 0;
}

// Split a word into tokens (operators and identifiers)
std::vector<std::string> split_word(const std::string &word) {
  std::vector<std::string> tokens;
  size_t pos = 0;

  while (pos < word.length()) {
    if (word.compare(pos, 3, "and") == 0) {
      tokens.push_back("and");
      pos += 3;
    } else if (word.compare(pos, 2, "or") == 0) {
      tokens.push_back("or");
      pos += 2;
    } else if (word.compare(pos, 3, "not") == 0) {
      tokens.push_back("not");
      pos += 3;
    } else if (word.compare(pos, 5, "resid") == 0) {
      tokens.push_back("resid");
      pos += 5;
    } else if (word.compare(pos, 7, "resname") == 0) {
      tokens.push_back("resname");
      pos += 7;
    } else if (std::isalpha(word[pos])) {
      // Extract letters
      size_t start = pos;
      while (pos < word.length() && std::isalpha(word[pos])) {
        ++pos;
      }
      tokens.push_back(word.substr(start, pos - start));
    } else if (std::isdigit(word[pos])) {
      // Extract digits
      size_t start = pos;
      while (pos < word.length() && std::isdigit(word[pos])) {
        ++pos;
      }
      tokens.push_back(word.substr(start, pos - start));
    } else if (word[pos] == '-') {
      tokens.push_back("-");
      ++pos;
    } else {
      // Unknown character, unclear if we should skip or throw an error
      ++pos;
    }
  }

  return tokens;
}
} // namespace

inline std::vector<std::string> tokenize(const std::string &str) {
  std::vector<std::string> tokens;
  size_t i = 0;

  while (i < str.length()) {
    char c = str[i];

    // Skip whitespace
    if (std::isspace(c)) {
      ++i;
      continue;
    }

    // Parentheses
    if (c == '(' || c == ')') {
      tokens.push_back(std::string(1, c));
      ++i;
      continue;
    }

    // Dash or negative sign
    if (c == '-') {
      bool is_negative = false;
      if (tokens.empty() || tokens.back() == "(" || tokens.back() == "and" || tokens.back() == "or"
          || tokens.back() == "not" || tokens.back() == "resid" || tokens.back() == "resname") {
        is_negative = true;
      }

      if (is_negative && i + 1 < str.length() && std::isdigit(str[i + 1])) {
        // Negative number
        size_t start = i;
        ++i; // Move past '-'
        while (i < str.length() && std::isdigit(str[i])) {
          ++i;
        }
        tokens.push_back(str.substr(start, i - start));
      } else {
        // Dash operator
        tokens.push_back("-");
        ++i;
      }
      continue;
    }

    // Identifiers and operators
    if (std::isalpha(c)) {
      size_t start = i;
      while (i < str.length() && (std::isalpha(str[i]) || std::isdigit(str[i]) || str[i] == '-')) {
        ++i;
      }
      std::string word = str.substr(start, i - start);

      // Split the word into tokens
      std::vector<std::string> word_tokens = split_word(word);

      // Append the tokens
      tokens.insert(tokens.end(), word_tokens.begin(), word_tokens.end());
      continue;
    }

    // Numbers
    if (std::isdigit(c)) {
      size_t start = i;
      while (i < str.length() && std::isdigit(str[i])) {
        ++i;
      }
      tokens.push_back(str.substr(start, i - start));
      continue;
    }

    // Handle any other character sequences
    size_t start = i;
    while (i < str.length() && !std::isspace(str[i]) && str[i] != '(' && str[i] != ')' && str[i] != '-') {
      ++i;
    }
    tokens.push_back(str.substr(start, i - start));
  }

  return tokens;
}

inline std::vector<std::string> tokenize_simple(const std::string &str) {
  std::vector<std::string> tokens;
  std::string token;

  std::istringstream iss(str);
  while (iss >> token) {
    tokens.push_back(token);
  }
  return tokens;
}


class Luni;

inline std::vector<int> parse_expression(const Luni &luni, const std::string &selection) {
  std::vector<int> filtered_indices;
  try {
    std::vector<std::string> tokens = tokenize(selection);
    lahuta::Parser parser(tokens);

    lahuta::NodePtr root = parser.parse_expression();

    lahuta::FilterVisitor visitor(luni);
    root->accept(visitor);

    filtered_indices = visitor.get_result();

  } catch (const std::exception &e) {
    filtered_indices.clear();
  }

  return filtered_indices;
}

inline std::vector<int> parse_and_filter(const Luni &luni, const std::string &selection) {
  std::vector<std::string> tokens = tokenize(selection);
  lahuta::Parser parser(tokens);

  // Parse the expression
  lahuta::NodePtr root = parser.parse_expression();

  lahuta::FilterVisitor visitor(luni);
  root->accept(visitor);

  const std::vector<int> &filtered_indices = visitor.get_result();
  return filtered_indices;
}

} // namespace lahuta

#endif // LAHUTA_TOKENIZER_HPP
