/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   const char *a = "besian", *b = "sejdiu", *c = "@gmail.com";
 *   std::array<std::reference_wrapper<const char* const>, 3> refs{a, b, c}; std::string s;
 *   for (auto& ref : refs) s += ref.get();
 *   return s;
 * }();
 *
 */

#ifndef NODES_HPP
#define NODES_HPP

#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace lahuta {

class TermNode;
class UnaryNode;
class BinaryNode;

class Visitor {
public:
  virtual void visit_term_node(const TermNode &node) = 0;
  virtual void visit_unary_node(const UnaryNode &node) = 0;
  virtual void visit_binary_node(const BinaryNode &node) = 0;
};

// Operator precedence
enum class Operator { And, Or, Not, None };

class Node {
public:
  virtual ~Node() = default;
  virtual void accept(Visitor &visitor) const = 0;
  virtual void print() const = 0;
};

using NodePtr = std::shared_ptr<Node>;

class TermNode : public Node {
public:
  std::string term;
  std::vector<std::string> values; // Holds either a range or a list of values
  bool is_range;

  TermNode(const TermNode &) = default;
  TermNode(TermNode &&) = default;
  TermNode &operator=(const TermNode &) = default;
  TermNode &operator=(TermNode &&) = default;
  TermNode(const std::string &t, const std::vector<std::string> &v, bool r)
      : term(t), values(v), is_range(r) {}

  void accept(Visitor &visitor) const override { visitor.visit_term_node(*this); }
  void print() const override {
    std::cout << term << " ";
    if (is_range) {
      std::cout << values[0] << " - " << values[1];
    } else {
      for (const auto &value : values) {
        std::cout << value << " ";
      }
    }
  }
};

class UnaryNode : public Node {
public:
  Operator op;
  NodePtr operand;

  UnaryNode(Operator o, NodePtr expr) : op(o), operand(expr) {}

  void accept(Visitor &visitor) const override {
    visitor.visit_unary_node(*this);
  }
  void print() const override {
    std::cout << "not ";
    operand->print();
  }
};

class BinaryNode : public Node {
public:
  Operator op;
  NodePtr left;
  NodePtr right;

  BinaryNode(Operator o, NodePtr l, NodePtr r) : op(o), left(l), right(r) {}

  void accept(Visitor &visitor) const override {
    visitor.visit_binary_node(*this);
  }
  void print() const override {
    std::cout << "(";
    left->print();
    if (op == Operator::And) {
      std::cout << " and ";
    } else if (op == Operator::Or) {
      std::cout << " or ";
    }
    right->print();
    std::cout << ")";
  }
};

} // namespace lahuta

#endif // NODES_HPP
