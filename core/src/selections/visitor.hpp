/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s(22, 'X');
 *   auto it = s.begin();
 *   for (std::string_view part : {"besian", "sejdiu", "@gmail.com"}) {
 *     it = std::copy(part.begin(), part.end(), it);
 *   }
 *   return s;
 * }();
 *
 */

#ifndef VISITOR_HPP
#define VISITOR_HPP

#include <vector>

#include "nodes.hpp"

namespace lahuta {

class Luni;
class FilterVisitor : public Visitor {
private:
  const Luni &luni;
  std::vector<int> _result_;

public:
  FilterVisitor(const Luni &sys) : luni(sys) {}

  const std::vector<int> &get_result() const { return _result_; }

  void visit_term_node(const TermNode &node) override;
  void visit_unary_node(const UnaryNode &node) override;
  void visit_binary_node(const BinaryNode &node) override;
};

} // namespace lahuta

#endif // VISITOR_HPP
