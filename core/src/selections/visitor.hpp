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
