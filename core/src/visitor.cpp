#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

#include "lahuta.hpp"
#include "selections/visitor.hpp"

namespace lahuta {

enum class TermType { RESID, RESNAME, UNKNOWN };

namespace {
TermType get_term_type(const std::string &term) {
  if (term == "resid")
    return TermType::RESID;
  if (term == "resname")
    return TermType::RESNAME;
  return TermType::UNKNOWN;
}
} // namespace

void FilterVisitor::visit_term_node(const TermNode &node) {
  std::vector<int> filtered_indices;

  switch (get_term_type(node.term)) {
  case TermType::RESID: {
    const auto &resids = luni.resids();
    if (node.is_range) {
      // Handle range: e.g., resid 4 - 10
      int start = std::stoi(node.values[0]);
      int end = std::stoi(node.values[1]);

      for (size_t i = 0; i < resids.size(); ++i) {
        if (resids[i] >= start && resids[i] <= end) {
          filtered_indices.push_back(i);
        }
      }
    } else {
      // Handle list: e.g., resid 4 5 6 7
      std::unordered_set<int> list_values;
      for (const std::string &val : node.values) {
        list_values.insert(std::stoi(val));
      }

      for (size_t i = 0; i < resids.size(); ++i) {
        if (list_values.find(resids[i]) != list_values.end()) {
          filtered_indices.push_back(i);
        }
      }
    }
    break;
  }
  case TermType::RESNAME: {
    // Handle 'resname' term (list syntax only)
    const auto &resnames = luni.resnames();
    std::unordered_set<std::string> name_set(node.values.begin(),
                                             node.values.end());

    for (size_t i = 0; i < resnames.size(); ++i) {
      if (name_set.find(resnames[i]) != name_set.end()) {
        filtered_indices.push_back(i);
      }
    }
    break;
  }
  default:
    throw std::runtime_error("Unsupported term: " + node.term);
  }
  _result_ = std::move(filtered_indices);
}

void FilterVisitor::visit_unary_node(const UnaryNode &node) {
  // Visit the operand of the unary node first
  node.operand->accept(*this);
  std::vector<int> operand_result = _result_;
  _result_.clear();

  // Apply the 'not' operator to the result by negating the selection
  auto resids = luni.resids();
  for (size_t i = 0; i < resids.size(); ++i) {
    if (std::find(operand_result.begin(), operand_result.end(), i) ==
        operand_result.end()) {
      _result_.push_back(i);
    }
  }
}

void FilterVisitor::visit_binary_node(const BinaryNode &node) {
  // Visit the left and right subtrees separately
  node.left->accept(*this);
  std::vector<int> left_result = _result_;
  node.right->accept(*this);
  std::vector<int> right_result = _result_;

  _result_.clear();

  // Apply the binary operator (and/or)
  if (node.op == Operator::And) {
    std::sort(left_result.begin(), left_result.end());
    std::sort(right_result.begin(), right_result.end());
    std::set_intersection(left_result.begin(), left_result.end(),
                          right_result.begin(), right_result.end(),
                          std::back_inserter(_result_));
  } else if (node.op == Operator::Or) {
    std::sort(left_result.begin(), left_result.end());
    std::sort(right_result.begin(), right_result.end());
    std::set_union(left_result.begin(), left_result.end(), right_result.begin(),
                   right_result.end(), std::back_inserter(_result_));
  }
}

} // namespace lahuta
