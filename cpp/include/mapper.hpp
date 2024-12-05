#ifndef LAHUTA_MAPPER_HPP
#define LAHUTA_MAPPER_HPP

#include "GraphMol/RWMol.h"
#include <iostream>
#include <string>
#include <vector>

#include "matcher.hpp"
#include "seq.hpp"

namespace lahuta {

class Mapper {
public:
  Mapper(const SeqData &query_, const SeqData &target_, const Matcher::result_t &res_)
      : query(query_), target(target_), res(res_) {}

  void map() {
    query_indices = get_indices('D');
    target_indices = get_indices('I');
  }

  void print() {
    std::cout << "MQ: " << query.mol->getNumAtoms() << " " << query_indices.size() << "\n";
    std::cout << "MT: " << target.mol->getNumAtoms() << " " << target_indices.size() << "\n";
    std::cout << target_indices[21] << " " << target_indices[22] << " " << target_indices[23] << "\n";
    std::cout << query_indices[31] << " " << query_indices[32] << " " << query_indices[33] << "\n";
  }

  std::optional<unsigned int> get_query_index(unsigned int i) const {
    if (i < res.qStartPos || i >= res.qEndPos) {
      return std::nullopt;
    }
    auto q_idx = i - res.qStartPos;
    return query_indices[q_idx];
  };

  std::optional<unsigned int> get_target_index(unsigned int i) const {
    if (i < res.dbStartPos || i >= res.dbEndPos) {
      return std::nullopt;
    }
    auto t_idx = i - res.dbStartPos;
    return target_indices[t_idx];
  };

private:
  std::vector<size_t> get_indices(const char exclude) const {
    std::vector<size_t> indices;
    for (size_t i = 0; i < res.backtrace.size(); ++i) {
      char c = res.backtrace[i];
      if (c != exclude) {
        indices.push_back(i);
      }
    }
    return indices;
  }


public:
  std::shared_ptr<RDKit::RWMol> mol = std::make_shared<RDKit::RWMol>();
  SeqData query;
  SeqData target;
  const Matcher::result_t &res;
  std::vector<size_t> query_indices;
  std::vector<size_t> target_indices;
};

} // namespace lahuta

#endif // LAHUTA_MAPPER_HPP
