#ifndef LAHUTA_MAPPER_HPP
#define LAHUTA_MAPPER_HPP

#include "GraphMol/RWMol.h"
#include <iostream>
#include <string>
#include <vector>

#include "extract.hpp"
#include "ob/bitvec.h"

namespace lahuta {

inline std::vector<size_t> get_non_deletion_indices(const std::string &cigar) {
  std::vector<size_t> indices;
  size_t current_pos = 0;
  size_t count = 0;

  for (size_t i = 0; i < cigar.size(); ++i) {
    char c = cigar[i];
    if (c >= '0' && c <= '9') {
      count = count * 10 + c - '0';
    } else {
      if (c != 'D') {
        size_t repeats = (count == 0) ? 1 : count;
        for (size_t j = 0; j < repeats; ++j) {
          indices.push_back(current_pos + j);
        }
        current_pos += repeats;
      } else {
        current_pos += (count == 0 ? 1 : count);
      }
      count = 0;
    }
  }
  return indices;
}

inline std::vector<size_t> get_non_insertion_indices(const std::string &cigar) {
  std::vector<size_t> indices;
  size_t current_pos = 0;
  size_t count = 0;

  for (size_t i = 0; i < cigar.size(); ++i) {
    char c = cigar[i];
    if (c >= '0' && c <= '9') {
      count = count * 10 + c - '0';
    } else {
      if (c != 'I') {
        size_t repeats = (count == 0) ? 1 : count;
        for (size_t j = 0; j < repeats; ++j) {
          indices.push_back(current_pos + j);
        }
        current_pos += repeats;
      } else {
        current_pos += (count == 0 ? 1 : count);
      }
      count = 0;
    }
  }
  return indices;
}

class Mapper {
public:
  Mapper(const SeqData &query_, const SeqData &target_, const Matcher::result_t &res_)
      : query(query_), target(target_), res(res_) {}

  void map() {
    query_indices = get_non_deletion_indices(res.backtrace);
    target_indices = get_non_insertion_indices(res.backtrace);
  }

  void print() {
    std::cout << "MQ: " << query.mol->getNumAtoms() << " " << query_indices.size() << "\n";
    std::cout << "MT: " << target.mol->getNumAtoms() << " " << target_indices.size() << "\n";
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
