#ifndef LAHUTA_MAPPER_HPP
#define LAHUTA_MAPPER_HPP

#include "GraphMol/RWMol.h"
#include <iostream>
#include <string>
#include <vector>

#include "lahuta.hpp"
#include "matcher.hpp"
#include "seq.hpp"
#include "topology.hpp"

namespace lahuta {

inline std::vector<unsigned int>
factorize(const std::vector<unsigned int> &input) { // , std::true_type /*sorted*/) {
  std::vector<unsigned int> output;
  output.reserve(input.size());

  if (input.empty()) return output;

  int current_value = input[0];
  int current_code = 0;
  for (const auto &val : input) {
    if (val != current_value) {
      // Found a new unique value, increment code
      current_value = val;
      current_code++;
    }
    output.push_back(current_code);
  }

  return output;
}

class LuniMapper {
public:
  enum class MappingType { Query, Target };

public:
  LuniMapper(SeqData &sd_, MappingType type) : sd(sd_), type_(type) {
    luni_ptr = std::make_unique<Luni>(Luni::build(sd.st));
  }

  void map(const Matcher::result_t &res) {
    map_backtrace(res);
    map_atoms();
  }

  std::optional<unsigned int> get_mapped_resid(unsigned int atom_index) const {
    return mapped_atom_indices[atom_index];
  }

  const Luni &get_luni() { return *luni_ptr; };

private:
  void map_backtrace(const Matcher::result_t &res) {
    switch (type_) {
      case MappingType::Query:
        indices = get_query_indices(res);
        start_pos = res.qStartPos;
        end_pos = res.qEndPos;
        break;
      case MappingType::Target:
        indices = get_target_indices(res);
        start_pos = res.dbStartPos;
        end_pos = res.dbEndPos;
        break;
    }
  }

  void map_atoms() {
    const Topology &top = luni_ptr->get_topology();
    mapped_atom_indices.resize(luni_ptr->n_atoms());
    int resid_idx_distance = 0; // we need the resid index, not the resid
    const Residues &residues = *top.residues;
    for (const auto &residue : residues) {
      if (residue.chain_id != sd.chain_name) continue;
      for (const auto &atom : residue.atoms) {
        mapped_atom_indices[atom->getIdx()] = get_aligned_resid(resid_idx_distance);
      }
      resid_idx_distance++;
    }
  }

  std::vector<unsigned int> get_indices(const Matcher::result_t &res, const char exclude) const {
    std::vector<unsigned int> indices;
    for (size_t i = 0; i < res.backtrace.size(); ++i) {
      char c = res.backtrace[i];
      if (c != exclude) {
        indices.push_back(i);
      }
    }
    return indices;
  }

  std::vector<unsigned int> get_query_indices(const Matcher::result_t &res) const {
    return get_indices(res, 'D');
  }
  std::vector<unsigned int> get_target_indices(const Matcher::result_t &res) const {
    return get_indices(res, 'I');
  }

  std::optional<unsigned int> get_aligned_resid(unsigned int i) const {
    if (i < start_pos || i >= end_pos) {
      return std::nullopt;
    }
    auto idx = i - start_pos;
    return indices[idx];
  }

  std::unique_ptr<Luni> luni_ptr;
  const SeqData &sd;
  MappingType type_;

  std::vector<std::optional<unsigned int>> mapped_atom_indices;
  std::vector<unsigned int> indices;
  int start_pos = 0;
  int end_pos = 0;
};

class Mapper {
public:
  Mapper(const SeqData &query_, const SeqData &target_) // , const Matcher::result_t &res_)
      : query(query_), target(target_) /* ,  res(res_) */ {}

  void map(const Matcher::result_t &res_) {
    res = std::make_unique<Matcher::result_t>(res_);
    query_indices = get_query_indices();
    target_indices = get_target_indices();
  }

  void print() {
    std::cout << "MQ: " << query.mol->getNumAtoms() << " " << query_indices.size() << "\n";
    std::cout << "MT: " << target.mol->getNumAtoms() << " " << target_indices.size() << "\n";
    std::cout << target_indices[21] << " " << target_indices[22] << " " << target_indices[23] << "\n";
    std::cout << query_indices[31] << " " << query_indices[32] << " " << query_indices[33] << "\n";
  }

  std::optional<unsigned int> get_query_index(unsigned int i) const {
    if (!res) return std::nullopt;
    if (i < res->qStartPos || i >= res->qEndPos) {
      return std::nullopt;
    }
    auto q_idx = i - res->qStartPos;
    return query_indices[q_idx];
  };

  std::optional<unsigned int> get_target_index(unsigned int i) const {
    if (!res) return std::nullopt;
    if (i < res->dbStartPos || i >= res->dbEndPos) {
      return std::nullopt;
    }
    auto t_idx = i - res->dbStartPos;
    return target_indices[t_idx];
  };

private:
  std::vector<size_t> get_indices(const char exclude) const {
    std::vector<size_t> indices;
    for (size_t i = 0; i < res->backtrace.size(); ++i) {
      char c = res->backtrace[i];
      if (c != exclude) {
        indices.push_back(i);
      }
    }
    return indices;
  }

  std::vector<size_t> get_query_indices() const { return get_indices('D'); }

  std::vector<size_t> get_target_indices() const { return get_indices('I'); }

public:
  std::shared_ptr<RDKit::RWMol> mol = std::make_shared<RDKit::RWMol>();
  const SeqData &query;
  const SeqData &target;
  /*const Matcher::result_t &res;*/
  std::unique_ptr<Matcher::result_t> res;
  std::vector<size_t> query_indices;
  std::vector<size_t> target_indices;
};

} // namespace lahuta

#endif // LAHUTA_MAPPER_HPP
