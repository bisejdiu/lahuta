#ifndef LAHUTA_MAPPER_HPP
#define LAHUTA_MAPPER_HPP

#include "GraphMol/Atom.h"
#include "GraphMol/RWMol.h"
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "lahuta.hpp"
#include "matcher.hpp"
#include "nn.hpp"
#include "seq.hpp"
#include "topology.hpp"

namespace lahuta {

class TopologyMapper {
public:
  enum class MappingType { Query = 1, Target };

public:
  TopologyMapper(SeqData &sd_, MappingType type) : sd(sd_), type_(type) {
    luni_ptr = std::make_unique<Luni>(Luni::build(sd.st));
  }

  void map(const Matcher::result_t &res) {
    map_backtrace(res);
    map_atoms();
  }

  std::optional<unsigned int> get_mapped_resid(unsigned int atom_index) const {
    return mapped_atom_indices[atom_index];
  }

  void find_and_log_residue(unsigned int atom_id) {
    static_assert(std::is_same_v<unsigned int, decltype(atom_id)>);
    auto &mol = luni_ptr->get_molecule();
    auto atom = mol.getAtomWithIdx(atom_id);
    auto ri = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
    auto r_ix = ri->getResidueIndex();
    auto r_id = ri->getResidueNumber();

    auto mapped_resid = get_mapped_resid(atom_id);
    if (mapped_resid) {
      std::cout << "Residue index: " << r_id << " " << r_ix << " -> " //
                << *mapped_resid << " " << ri->getName()              //
                << "-" << ri->getResidueName() << "-" << ri->getResidueNumber() << std::endl;
    }
  }

  std::vector<const RDKit::Atom *> get_entity_atoms(const EntityID &entity) const {
    switch (get_entity_type(entity)) {
      case EntityType::Group: {
        const GroupEntity &group = luni_ptr->get_entity<GroupEntity>(entity);
        return group.atoms;
      }
      case EntityType::Atom: {
        const RDKit::Atom &atom = luni_ptr->get_entity<RDKit::Atom>(entity);
        return {&atom};
      }
      case EntityType::Ring: {
        const RingEntity &ring = luni_ptr->get_entity<RingEntity>(entity);
        return ring.atoms;
      }
    }
  }

  const Luni &get_luni() const { return *luni_ptr; };

private:
  void map_backtrace(const Matcher::result_t &res) {
    switch (type_) {
      case MappingType::Query:
        indices = get_indices(res, 'D');
        start_pos = res.qStartPos;
        end_pos = res.qEndPos;
        break;
      case MappingType::Target:
        indices = get_indices(res, 'I');
        start_pos = res.dbStartPos;
        end_pos = res.dbEndPos;
        break;
      default:
        throw std::invalid_argument("Unsupported MappingType");
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

struct EquivalencyConfig {
  bool contact_type = false;
  bool number_of_atoms = false;
  bool atom_name = false;
  bool element = false;
  bool resname = false;
};

class TopologicalEquivalency {
public:
  TopologicalEquivalency(
      const TopologyMapper &lm1, const TopologyMapper &lm2,
      std::optional<EquivalencyConfig> config = std::nullopt)
      : lm1_(lm1), lm2_(lm2), cfg(config.value_or(EquivalencyConfig())) {}

  bool evaluate(const RDKit::Atom *a1, const RDKit::Atom *a2) const {
    auto m_a1 = lm1_.get_mapped_resid(a1->getIdx());
    auto m_a2 = lm2_.get_mapped_resid(a2->getIdx());

    if (!m_a1 || !m_a2 || (m_a1.value() != m_a2.value())) return false;
    if (cfg.element && (a1->getSymbol() != a2->getSymbol())) return false;
    if (cfg.atom_name && compare_atom_names(a1, a2)) return false;
    if (cfg.resname && compare_resnames(a1, a2)) return false;

    return true;
  }

  bool evaluate(const std::vector<const RDKit::Atom *> a1, const std::vector<const RDKit::Atom *> a2) const {
    if (a1.empty() || a2.empty()) return false;
    if (cfg.number_of_atoms && (a1.size() != a2.size())) return false;
    if (!cfg.number_of_atoms) return evaluate(a1.front(), a2.front()); // FIX: check the mappable atom

    for (auto ix = 0; ix < a1.size(); ix++) {
      if (!evaluate(a1[ix], a2[ix])) return false;
    }

    return true;
  }

  bool should_consider(const std::vector<const RDKit::Atom *> &atoms, TopologyMapper::MappingType mt) const {
    if (atoms.empty()) return false;

    return cfg.number_of_atoms
               ? std::all_of(atoms.begin(), atoms.end(), [&](const auto *at) { return is_mappable(mt, at); })
               : std::any_of(atoms.begin(), atoms.end(), [&](const auto *at) { return is_mappable(mt, at); });
  }

  bool is_mappable(TopologyMapper::MappingType mt, const RDKit::Atom *atom) const {
    std::optional<unsigned int> val;
    switch (mt) {
      case TopologyMapper::MappingType::Query:
        val = lm1_.get_mapped_resid(atom->getIdx());
        break;
      case TopologyMapper::MappingType::Target:
        val = lm2_.get_mapped_resid(atom->getIdx());
        break;
      default:
        throw std::invalid_argument("Unsupported MappingType");
    }

    return val.has_value();
  }

  const TopologyMapper &get_lm(TopologyMapper::MappingType type) const {
    switch (type) {
      case TopologyMapper::MappingType::Query:
        return lm1_;
      case TopologyMapper::MappingType::Target:
        return lm2_;
      default:
        throw std::invalid_argument("Unsupported MappingType");
    }
  }

private:
  bool compare_atom_names(const RDKit::Atom *a1, const RDKit::Atom *a2) const noexcept {
    auto r1 = static_cast<const RDKit::AtomPDBResidueInfo *>(a1->getMonomerInfo());
    auto r2 = static_cast<const RDKit::AtomPDBResidueInfo *>(a2->getMonomerInfo());
    if (r1->getName() != r2->getName()) return false;
    return true;
  }

  bool compare_resnames(const RDKit::Atom *a1, const RDKit::Atom *a2) const noexcept {
    auto r1 = static_cast<const RDKit::AtomPDBResidueInfo *>(a1->getMonomerInfo());
    auto r2 = static_cast<const RDKit::AtomPDBResidueInfo *>(a2->getMonomerInfo());
    if (r1->getResidueName() != r2->getResidueName()) return false;
    return true;
  }

  const TopologyMapper &lm1_;
  const TopologyMapper &lm2_;
  EquivalencyConfig cfg;
};

class LahutaMapper {
  using Mapping = TopologyMapper::MappingType;
  using AtomVec = std::vector<RDKit::Atom *>;

public:
  LahutaMapper(SeqData &query, SeqData &target) : qm(query, Mapping::Query), tm(target, Mapping::Target) {}

  void map(const Matcher::result_t &res, std::optional<EquivalencyConfig> config = std::nullopt) {
    cfg = config.value_or(EquivalencyConfig());
    qm.map(res);
    tm.map(res);

    te = std::make_unique<TopologicalEquivalency>(qm, tm, cfg);
  }

  bool evaluate(const Contact &c1, Mapping m1, const Contact &c2, Mapping m2) const {
    if (m1 == m2) throw std::runtime_error("Mapping not sensible between same object type");
    if (!te) throw std::runtime_error("No alignment result has been mapped.");

    if (cfg.contact_type && c1.type != c2.type) return false;

    auto [c1_e1, c1_e2] = get_entity_atoms(c1, m1);
    auto [c2_e1, c2_e2] = get_entity_atoms(c2, m2);

    auto pair1 = te->evaluate(c1_e1, c2_e1);
    auto pair2 = te->evaluate(c1_e2, c2_e2);
    return pair1 && pair2;
  }

  void evaluate(const Contacts &c1, Mapping m1, const Contacts &c2, Mapping m2) const {
    if (m1 == m2) throw std::runtime_error("Mapping not sensible between same object type");
    if (!te) throw std::runtime_error("No alignment result has been mapped.");

    // TODO: Optimize! This is O(n^2). Interactions can be sorted, so we may be able to do better. The
    // should_consider checks help, but they only mitigate the issue. We can also avoid some overhead with
    // returning tuples by defining temporary vectors and using them to store atoms, but this, too, is not an
    // ideal solution. Ideally, we should be able to tell when we have exceeded the equivalence point in the
    // alignment and stop checking.
    for (const auto &i : c1.interactions) {
      auto [c1_e1, c1_e2] = get_entity_atoms(i, m1);
      if (!te->should_consider(c1_e1, m1) || !te->should_consider(c1_e2, m1)) continue;
      for (const auto &j : c2.interactions) {
        auto [c2_e1, c2_e2] = get_entity_atoms(j, m2);
        if (!te->should_consider(c2_e1, m2) || !te->should_consider(c2_e2, m2)) continue;

        if (cfg.contact_type && i.type != j.type) continue;
        auto pair1 = te->evaluate(c1_e1, c2_e1);
        auto pair2 = te->evaluate(c1_e2, c2_e2);

        if (pair1 && pair2) {
          log_contact(qm, i);
          log_contact(tm, j);
          std::cout << std::endl;
          break;
        }
      }
    }
  }

  void log_contact(const TopologyMapper &lm, const Contact c) const noexcept {
    std::string e1_atoms = stringify_atom_ids(lm.get_entity_atoms(c.entity1));
    std::string e2_atoms = stringify_atom_ids(lm.get_entity_atoms(c.entity2));

    std::cout << "Interaction between entity "                                              //
              << entity_type_to_string(get_entity_type(c.entity1)) << " ("                  //
              << get_entity_index(c.entity1) << ") and entity "                             //
              << entity_type_to_string(get_entity_type(c.entity2)) << " ("                  //
              << get_entity_index(c.entity2) << ") with distance " << c.distance << " --- " //
              << e1_atoms << " --- " << e2_atoms << "\n";
  }

  std::string stringify_atom_ids(const std::vector<const RDKit::Atom *> &atoms) const {
    std::string atom_ids;
    for (const auto *atom : atoms) {
      atom_ids += std::to_string(atom->getIdx()) + " ";
    }
    return atom_ids;
  }

  std::tuple<std::vector<const RDKit::Atom *>, std::vector<const RDKit::Atom *>>
  get_entity_atoms(const Contact &c, Mapping m) const {

    switch (m) {
      case Mapping::Query:
        return std::make_tuple(qm.get_entity_atoms(c.entity1), qm.get_entity_atoms(c.entity2));
      case Mapping::Target:
        return std::make_tuple(tm.get_entity_atoms(c.entity1), tm.get_entity_atoms(c.entity2));
      default:
        throw std::invalid_argument("Unsupported MappingType");
    }
  }

  const Luni &get_luni(Mapping type) const {
    switch (type) {
      case Mapping::Query:
        return qm.get_luni();
      case Mapping::Target:
        return tm.get_luni();
      default:
        throw std::invalid_argument("Unsupported MappingType");
    }
  }

private:
  TopologyMapper qm;
  TopologyMapper tm;
  EquivalencyConfig cfg;
  std::unique_ptr<TopologicalEquivalency> te;
};

} // namespace lahuta

#endif // LAHUTA_MAPPER_HPP
