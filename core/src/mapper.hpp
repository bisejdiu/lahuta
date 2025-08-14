#ifndef LAHUTA_MAPPER_HPP
#define LAHUTA_MAPPER_HPP

#include "GraphMol/Atom.h"
#include <entities/records.hpp>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "lahuta.hpp"
#include "logging.hpp"
#include "matcher.hpp"
#include "seq.hpp"
#include "topology.hpp"
#include "backtrace.hpp"
#include "entities/contact.hpp"

#include <algorithm>
#include <memory>
#include <optional>
#include <unordered_map>
#include <vector>
// get the custom has function for pairs
#include <functional>
#include <utility>

namespace lahuta {

constexpr std::size_t GOLDEN_RATIO = 0x9e3779b97f4a7c15ULL;

class TopologyMapper {
public:
  enum class MappingType { Query = 1, Target };

public:
  TopologyMapper(SeqData &sd_, MappingType type) : sd(sd_), type_(type) {
    luni_ptr = std::make_unique<Luni>(Luni::create(sd.st));
    auto ok = luni_ptr->build_topology();
    if (!ok) {
      Logger::get_logger()->warn("Failed building the topology!");
    }
  }

  void map(const Matcher::result_t &res) { map_backtrace_and_atoms(res); }

  std::optional<unsigned int> get_mapped_resid(unsigned int atom_index) const { return _map_[atom_index]; }

  // FIX: get_entity_atoms does not belong here
  std::vector<const RDKit::Atom *> get_entity_atoms(const EntityID &entity) const {
    switch (entity.kind()) {
      case Kind::Group: {
        const auto &group = luni_ptr->get_topology().records<GroupRec>()[entity.index()];
        auto atom_indices = group.atoms;
        std::vector<const RDKit::Atom *> atoms;
        atoms.reserve(atom_indices.size());
        for (const auto atom_index : atom_indices) {
          // atoms.push_back(luni_ptr->get_molecule().getAtomWithIdx(atom_index));
          atoms.push_back(&atom_index.get());
        }
        return atoms;
      }
      case Kind::Atom: {
        const auto atom_idx = luni_ptr->get_topology().records<AtomRec>()[entity.index()].atom.get().getIdx();
        return {luni_ptr->get_molecule().getAtomWithIdx(atom_idx)};
      }
      case Kind::Ring: {
        const auto &ring = luni_ptr->get_topology().records<RingRec>()[entity.index()];
        std::vector<const RDKit::Atom *> atoms;
        atoms.reserve(ring.atoms.size());
        for (const auto atom_index : ring.atoms) {
          // atoms.push_back(luni_ptr->get_molecule().getAtomWithIdx(atom_index));
          atoms.push_back(&atom_index.get());
        }
        return atoms;
      }
    }
  }

  const Luni &get_luni() const { return *luni_ptr; };

private:
  void map_backtrace_and_atoms(const Matcher::result_t &res) {
    auto &mol = luni_ptr->get_molecule();
    _map_.assign(mol.getNumAtoms(), std::nullopt);

    int start = (type_ == MappingType::Query) ? res.qStartPos : res.dbStartPos;
    int end   = (type_ == MappingType::Query) ? res.qEndPos   : res.dbEndPos;

    BacktraceParser parser{res.backtrace};

    int residue_index = 0;
    for (const auto &residue : luni_ptr->get_topology().get_residues()) {
      if (residue.chain_id != sd.chain_name) continue;

      if (residue_index < start || residue_index > end) {
        residue_index++;
        continue;
      }

      bool residue_is_mapped = false;
      char letter = parser.current();

      if (letter == 'M') {
        residue_is_mapped = true;
      } else if (type_ == MappingType::Query && letter == 'D') {
        parser.skip('D');
        if (parser.has_next() && parser.current() == 'M') residue_is_mapped = true;
      } else if (type_ == MappingType::Target && letter == 'I') {
        parser.skip('I');
        if (parser.has_next() && parser.current() == 'M') residue_is_mapped = true;
      }

      if (residue_is_mapped) {
        for (const auto atom : residue.atoms) {
          _map_[atom->getIdx()] = parser.position();
        }
      }

      residue_index++;
      parser.get_next();
    }
  }

  std::unique_ptr<Luni> luni_ptr;
  const SeqData &sd;
  MappingType type_;
  std::vector<std::optional<int>> _map_;
};

struct ContactEquivKey {
  // sorted (ascending) list of mapped residue IDs
  std::vector<unsigned int> e1_mapped_ids;
  std::vector<unsigned int> e2_mapped_ids;

  InteractionType contact_type = InteractionType::Generic;
};

struct ContactEquivKeyHash {
  std::size_t operator()(const ContactEquivKey &key) const {

    std::size_t h = std::hash<int>()(static_cast<int>(key.contact_type));

    for (auto val : key.e1_mapped_ids) {
      h ^= std::hash<unsigned int>()(val) + GOLDEN_RATIO + (h << 6) + (h >> 2);
    }

    for (auto val : key.e2_mapped_ids) {
      h ^= std::hash<unsigned int>()(val) + GOLDEN_RATIO + (h << 6) + (h >> 2);
    }

    return h;
  }
};

struct ContactEquivKeyEqual {
  bool operator()(const ContactEquivKey &lhs, const ContactEquivKey &rhs) const {
    // FIX: here we force the contact type to be the same
    if (lhs.contact_type != rhs.contact_type)                 return false;
    if (lhs.e1_mapped_ids.size() != rhs.e1_mapped_ids.size()) return false;
    if (lhs.e2_mapped_ids.size() != rhs.e2_mapped_ids.size()) return false;
    if (lhs.e1_mapped_ids != rhs.e1_mapped_ids)               return false;
    if (lhs.e2_mapped_ids != rhs.e2_mapped_ids)               return false;
    return true;
  }
};

struct EquivalencyConfig {
  bool contact_resolution = false; // Same resolution
  bool contact_type = true;        // Same interaction type (takes precedence over hbond_type)
  bool hbond_type = false;         // Same hydrogen bond type
  bool number_of_atoms = false;    // Number of atoms must be the same
  bool atom_name = false;          // Same atom name
  bool element = false;            // Same element
  bool resname = false;            // Same residue name
  //

  bool one_contact_per_residue = false; // Only one contact per residue // FIX: not used
  // NOTE: I think we should either have atom_name set to true or one_contact_per_residue.
  // if we allow both to be false, then we'd get mappings for all pairwise contacts between residues.
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
    if (!cfg.number_of_atoms)
      return evaluate(a1.front(), a2.front()); // FIX: not first, but get the mappable atom

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

public:
  LahutaMapper(SeqData &query, SeqData &target) : qm(query, Mapping::Query), tm(target, Mapping::Target) {}

  void map(const Matcher::result_t &res, std::optional<EquivalencyConfig> config = std::nullopt) {
    cfg = config.value_or(EquivalencyConfig());
    qm.map(res);
    tm.map(res);

    te = std::make_unique<TopologicalEquivalency>(qm, tm, cfg);
  }

  int evaluate(const ContactSet &c1, Mapping m1, const ContactSet &c2, Mapping m2) {
    if (m1 == m2) throw std::runtime_error("Mapping not sensible between same object type");
    // FIX: could also make an initialization parameter
    if (!te) throw std::runtime_error("No alignment result has been mapped.");

    // build an index for c1
    std::unordered_map<ContactEquivKey, std::vector<const Contact *>, ContactEquivKeyHash, ContactEquivKeyEqual> index;

    index.reserve(c1.data().size());
    for (auto &contact1 : c1.data()) {
      ContactEquivKey key = make_equiv_key(contact1, m1);
      index[key].push_back(&contact1);
    }

    int count = 0;
    for (auto &contact2 : c2.data()) {

      auto [e1_atoms, e2_atoms] = get_entity_atoms(contact2, m2);

      if (handle_contact_type_skippable(contact2)) {
        continue;
      }

      ContactEquivKey key2 = make_equiv_key(contact2, m2);

      // Lookup in the index
      auto it = index.find(key2);
      if (it == index.end()) continue;

      for (auto *c1_candidate : it->second) {
        if (cfg.contact_resolution && !is_same_resolution(*c1_candidate, contact2)) continue;

        if (handle_contact_type(*c1_candidate, contact2)) {
          continue;
        }

        // FIX: finalize all checks here
        //
        // c1_candidate is from c1, contact2 is from c2
        auto [c1_e1, c1_e2] = get_entity_atoms(*c1_candidate, m1);
        auto [c2_e1, c2_e2] = get_entity_atoms(contact2, m2);

        bool pair1 = te->evaluate(c1_e1, c2_e1);
        bool pair2 = te->evaluate(c1_e2, c2_e2);

        if (pair1 && pair2) {

          auto c1_e1_a = c1_e1.front();
          auto c1_e2_a = c1_e2.front();
          auto c2_e1_a = c2_e1.front();
          auto c2_e2_a = c2_e2.front();

          auto m_c1_e1 = te->get_lm(Mapping::Query). get_mapped_resid(c1_e1.front()->getIdx());
          auto m_c1_e2 = te->get_lm(Mapping::Query). get_mapped_resid(c1_e2.front()->getIdx());
          auto m_c2_e1 = te->get_lm(Mapping::Target).get_mapped_resid(c2_e1.front()->getIdx());
          auto m_c2_e2 = te->get_lm(Mapping::Target).get_mapped_resid(c2_e2.front()->getIdx());

          auto r1_a1 = static_cast<const RDKit::AtomPDBResidueInfo *>(c1_e1.front()->getMonomerInfo());
          auto r1_a2 = static_cast<const RDKit::AtomPDBResidueInfo *>(c1_e2.front()->getMonomerInfo());
          auto r2_a1 = static_cast<const RDKit::AtomPDBResidueInfo *>(c2_e1.front()->getMonomerInfo());
          auto r2_a2 = static_cast<const RDKit::AtomPDBResidueInfo *>(c2_e2.front()->getMonomerInfo());

          /*if (r1_a1->getName() != r2_a1->getName() || r1_a2->getName() != r2_a2->getName()) {*/
          /*  continue;*/
          /*}*/

          std::cout << "info: "
                    << r1_a1->getResidueName() << " - " << r1_a2->getResidueName() << " "
                    << r2_a1->getResidueName() << " - " << r2_a2->getResidueName() << " : "
                    << std::setw(3) << r1_a1->getName() << " - " << std::setw(3) << r1_a2->getName() << " "
                    << std::setw(3) << r2_a1->getName() << " - " << std::setw(3) << r2_a2->getName() << " : "
                    << c1_e1_a->getIdx() << " - " << c1_e2_a->getIdx() << " "
                    << c2_e1_a->getIdx() << " - " << c2_e2_a->getIdx() << " : "
                    << m_c1_e1.value() << " - " << m_c1_e2.value() << " "
                    << m_c2_e1.value() << " - " << m_c2_e2.value() << " : "
                    << (int)c1_candidate->type << " - " << (int)contact2.type << " : "
                    << std::endl;

          count++;
          break; // break here means we consider only one contact per residue
        }
      }
    }

    return count;
  }

private:
  bool handle_contact_type_skippable(const Contact &c) const {
    if (!cfg.hbond_type
        && (c.type == InteractionType::HydrogenBond || c.type == InteractionType::WeakHydrogenBond)) {
      return false;
    }

    return false;
  }

  ContactEquivKey make_equiv_key(const Contact &c, Mapping m) const {
    ContactEquivKey key;

    auto [atoms1, atoms2] = get_entity_atoms(c, m);

    auto e1_ids = gather_mapped_ids(atoms1, m);
    auto e2_ids = gather_mapped_ids(atoms2, m);

    if (e1_ids.size() > 1) std::sort(e1_ids.begin(), e1_ids.end());
    if (e2_ids.size() > 1) std::sort(e2_ids.begin(), e2_ids.end());

    key.e1_mapped_ids = std::move(e1_ids);
    key.e2_mapped_ids = std::move(e2_ids);
    key.contact_type = c.type;

    return key;
  }

  std::vector<unsigned int>
  gather_mapped_ids(const std::vector<const RDKit::Atom *> &atoms, Mapping m) const {
    std::vector<unsigned int> ids;
    ids.reserve(atoms.size());

    for (auto *atom : atoms) {
      std::optional<unsigned int> resid;
      if (m == Mapping::Query) {
        resid = qm.get_mapped_resid(atom->getIdx());
      } else {
        resid = tm.get_mapped_resid(atom->getIdx());
      }
      if (resid.has_value()) {
        ids.push_back(resid.value());
      }
    }
    return ids;
  }

  bool handle_contact_type(const Contact &i, const Contact &j) {
    static auto is_hbond_type = [](const InteractionType type) {
      return type == InteractionType::HydrogenBond || type == InteractionType::WeakHydrogenBond;
    };

    if (cfg.hbond_type && is_hbond_type(i.type) && is_hbond_type(j.type)) {
      return i.type != j.type; // must match exactly
    }
    if (!cfg.hbond_type && is_hbond_type(i.type) && is_hbond_type(j.type)) {
      return false; // ignoring the difference in hbond
    }
    if (cfg.contact_type) {
      return i.type != j.type;
    }

    return false;
  }

  bool is_same_resolution(const Contact &c1, const Contact &c2) const {
    if (c1.lhs.kind() != c2.lhs.kind()) return false;
    if (c1.rhs.kind() != c2.rhs.kind()) return false;
    return true;
  }

  // Provide access to the relevant atoms for a contact
  std::tuple<std::vector<const RDKit::Atom *>, std::vector<const RDKit::Atom *>>
  get_entity_atoms(const Contact &c, Mapping m) const {
    switch (m) {
      case Mapping::Query:
        return std::make_tuple(qm.get_entity_atoms(c.lhs), qm.get_entity_atoms(c.rhs));
      case Mapping::Target:
        return std::make_tuple(tm.get_entity_atoms(c.lhs), tm.get_entity_atoms(c.rhs));
      default:
        throw std::invalid_argument("Unsupported MappingType");
    }
  }

public:
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
