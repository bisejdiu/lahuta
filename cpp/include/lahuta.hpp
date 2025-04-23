#ifndef LAHUTA_HPP
#define LAHUTA_HPP

#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <gemmi/mmread_gz.hpp>
#include <gemmi/model.hpp>
#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/BondIterators.h>

#include "convert.hpp"
#include "logging.hpp"
#include "topology.hpp"

constexpr const char* LAHUTA_VERSION = "0.50.0";

namespace lahuta {

// NOTE: rename to Lahuta?
class Luni {
public:
  Luni(const Luni&)            = delete;
  Luni& operator=(const Luni&) = delete;

  Luni(Luni&&)            = default;
  Luni& operator=(Luni&&) = default;

  explicit Luni(std::string file_name);
  explicit Luni(std::string file_name, bool test);

  static Luni create(const IR &ir);
  static Luni create(const gemmi::Structure &st);
  static Luni create(std::shared_ptr<RDKit::RWMol> mol) { return Luni(mol); }

  bool build_topology(std::optional<TopologyBuildingOptions> tops = std::nullopt); 

  std::string get_file_name() const { return file_name_; };
  const Topology &get_topology() const { return *get_topology_ptr(); }
  bool has_topology_built() { return topology.has_value(); };

  RDKit::RWMol &get_molecule() { return *mol; }
  const RDKit::RWMol &get_molecule() const { return *mol; }

  auto *get_info(int idx) { return static_cast<RDKit::AtomPDBResidueInfo *>(get_atom(idx)->getMonomerInfo()); }
  const auto *get_info(int idx) const { return static_cast<const RDKit::AtomPDBResidueInfo *>(get_atom(idx)->getMonomerInfo()); }

  RDKit::Conformer &get_conformer(int id = -1) { return mol->getConformer(id); }
  const RDKit::Conformer &get_conformer(int id = -1) const { return mol->getConformer(id); }

  const AtomEntityCollection  &get_atom_types() const { return get_topology_ptr()->get_atom_types(); }
  const RingEntityCollection  &get_rings()      const { return get_topology_ptr()->get_rings(); }
  const GroupEntityCollection &get_features()   const { return get_topology_ptr()->get_features(); }

  /// filter the molecule based on the atom indices
  Luni filter(std::vector<int> &atom_indices) const;

  /// EntityID -> AtomEntity/GroupEntity/RingEntity
  template <typename T> const T &get_entity(EntityID id) const;

  /// AtomEntity/GroupEntity/RingEntity -> EntityID
  const std::vector<EntityID> &get_or_create_atom_entities();
  const std::vector<EntityID> &get_or_create_ring_entities();
  const std::vector<EntityID> &get_or_create_group_entities();

  /// Can be called using the topology
  void assign_molstar_atom_types()  { 
    if (topology) { topology->assign_molstar_typing(); } 
    else { Logger::get_logger()->error("Topology not built. Cannot assign Molstar atom types."); }
  }

  void assign_arpeggio_atom_types() {
    if (topology) { topology->assign_arpeggio_atom_types(); } 
    else { Logger::get_logger()->error("Topology not built. Cannot assign Arpeggio atom types."); }
  }

  //! Returns the atoms of the molecule.
  const auto atoms() const { return mol->atoms(); }

  //! Returns the number of atoms in the molecule.
  const auto n_atoms() const { return mol->getNumAtoms(); }

  //! Returns the names of the atoms.
  const std::vector<std::string> names() const;

  //! Returns the symbols of the atoms.
  const std::vector<std::string> symbols() const;

  //! Returns the residue indices of the atoms.
  const std::vector<int> indices() const;

  //! Returns the atomic numbers of the atoms.
  const std::vector<int> atomic_numbers() const;

  //! Returns the elements of the atoms.
  const std::vector<std::string> elements() const;

  //! Returns the residue names of the atoms.
  const std::vector<std::string> resnames() const;

  //! Returns the residue ids of the atoms.
  const std::vector<int> resids() const;

  //! Returns the residue indices of the atoms.
  const std::vector<int> resindices() const;

  //! Returns the chain labels of the atoms.
  const std::vector<std::string> chainlabels() const;

  const std::vector<RDGeom::Point3D> &positions(int confId = -1) const {
    return get_conformer(confId).getPositions();
  }

  friend class Contacts;

  // FIX: add helper functions to get topology information
  const RDKit::Atom *get_atom(int idx) const { return mol->getAtomWithIdx(idx); }
  RDKit::Atom *get_atom(int idx) { return mol->getAtomWithIdx(idx); }

  /// very rough estimate of the memory size
  size_t total_size() const;

private:
  explicit Luni(std::shared_ptr<RDKit::RWMol> valid_mol) : mol(valid_mol), topology(valid_mol) {}

  auto match_smarts_string(std::string sm, std::string atype = "", bool log_values = false) const;
  const Topology* get_topology_ptr() const;
  void read_structure();

  template <typename T>
  std::vector<T> atom_attrs(std::function<T(const RDKit::Atom *)> func) const;

  template <typename T>
  std::vector<std::reference_wrapper<const T>>
  atom_attrs_ref(std::function<const T &(const RDKit::Atom *)> func) const;


  std::shared_ptr<RDKit::RWMol> mol = std::make_shared<RDKit::RWMol>();
  std::optional<Topology> topology;
  std::unordered_map<EntityType, std::vector<EntityID>> entities;

  std::string file_name_;
  std::vector<int> filtered_indices;
  bool is_in_filtered_state = false; // ambitious name, for know it's just a flag
};

} // namespace lahuta

#endif // LAHUTA_HPP
