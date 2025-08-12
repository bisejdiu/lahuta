#ifndef LAHUTA_HPP
#define LAHUTA_HPP

#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "convert.hpp"
#include "logging.hpp"
#include "topology.hpp"

constexpr const char *LAHUTA_VERSION = "0.75.0";

// clang-format off
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

  // FIX: no need to use optional here
  bool build_topology(std::optional<TopologyBuildingOptions> tops = std::nullopt); 

  std::string get_file_name() const { return file_name_; };
  const Topology &get_topology() const {
    if (!topology) {
      throw std::logic_error("Topology not built");
    }
    return *topology;
  }

  std::unique_ptr<Topology> release_topology() {
    if (!topology) {
      Logger::get_logger()->error("Topology not initialized. Cannot release topology.");
      return nullptr;
    }
    topology_built_ = false;
    return std::move(topology);
  }

  bool has_topology_built() const { return topology_built_; }

  /// filter the molecule based on the atom indices
  Luni filter(std::vector<int> &atom_indices) const;

  /// Can be called using the topology
  void assign_molstar_atom_types()  {
    if (topology) { topology->assign_molstar_typing(); }
    else { Logger::get_logger()->error("Topology not initialized. Cannot assign Molstar atom types."); }
  }

  void assign_arpeggio_atom_types() {
    if (topology) { topology->assign_arpeggio_atom_types(); }
    else { Logger::get_logger()->error("Topology not initialized. Cannot assign Arpeggio atom types."); }
  }

  /// Enable or disable a specific computation in the topology
  void enable_topology_computation(TopologyComputation comp, bool enabled) {
    ensure_topology_initialized();
    if (topology) {
      topology->enable_computation(comp, enabled);
    }
  }

  /// Enable only the specified computations (disabling all others)
  void enable_topology_only(TopologyComputation comps) {
    ensure_topology_initialized();
    if (topology) {
      topology->enable_only(comps);
    }
  }

  /// Check if a specific computation is enabled
  bool is_topology_computation_enabled(TopologyComputation comp) const {
    if (topology) {
      return topology->is_computation_enabled(comp);
    }
    Logger::get_logger()->error("Topology not initialized. Cannot check computation status.");
    return false;
  }

  /// Execute a specific computation with its dependencies
  bool execute_topology_computation(TopologyComputation comp) {
    if (topology) {
      return topology->execute_computation(comp); 
    }
    Logger::get_logger()->error("Topology not initialized. Cannot execute computation.");
    return false;
  }

  /// Set the cutoff for neighbor search
  void set_neighbor_search_cutoff(double cutoff) {
    ensure_topology_initialized();
    if (topology) {
      topology->set_cutoff(cutoff);
    }
  }

  /// Set the atom typing method
  void set_atom_typing_method(ContactComputerType method) {
    ensure_topology_initialized();
    if (topology) {
      topology->set_atom_typing_method(method);
    }
  }

  const auto n_atoms() const { return mol->getNumAtoms(); }
  const std::vector<int> indices()    const;
  const std::vector<int> resids()     const;
  const std::vector<int> resindices() const;
  const std::vector<int> atomic_numbers() const;
  const std::vector<std::string> names()    const;
  const std::vector<std::string> symbols()  const;
  const std::vector<std::string> elements() const;
  const std::vector<std::string> resnames() const;
  const std::vector<std::string> chainlabels() const;

  const auto &get_molecule() const { return *mol; }
  const auto &get_conformer(int id = -1) const { return mol->getConformer(id); }
  const auto *get_atom(int idx) const { return mol->getAtomWithIdx(idx); }

  const auto *get_info(int idx) const {
    return static_cast<const RDKit::AtomPDBResidueInfo *>(get_atom(idx)->getMonomerInfo());
  }

  const std::vector<RDGeom::Point3D> &positions(int confId = -1) const {
    return get_conformer(confId).getPositions();
  }

  friend class Contacts;

  /// very rough estimate of the memory size
  size_t total_size() const;

private:
  explicit Luni(std::shared_ptr<RDKit::RWMol> valid_mol) 
    : mol(valid_mol), topology(std::make_unique<Topology>(valid_mol)), topology_built_(false) {}

  void ensure_topology_initialized() {
    if (!topology) {
      Logger::get_logger()->debug("Initializing topology for configuration");
      topology = std::make_unique<Topology>(mol);
    }
  }

  auto match_smarts_string(std::string sm, std::string atype = "", bool log_values = false) const;

  template <typename T>
  std::vector<T> atom_attrs(std::function<T(const RDKit::Atom *)> func) const;

  template <typename T>
  std::vector<std::reference_wrapper<const T>>
  atom_attrs_ref(std::function<const T &(const RDKit::Atom *)> func) const;

  // FIX: It should not be necessary to have a default constructed RWMol here
  std::shared_ptr<RDKit::RWMol> mol = std::make_shared<RDKit::RWMol>();
  std::unique_ptr<Topology> topology;
  bool topology_built_ = false;

  std::string file_name_;
  std::vector<int> filtered_indices;
  bool is_in_filtered_state = false; // ambitious name, for know it's just a flag
};

} // namespace lahuta

#endif // LAHUTA_HPP
