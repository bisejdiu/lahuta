#ifndef LAHUTA_HPP
#define LAHUTA_HPP

#include <atomic>
#include <condition_variable>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <vector>

#include <gemmi/model.hpp>
#include <rdkit/Geometry/point.h>

#include "logging.hpp"
#include "topology.hpp"

// clang-format off
namespace lahuta {

struct ModelParserResult;
struct IR;
enum class ModelTopologyMethod;

// NOTE: rename to Lahuta?
class Luni {
public:
  Luni(const Luni&)            = delete;
  Luni& operator=(const Luni&) = delete;

  Luni(Luni&&)            noexcept;
  Luni& operator=(Luni&&) noexcept;

  explicit Luni(std::string file_name);

  // Tag type for model-file input path
  struct ModelFileTag { explicit ModelFileTag() = default; };
  static inline constexpr ModelFileTag ModelFile{};
  explicit Luni(std::string file_name, ModelFileTag);

  static Luni create(const IR &ir);
  static Luni create(const gemmi::Structure &st);
  static Luni create(std::shared_ptr<RDKit::RWMol> mol) { return Luni(mol); }
  static Luni create(std::shared_ptr<RDKit::RWMol> mol, TopologyBuildMode mode) {
    Luni l(std::move(mol));
    if (mode == TopologyBuildMode::Model) l.model_origin_ = true;
    return l;
  }
  static Luni from_model_file(std::string file_name) { return Luni(std::move(file_name), ModelFile); }
  static Luni from_model_data(const ModelParserResult &data);
  static Luni from_model_data(const ModelParserResult &data, ModelTopologyMethod method);

  // Concurrency: serialized per Luni instance. If multiple threads call concurrently
  // with different TopologyBuildingOptions, the first successful one wins
  bool build_topology(std::optional<TopologyBuildingOptions> tops = std::nullopt) const;

  std::string get_file_name() const { return file_name_; };

  std::shared_ptr<const Topology> get_topology() const {
    if (!topology) {
      Logger::get_logger()->error("Topology not initialized. Cannot get shared topology.");
      return nullptr;
    }
    return topology;
  }

  std::shared_ptr<Topology> get_topology() {
    if (!topology) {
      Logger::get_logger()->error("Topology not initialized. Cannot get shared topology.");
      return nullptr;
    }
    return topology;
  }

  bool has_topology_built() const { return topology_built_.load(std::memory_order_acquire); }

  /// filter the molecule based on the atom indices
  Luni filter(std::vector<int> &atom_indices) const;

  /// Enable or disable a specific computation in the topology
  void enable_computation(TopologyComputation comp, bool enabled) const {
    ensure_topology_initialized();
    if (topology) {
      topology->enable_computation(comp, enabled);
    }
  }

  /// Enable only the specified computations (disabling all others)
  void enable_only(TopologyComputation comps) const {
    ensure_topology_initialized();
    if (topology) {
      topology->enable_only(comps);
    }
  }

  /// Check if a specific computation is enabled
  bool is_computation_enabled(TopologyComputation comp) const {
    if (topology) {
      return topology->is_computation_enabled(comp);
    }
    Logger::get_logger()->error("Topology not initialized. Cannot check computation status.");
    return false;
  }

  // Whether this system originated from a model input
  bool is_model_origin() const { return model_origin_; }

  /// Execute a specific computation with its dependencies
  bool execute_computation(TopologyComputation comp) {
    if (topology) {
      return topology->execute_computation(comp); 
    }
    Logger::get_logger()->error("Topology not initialized. Cannot execute computation.");
    return false;
  }

  /// Set the cutoff for neighbor search
  void set_search_cutoff_for_bonds(double cutoff) const {
    ensure_topology_initialized();
    if (topology) {
      topology->set_cutoff(cutoff);
    }
  }

  /// Set the atom typing method
  void set_atom_typing_method(AtomTypingMethod method) const {
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
  const auto &get_positions(int confId = -1) const { return get_conformer(confId).getPositions(); }
  const auto *get_atom(int idx) const { return mol->getAtomWithIdx(idx); }

  const auto *get_info(int idx) const {
    return static_cast<const RDKit::AtomPDBResidueInfo *>(get_atom(idx)->getMonomerInfo());
  }

  friend class Contacts;

private:
  struct TopologyBuildState {
    std::mutex mutex;
    std::condition_variable cv;
    bool building = false;
  };

  explicit Luni(std::shared_ptr<RDKit::RWMol> valid_mol)
    : mol(valid_mol), topology(std::make_shared<Topology>(valid_mol)), topology_built_(false),
      topology_state_(std::make_shared<TopologyBuildState>()) {}

  void ensure_topology_initialized() const {
    std::call_once(topology_init_once_, [this]() {
      if (!topology) {
        Logger::get_logger()->debug("Initializing topology for configuration");
        topology = std::make_shared<Topology>(mol);
      }
    });
  }

  auto match_smarts_string(std::string sm, std::string atype = "", bool log_values = false) const;

  template <typename T>
  std::vector<T> atom_attrs(std::function<T(const RDKit::Atom *)> func) const;

  template <typename T>
  std::vector<std::reference_wrapper<const T>>
  atom_attrs_ref(std::function<const T &(const RDKit::Atom *)> func) const;

  // FIX: It should not be necessary to have a default constructed RWMol here
  std::shared_ptr<RDKit::RWMol> mol = std::make_shared<RDKit::RWMol>();
  mutable std::shared_ptr<Topology> topology;
  mutable std::atomic<bool> topology_built_{false};
  bool model_origin_ = false; // flag controlling model code path

  mutable std::shared_ptr<TopologyBuildState> topology_state_ = std::make_shared<TopologyBuildState>();
  mutable std::once_flag topology_init_once_;

  std::string file_name_;
  std::vector<int> filtered_indices;
  bool is_in_filtered_state = false;
};

} // namespace lahuta

#endif // LAHUTA_HPP
