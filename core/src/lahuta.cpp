/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto f = [](auto arg) {
 *     using T = std::decay_t<decltype(arg)>;
 *     if constexpr (std::disjunction_v<std::is_same<T, const char*>, std::is_same<T, std::string_view>>) return std::string(arg);
 *   };
 *   return f("besian") + f("sejdiu") + f("@gmail.com");
 * }();
 *
 */

#include <stdexcept>
#include <string>

#include <gemmi/gz.hpp>
#include <rdkit/GraphMol/Conformer.h>

#include "analysis/system/model_loader.hpp"
#include "io/convert.hpp"
#include "io/read_file.hpp"
#include "lahuta.hpp"
#include "models/factory.hpp"
#include "models/parser.hpp"
#include "models/topology.hpp"
#include "selections/mol_filters.hpp"

namespace lahuta {

std::shared_ptr<Luni::TopologyBuildState> Luni::ensure_topology_state() const {
  if (topology_state_) return topology_state_;

  std::lock_guard<std::mutex> lk(topology_state_init_mutex_);
  if (!topology_state_) topology_state_ = std::make_shared<TopologyBuildState>();
  return topology_state_;
}

Luni::TopologyWriteGuard::TopologyWriteGuard(std::shared_ptr<TopologyBuildState> state)
    : state_(std::move(state)), lock_() {
  if (!state_) return;

  lock_ = std::unique_lock<std::mutex>(state_->modify_mutex, std::defer_lock);
  std::unique_lock<std::mutex> state_lock(state_->mutex);
  state_->cv.wait(state_lock, [this] { return !state_->building; });
  lock_.lock();
}

Luni::TopologyWriteGuard Luni::acquire_topology_write_guard() const {
  return TopologyWriteGuard(ensure_topology_state());
}

// clang-format off
Luni::Luni(Luni&& other) noexcept
  : mol             (std::move(other.mol)),
    topology        (std::move(other.topology)),
    model_origin_   (other.model_origin_),
    topology_state_ (std::move(other.topology_state_)),
    file_name_      (std::move(other.file_name_)),
    filtered_indices(std::move(other.filtered_indices)),
    is_in_filtered_state(other.is_in_filtered_state) {

    topology_built_.store(other.topology_built_.load(std::memory_order_acquire), std::memory_order_release);
}

Luni& Luni::operator=(Luni&& other) noexcept {
  if (this == &other) return *this;
  mol              = std::move(other.mol);
  topology         = std::move(other.topology);
  model_origin_    = other.model_origin_;
  topology_state_  = std::move(other.topology_state_);
  file_name_       = std::move(other.file_name_);
  filtered_indices = std::move(other.filtered_indices);
  is_in_filtered_state = other.is_in_filtered_state;

  topology_built_.store(other.topology_built_.load(std::memory_order_acquire), std::memory_order_release);
  return *this;
}

// clang-format on
Luni::Luni(std::string file_name) : file_name_(file_name) {
  Logger::get_logger()->debug("Processing file: {}", file_name_);
  mol = read_and_make_molecule(gemmi::MaybeGzipped(file_name_));

  /*auto st = gemmi::read_structure_gz(file_name_);*/
  /*mol = std::make_shared<RDKit::RWMol>();*/
  /*RDKit::Conformer *conformer = new RDKit::Conformer();*/
  /*create_RDKit_repr(*mol, st, *conformer);*/
  /*mol->updatePropertyCache(false);*/
  /*mol->addConformer(conformer, true);*/
}

Luni Luni::create(const gemmi::Structure &st) {
  auto mol = create_RDKit(st);
  return Luni(mol);
}

Luni Luni::create(const IR &ir) {
  auto mol = std::make_shared<RDKit::RWMol>();
  IR_to_RWMol(*mol, ir);
  return Luni(mol);
}

Luni::Luni(std::string file_name, ModelFileTag) : file_name_(file_name) {
  // Makes sure model pools are initialized for single threaded use by default.
  // The pipeline path will reinitialize with the correct thread count.
  try {
    // Pool initialization is idempotent
    InfoPoolFactory::initialize(1);
    BondPoolFactory::initialize(1);
    AtomPoolFactory::initialize(1);
  } catch (...) {
  }

  model_origin_ = true;
  try {
    auto parsed = analysis::load_model_parser_result(file_name_);
    mol         = analysis::build_model_molecule(parsed, ModelTopologyMethod::CSR);
  } catch (const std::exception &e) {
    Logger::get_logger()->critical("Exception processing file {}: {}", file_name_, e.what());
  } catch (...) {
    Logger::get_logger()->critical("Unknown exception processing file: {}", file_name_);
  }
}

Luni Luni::from_model_data(const ModelParserResult &data) {
  return Luni::from_model_data(data, ModelTopologyMethod::CSR);
}

Luni Luni::from_model_data(const ModelParserResult &data, ModelTopologyMethod method) {
  auto mol = std::make_shared<RDKit::RWMol>();
  if (!build_model_topology(mol, data, method)) {
    throw std::runtime_error("Failed to build model topology");
  }
  return Luni::create(mol, TopologyBuildMode::Model);
}

bool Luni::build_topology(std::optional<TopologyBuildingOptions> tops) const {
  auto state = ensure_topology_state();
  std::unique_lock<std::mutex> state_lock(state->mutex);
  while (state->building) {
    state->cv.wait(state_lock, [state]() { return !state->building; });
  }

  if (topology_built_.load(std::memory_order_acquire)) {
    bool complete = true;
    if (!model_origin_) {
      ensure_topology_initialized();
      complete = topology->has_computed(TopologyComputation::Complete);
    }
    if (complete) return true;
  }

  state->building = true;
  std::unique_lock<std::mutex> modify_lock(state->modify_mutex);
  state_lock.unlock();

  bool built = false;
  try {
    ensure_topology_initialized();

    if (!topology_built_.load(std::memory_order_acquire)) {
      if (tops) {
        built = topology->build(*tops);
      } else {
        TopologyBuildingOptions opts;
        if (model_origin_) opts.mode = TopologyBuildMode::Model;
        built = topology->build(opts);
      }
      if (built) {
        topology_built_.store(true, std::memory_order_release);
      }
    } else {
      built = true;
    }

    if (built && !model_origin_) {
      if (!topology->has_computed(TopologyComputation::Complete)) {
        built = topology->ensure_computed(TopologyComputation::Complete);
      }
    }
  } catch (const std::exception &e) {
    Logger::get_logger()->error("Error building topology: {}", e.what());
    built = false;
  }

  state_lock.lock();
  state->building = false;
  state_lock.unlock();
  modify_lock.unlock();
  state->cv.notify_all();

  return built;
}

bool Luni::build_topology(const TopologyBuildingOptions &tops, TopologyComputation include) const {
  auto state = ensure_topology_state();
  std::unique_lock<std::mutex> state_lock(state->mutex);

  while (state->building) {
    state->cv.wait(state_lock, [state]() { return !state->building; });
  }

  if (topology_built_.load(std::memory_order_acquire)) {
    bool complete = (include == TopologyComputation::None);
    if (!complete) {
      ensure_topology_initialized();
      complete = topology->has_computed(include);
    }
    if (complete) return true;
  }

  state->building = true;
  std::unique_lock<std::mutex> modify_lock(state->modify_mutex);
  state_lock.unlock();

  bool ok = false;
  try {
    ensure_topology_initialized();

    if (tops.mode == TopologyBuildMode::Model) {
      if (!topology_built_.load(std::memory_order_acquire)) {
        ok = topology->build(tops);
        if (ok) topology_built_.store(true, std::memory_order_release);
      } else {
        ok = true;
      }
    } else {
      bool initialized_ok = true;
      if (!topology_built_.load(std::memory_order_acquire)) {
        initialized_ok = topology->initialize(tops);
      }

      if (initialized_ok) {
        if (include == TopologyComputation::None) {
          ok = true;
        } else {
          ok = topology->ensure_computed(include);
          if (ok) topology_built_.store(true, std::memory_order_release);
        }
      } else {
        ok = false;
      }
    }
  } catch (const std::exception &e) {
    Logger::get_logger()->error("Error building topology: {}", e.what());
    ok = false;
  }

  state_lock.lock();
  state->building = false;
  state_lock.unlock();
  modify_lock.unlock();
  state->cv.notify_all();

  return ok;
}

Luni Luni::reset_topology() const {
  if (file_name_.empty()) {
    throw std::runtime_error("reset_topology requires a file-backed LahutaSystem");
  }

  if (model_origin_) {
    return Luni::from_model_file(file_name_);
  }

  return Luni(file_name_);
}

const std::vector<std::string> Luni::symbols() const {
  return atom_attrs<std::string>([](const RDKit::Atom *atom) { return atom->getSymbol(); });
}

const std::vector<std::string> Luni::names() const {
  return atom_attrs<std::string>([](const RDKit::Atom *atom) -> const std::string & {
    auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
    return info->getName();
  });
}

const std::vector<int> Luni::indices() const {
  return atom_attrs<int>([](const RDKit::Atom *atom) { return atom->getIdx(); });
}

const std::vector<int> Luni::atomic_numbers() const {
  return atom_attrs<int>([](const RDKit::Atom *atom) { return atom->getAtomicNum(); });
}

const std::vector<std::string> Luni::elements() const {
  const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
  return atom_attrs<std::string>(
      [&tbl](const RDKit::Atom *atom) { return tbl->getElementSymbol(atom->getAtomicNum()); });
}

const std::vector<std::string> Luni::resnames() const {
  return atom_attrs<std::string>([](const RDKit::Atom *atom) -> const std::string & {
    auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
    return info->getResidueName();
  });
}

const std::vector<int> Luni::resids() const {
  return atom_attrs<int>([](const RDKit::Atom *atom) -> int {
    auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
    return info->getResidueNumber();
  });
}

const std::vector<std::string> Luni::chainlabels() const {
  return atom_attrs<std::string>([](const RDKit::Atom *atom) -> const std::string & {
    auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
    return info->getChainId();
  });
}

template <typename T>
std::vector<T> Luni::atom_attrs(std::function<T(const RDKit::Atom *)> func) const {
  std::vector<T> attrs;
  attrs.reserve(mol->getNumAtoms());
  for (const auto atom : mol->atoms()) {
    attrs.push_back(func(atom));
  }
  return attrs;
}

template <typename T>
std::vector<std::reference_wrapper<const T>>
Luni::atom_attrs_ref(std::function<const T &(const RDKit::Atom *)> func) const {
  std::vector<std::reference_wrapper<const T>> attributes;
  attributes.reserve(mol->getNumAtoms());
  for (const auto atom : mol->atoms()) {
    attributes.push_back(std::cref(func(atom)));
  }
  return attributes;
}

auto Luni::match_smarts_string(std::string sm, std::string atype, bool log_values) const {

  if (!mol->getRingInfo()->isInitialized()) {
    RDKit::MolOps::symmetrizeSSSR(*mol);
  }
  std::vector<RDKit::MatchVectType> match_list;
  auto sm_mol = RDKit::SmartsToMol(sm);
  mol->updatePropertyCache(false);
  RDKit::SubstructMatch(*mol, *sm_mol, match_list);
  return match_list;
};

Luni Luni::filter(std::vector<int> &atom_indices) const {

  auto new_mol = filter_with_bonds(*mol, atom_indices);
  new_mol.updatePropertyCache(false);

  auto new_luni                 = Luni::create(std::make_shared<RDKit::RWMol>(new_mol));
  new_luni.is_in_filtered_state = true;
  return new_luni;
}

} // namespace lahuta
