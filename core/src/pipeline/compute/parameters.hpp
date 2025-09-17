#pragma once

#include <memory>
#include <string>

#include "analysis/contacts/provider.hpp"
#include "compute/parameters.hpp"
#include "db/db.hpp"
#include "entities/interaction_types.hpp"
#include "topology_flags.hpp"

// clang-format off
namespace lahuta::pipeline::compute {
using namespace lahuta::topology::compute;

namespace param_ids {
  // Using ids >= 30 to avoid clashing with topology (1..6)
  constexpr ParameterInterface::TypeId SYSTEM_READ      = 30;
  constexpr ParameterInterface::TypeId BUILD_TOPOLOGY   = 31;
  constexpr ParameterInterface::TypeId CONTACTS         = 32;
  constexpr ParameterInterface::TypeId DYNAMIC_TASK     = 33;
  constexpr ParameterInterface::TypeId MODEL_FETCH      = 34;
  constexpr ParameterInterface::TypeId ENSURE_TYPING    = 35;
}

struct SystemReadParams : public ParameterBase<SystemReadParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = param_ids::SYSTEM_READ;
  bool is_model = false;
};

struct BuildTopologyParams : public ParameterBase<BuildTopologyParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = param_ids::BUILD_TOPOLOGY;
  TopologyComputation flags = TopologyComputation::All;
  ContactComputerType atom_typing_method = ContactComputerType::Molstar;
};

struct ContactsParams : public ParameterBase<ContactsParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = param_ids::CONTACTS;
  analysis::contacts::ContactProvider provider = analysis::contacts::ContactProvider::MolStar;
  InteractionType type = InteractionType::All;
  std::string channel = "contacts";
  bool json = true; // true -> JSON, false -> TEXT
};

// Ensure that the topology's atom typing matches desired mode. desired = None means no preference
struct EnsureTypingParams : public ParameterBase<EnsureTypingParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = param_ids::ENSURE_TYPING;
  ContactComputerType desired = ContactComputerType::None;
};

struct DynamicTaskParams : public ParameterBase<DynamicTaskParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = param_ids::DYNAMIC_TASK;
};

struct ModelFetchParams : public ParameterBase<ModelFetchParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = param_ids::MODEL_FETCH;
  std::shared_ptr<LMDBDatabase> db; // shared LMDB handle. Reader will be thread_local inside kernel
};

} // namespace lahuta::pipeline::compute
