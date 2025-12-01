#ifndef LAHUTA_PIPELINE_COMPUTE_PARAMETERS_HPP
#define LAHUTA_PIPELINE_COMPUTE_PARAMETERS_HPP

#include <cstdint>
#include <optional>
#include <string>

#include "analysis/contacts/provider.hpp"
#include "compute/parameters.hpp"
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
  constexpr ParameterInterface::TypeId ENSURE_TYPING    = 34;
}

enum class ContactsOutputFormat : uint8_t {
  Json,
  Text,
  Binary,
};

struct SystemReadParams : public ParameterBase<SystemReadParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = param_ids::SYSTEM_READ;
  bool is_model = false;
};

struct BuildTopologyParams : public ParameterBase<BuildTopologyParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = param_ids::BUILD_TOPOLOGY;
  TopologyComputation flags = TopologyComputation::All;
  AtomTypingMethod atom_typing_method = AtomTypingMethod::Molstar;
};

struct ContactsParams : public ParameterBase<ContactsParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = param_ids::CONTACTS;
  analysis::contacts::ContactProvider provider = analysis::contacts::ContactProvider::MolStar;
  InteractionTypeSet type = InteractionTypeSet::all();
  std::string channel = "contacts";
  ContactsOutputFormat format = ContactsOutputFormat::Json;
};

// Ensure that the topology's atom typing matches desired mode. desired = std::nullopt means no preference
struct EnsureTypingParams : public ParameterBase<EnsureTypingParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = param_ids::ENSURE_TYPING;
  std::optional<AtomTypingMethod> desired = std::nullopt;
};

struct DynamicTaskParams : public ParameterBase<DynamicTaskParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = param_ids::DYNAMIC_TASK;
};

} // namespace lahuta::pipeline::compute

#endif // LAHUTA_PIPELINE_COMPUTE_PARAMETERS_HPP
