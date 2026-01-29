#ifndef LAHUTA_PIPELINE_TASK_COMPUTE_PARAMETERS_HPP
#define LAHUTA_PIPELINE_TASK_COMPUTE_PARAMETERS_HPP

#include <cstdint>
#include <memory>
#include <optional>
#include <string>

#include "analysis/contacts/provider.hpp"
#include "analysis/sasa/records.hpp"
#include "analysis/sasa/sasa.hpp"
#include "compute/parameters.hpp"
#include "entities/interaction_types.hpp"
#include "topology_flags.hpp"

namespace lahuta::pipeline {
namespace C = lahuta::compute;

namespace param_ids {
// Using ids >= 30 to avoid clashing with topology (1..6)
constexpr C::ParameterInterface::TypeId SYSTEM_READ    = 30;
constexpr C::ParameterInterface::TypeId BUILD_TOPOLOGY = 31;
constexpr C::ParameterInterface::TypeId CONTACTS       = 32;
constexpr C::ParameterInterface::TypeId DYNAMIC_TASK   = 33;
constexpr C::ParameterInterface::TypeId ENSURE_TYPING  = 34;
constexpr C::ParameterInterface::TypeId SASA_SR        = 35;
} // namespace param_ids

enum class ContactsOutputFormat : uint8_t {
  Json,
  Text,
  Binary,
};

struct SystemReadParams : public C::ParameterBase<SystemReadParams> {
  static constexpr C::ParameterInterface::TypeId TYPE_ID = param_ids::SYSTEM_READ;

  bool is_model = false;
};

struct BuildTopologyParams : public C::ParameterBase<BuildTopologyParams> {
  static constexpr C::ParameterInterface::TypeId TYPE_ID = param_ids::BUILD_TOPOLOGY;

  TopologyComputation flags           = TopologyComputation::All;
  AtomTypingMethod atom_typing_method = AtomTypingMethod::Molstar;
};

struct ContactsParams : public C::ParameterBase<ContactsParams> {
  static constexpr C::ParameterInterface::TypeId TYPE_ID = param_ids::CONTACTS;

  analysis::ContactProvider provider = analysis::ContactProvider::MolStar;
  InteractionTypeSet type            = InteractionTypeSet::all();
  std::string channel                = "contacts";
  ContactsOutputFormat format        = ContactsOutputFormat::Json;
};

struct SasaSrParams : public C::ParameterBase<SasaSrParams> {
  static constexpr C::ParameterInterface::TypeId TYPE_ID = param_ids::SASA_SR;

  analysis::SasaParams params;
  std::string channel = std::string(analysis::SasaSrOutputChannel);
  bool include_total = false;
  std::shared_ptr<analysis::SasaSrCounters> counters;
};

// Ensure that the topology's atom typing matches desired mode. desired = std::nullopt means no preference
struct EnsureTypingParams : public C::ParameterBase<EnsureTypingParams> {
  static constexpr C::ParameterInterface::TypeId TYPE_ID = param_ids::ENSURE_TYPING;

  std::optional<AtomTypingMethod> desired = std::nullopt;
};

struct DynamicTaskParams : public C::ParameterBase<DynamicTaskParams> {
  static constexpr C::ParameterInterface::TypeId TYPE_ID = param_ids::DYNAMIC_TASK;
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_TASK_COMPUTE_PARAMETERS_HPP
