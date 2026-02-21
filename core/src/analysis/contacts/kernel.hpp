/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<const char*, 3> parts{"besian", "sejdiu", "@gmail.com"}; std::string s;
 *   for (auto p : parts) s.append(p, std::strlen(p));
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_ANALYSIS_CONTACTS_KERNEL_HPP
#define LAHUTA_ANALYSIS_CONTACTS_KERNEL_HPP

#include <memory>
#include <string>
#include <vector>

#include "analysis/contacts/records.hpp"
#include "compute/result.hpp"
#include "compute/topology_snapshot.hpp"
#include "contacts/arpeggio/provider.hpp"
#include "contacts/engine.hpp"
#include "contacts/getcontacts/provider.hpp"
#include "contacts/molstar/provider.hpp"
#include "logging/logging.hpp"
#include "pipeline/task/compute/context.hpp"
#include "pipeline/task/compute/parameters.hpp"
#include "pipeline/task/context.hpp"
#include "pipeline/task/emission.hpp"
#include "serialization/formats.hpp"
#include "serialization/specializations/contacts.hpp"
#include "topology.hpp"

namespace lahuta::analysis {
namespace C = lahuta::compute;
namespace P = lahuta::pipeline;

namespace detail {
template <class Provider>
auto compute_with_provider(const P::ContactsParams &p, const C::TopologySnapshot &ts) -> lahuta::ContactSet {
  lahuta::InteractionEngine<Provider> engine;
  if (p.type.is_all()) return engine.compute(ts);
  return engine.compute(ts, std::optional<lahuta::InteractionTypeSet>{p.type});
}
} // namespace detail

// Computes contacts using MolStar, Arpeggio, or GetContacts providers and
// serializes results to JSON or binary format.
struct ContactsKernel {
  static C::ComputationResult execute(C::DataContext<P::PipelineContext, C::Mut::ReadWrite> &context,
                                      const P::ContactsParams &p) {
    auto &data = context.data();

    try {
      ContactsRecord res;
      res.file_path     = data.item_path;
      res.provider      = p.provider;
      res.contact_types = p.type;
      res.success       = false;
      res.num_contacts  = 0;
      res.frame_index   = static_cast<std::size_t>(data.conformer_id);
      res.topology      = nullptr;

      std::shared_ptr<const Topology> top;
      if (data.ctx) top = data.ctx->topology();
      if (!top) return C::ComputationResult(C::ComputationError("Contacts requires topology in context"));

      // Correctness guard: ensure atom typing matches provider before computing
      try {
        using namespace lahuta::topology;
        auto &label  = AtomTypingComputation<>::label;
        auto &eng    = const_cast<Topology &>(*top).get_engine();
        auto *params = eng.get_parameters<AtomTypingParams>(label);

        auto current_mode  = params ? params->mode : AtomTypingMethod::Molstar;
        auto required_mode = typing_for_provider(p.provider);
        if (current_mode != required_mode) {
          Logger::get_logger()->debug("[contacts:typing] Switching atom typing to {}",
                                      contact_provider_name(p.provider));
          auto &topo_mut = const_cast<Topology &>(*top);

          topo_mut.assign_typing(required_mode);
        }
      } catch (const std::exception &e) {
        Logger::get_logger()->warn("[contacts:typing] Typing guard failed: {}", e.what());
      }

      auto ts = data.ctx ? C::require_topology_snapshot(*data.ctx) : C::snapshot_of(*top);

      Logger::get_logger()->debug(
          "[contacts:compute] Starting contact computation using {} provider for {} contacts",
          contact_provider_name(p.provider),
          interaction_type_set_to_string(p.type));

      switch (p.provider) {
        case ContactProvider::Arpeggio:
          res.contacts = detail::compute_with_provider<ArpeggioContactProvider>(p, ts);
          break;
        case ContactProvider::GetContacts:
          res.contacts = detail::compute_with_provider<GetContactsProvider>(p, ts);
          break;
        case ContactProvider::MolStar:
          res.contacts = detail::compute_with_provider<MolStarContactProvider>(p, ts);
          break;
        default:
          return C::ComputationResult(C::ComputationError("Contacts unsupported provider"));
      }

      Logger::get_logger()->debug(
          "[contacts:compute] Completed contact computation using {} provider, found {} contacts",
          contact_provider_name(p.provider),
          res.contacts.size());

      res.num_contacts = res.contacts.size();
      res.success      = true;
      res.topology     = top;
      if (data.frame) res.frame_index = data.frame->index();
      if (auto fm = data.ctx->frame_metadata(); fm && fm->source_file) {
        res.trajectory_file = fm->source_file;
      }

      std::string payload;
      switch (p.format) {
        case P::ContactsOutputFormat::Json:
          payload = serialization::Serializer<fmt::json, ContactsRecord>::serialize(res);
          break;
        case P::ContactsOutputFormat::Binary:
          payload = serialization::Serializer<fmt::binary, ContactsRecord>::serialize(res);
          break;
        case P::ContactsOutputFormat::JsonCompact:
          payload = serialization::Serializer<fmt::json_compact, ContactsRecord>::serialize(res);
          break;
        default:
          throw std::runtime_error("Contacts unsupported output format");
      }

      P::EmissionList out;
      out.push_back({p.channel, std::move(payload)});

      if (data.ctx) {
        data.ctx->set_text("contacts_file", res.file_path);
        data.ctx->set_text("contacts_provider", std::string(contact_provider_name(p.provider)));
        data.ctx->set_text("contacts_success", res.success ? "1" : "0");
        data.ctx->set_text("contacts_count", std::to_string(res.num_contacts));
      }

      return C::ComputationResult(std::move(out));
    } catch (const std::exception &e) {
      return C::ComputationResult(C::ComputationError(std::string("Contacts failed: ") + e.what()));
    } catch (...) {
      return C::ComputationResult(C::ComputationError("Contacts failed"));
    }
  }
};

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_CONTACTS_KERNEL_HPP
