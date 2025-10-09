#pragma once

#include <memory>
#include <string>
#include <vector>

#include "analysis/contacts/records.hpp"
#include "compute/result.hpp"
#include "compute/topology_snapshot.hpp"
#include "contacts/arpeggio/provider.hpp"
#include "contacts/engine.hpp"
#include "contacts/molstar/provider.hpp"
#include "logging.hpp"
#include "pipeline/compute/context.hpp"
#include "pipeline/compute/parameters.hpp"
#include "pipeline/dynamic/types.hpp"
#include "serialization/formats.hpp"
#include "serialization/serializer.hpp"
#include "topology.hpp"

// clang-format off
namespace lahuta::analysis::contacts {
using namespace lahuta::pipeline::compute;

// Computes contacts using MolStar or Arpeggio providers and serializes results to JSON or text format.
struct ContactsKernel {
  static ComputationResult execute(DataContext<PipelineContext, Mut::ReadWrite>& context, const ContactsParams& p) {
    auto& data = context.data();

    try {
      ContactsRecord res;
      res.file_path    = data.item_path;
      res.provider     = p.provider;
      res.contact_type = p.type;
      res.success      = false;
      res.num_contacts = 0;
      res.topology     = nullptr;

      std::shared_ptr<const Topology> top;
      if (data.ctx) top = data.ctx->get_object<const Topology>(pipeline::CTX_TOPOLOGY_KEY);
      if (!top) return ComputationResult(ComputationError("Contacts requires topology in context"));

      // Correctness guard: ensure atom typing matches provider before computing
      try {
        auto& label  = topology::AtomTypingComputation<>::label;
        auto& eng    = const_cast<Topology&>(*top).get_engine();
        auto* params = eng.get_parameters<topology::AtomTypingParams>(label);

        auto current_mode = params ? params->mode : AtomTypingMethod::Molstar;
        auto required_mode = typing_for_provider(p.provider);
        if (current_mode != required_mode) {
          Logger::get_logger()->info("ContactsKernel: switching atom typing to {} for contacts computation", contact_provider_name(p.provider));
          auto& topo_mut = const_cast<Topology&>(*top);

          topo_mut.assign_typing(required_mode);
        }
      } catch (const std::exception& e) {
        Logger::get_logger()->error("ContactsKernel: typing guard failed: {}", e.what());
      }

      auto ts = data.ctx ? compute::require_topology_snapshot(*data.ctx)
                         : compute::snapshot_of(*top);

      if (p.provider == ContactProvider::Arpeggio) {
        InteractionEngine<ArpeggioContactProvider> engine;
        res.contacts = (p.type == InteractionType::All)
          ? engine.compute(ts)
          : engine.compute(ts, p.type);
      } else {
        InteractionEngine<MolStarContactProvider> engine;
        res.contacts = (p.type == InteractionType::All)
          ? engine.compute(ts)
          : engine.compute(ts, p.type);
      }
      res.num_contacts = res.contacts.size();
      res.success  = true;
      res.topology = top;

      std::string payload;
      switch (p.format) {
        case ContactsOutputFormat::Binary:
          payload = serialization::Serializer<fmt::binary, ContactsRecord>::serialize(res);
          break;
        case ContactsOutputFormat::Json:
          payload = serialization::Serializer<fmt::json, ContactsRecord>::serialize(res);
          break;
        case ContactsOutputFormat::Text:
          payload = serialization::Serializer<fmt::text, ContactsRecord>::serialize(res);
          break;
        default:
          throw std::runtime_error("ContactsKernel: unsupported output format");
      }

      pipeline::dynamic::EmissionList out;
      out.push_back({p.channel, std::move(payload)});

      if (data.ctx) {
        data.ctx->set_text("contacts_file",     res.file_path);
        data.ctx->set_text("contacts_provider", (p.provider == ContactProvider::Arpeggio) ? "arpeggio" : "molstar");
        data.ctx->set_text("contacts_success",  res.success ? "1" : "0");
        data.ctx->set_text("contacts_count",    std::to_string(res.num_contacts));
      }

      return ComputationResult(std::move(out));
    } catch (const std::exception& e) {
      return ComputationResult(ComputationError(std::string("Contacts failed: ") + e.what()));
    } catch (...) {
      return ComputationResult(ComputationError("Contacts failed"));
    }
  }
};

} // namespace lahuta::analysis::contacts
