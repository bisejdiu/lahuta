#pragma once

#include "analysis/system/records.hpp"
#include "compute/result.hpp"
#include "lahuta.hpp"
#include "pipeline/compute/context.hpp"
#include "pipeline/compute/parameters.hpp"
#include <logging.hpp>
#include <models/topology.hpp>

#include <memory>
#include <string>

// clang-format off
namespace lahuta::analysis::system {
using namespace lahuta::pipeline::compute;

struct SystemReadKernel {
  static ComputationResult execute(DataContext<PipelineContext, Mut::ReadWrite>& context, const SystemReadParams& p) {
    try {
      auto& data = context.data();
      std::shared_ptr<Luni> sys;

      if (p.is_model) {
        // Check if model data read from LMDB is available in the context
        std::shared_ptr<const analysis::system::ModelRecord> md = data.ctx ? data.ctx->get_object<analysis::system::ModelRecord>("model_data") : nullptr;
        if (md) {
          auto mol = std::make_shared<RDKit::RWMol>();
          lahuta::build_model_topology(mol, md->data, ModelTopologyMethod::CSR);
          auto s = Luni::create(mol);
          sys = std::make_shared<Luni>(std::move(s));
        } else {
          auto s = Luni::from_model_file(data.item_path);
          sys = std::make_shared<Luni>(std::move(s));
        }
      } else {
        sys = std::make_shared<Luni>(data.item_path);
      }

      if (data.ctx) {
        data.ctx->set_object<Luni>("system", sys);
        // Propagate system construction method
        data.ctx->set_text("system_mode", p.is_model ? std::string("model") : std::string("structure"));
      }
      return ComputationResult(true);
    } catch (const std::exception& e) {
      return ComputationResult(ComputationError(std::string("SystemRead failed: ") + e.what()));
    } catch (...) {
      return ComputationResult(ComputationError("SystemRead failed"));
    }
  }
};

} // namespace lahuta::analysis::system
