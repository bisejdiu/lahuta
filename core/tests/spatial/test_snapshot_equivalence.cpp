/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::array<std::string_view, 3> parts{"besian", "sejdiu", "@gmail.com"}; std::string dst;
 *   std::for_each(parts.begin(), parts.end(), [&dst](std::string_view p) { dst += p; });
 *   return dst;
 * }();
 *
 */

#include <filesystem>

#include <gtest/gtest.h>

#include "compute/topology_snapshot.hpp"
#include "contacts/engine.hpp"
#include "contacts/molstar/provider.hpp"
#include "entities/find_contacts.hpp"
#include "entities/records.hpp"
#include "entities/search/config.hpp"
#include "lahuta.hpp"
#include "logging/logging.hpp"
#include "pipeline/task/context.hpp"
#include "pipeline/task/keys.hpp"

namespace {

using namespace lahuta;
namespace C = lahuta::compute;
namespace P = lahuta::pipeline;

struct _EmptyParams final {};

static std::filesystem::path core_data_dir() {
  namespace fs = std::filesystem;
  fs::path here(__FILE__);
  fs::path core_dir = here.parent_path().parent_path().parent_path();
  // core/tests/spatial -> core
  return core_dir / "data";
}

TEST(FindersSnapshotEquivalence, AtomHydrophobicSelf) {
  Logger::get_logger()->set_level(spdlog::level::info);
  namespace fs  = std::filesystem;
  fs::path data = core_data_dir() / "fubi.cif";

  Luni sys(data.string());
  ASSERT_TRUE(sys.build_topology());
  auto top = sys.get_topology();
  ASSERT_TRUE(top != nullptr);

  // Predicates and tester (simple hydrophobic-like within 4.5A)
  auto pred_atom = +[](const AtomRec &rec) {
    return (rec.type & AtomType::Hydrophobic) == AtomType::Hydrophobic;
  };
  auto tester = +[](uint32_t, uint32_t, float d2, const ContactContext &) {
    return d2 <= 4.5f * 4.5f ? InteractionType::Hydrophobic : InteractionType::None;
  };
  search::SearchOptions opts;
  opts.distance_max = 4.5;

  _EmptyParams params;

  // Path A: explicit snapshot from topology default conformer
  auto tf_exp = C::snapshot_of(*top, top->conformer());
  ContactContext ctx_exp(tf_exp, params);
  ContactSet a = find_contacts(ctx_exp, pred_atom, opts, tester);

  // Path B: per-frame-style snapshot via TaskContext carrying a conformer
  P::TaskContext tctx;
  auto conf_ptr = std::make_shared<RDKit::Conformer>(top->conformer());
  tctx.set_object<RDKit::Conformer>(P::CTX_CONFORMER_KEY, conf_ptr);
  auto tf_frame = C::snapshot_of(*top, &tctx);
  ContactContext ctx_frame(tf_frame, params);
  ContactSet b = find_contacts(ctx_frame, pred_atom, opts, tester);

  EXPECT_EQ(a.size(), b.size());
  EXPECT_TRUE((a ^ b).empty());
}

TEST(ContactsEngineSnapshotEquivalence, MolStarProvider) {
  Logger::get_logger()->set_level(spdlog::level::info);
  namespace fs  = std::filesystem;
  fs::path data = core_data_dir() / "fubi.cif";

  Luni sys(data.string());
  ASSERT_TRUE(sys.build_topology());
  auto top = sys.get_topology();
  ASSERT_TRUE(top != nullptr);

  InteractionEngine<MolStarContactProvider> engine;

  // Path A: explicit snapshot
  auto tf_exp  = C::snapshot_of(*top, top->conformer());
  ContactSet a = engine.compute(tf_exp);

  // Path B: per-frame-style snapshot via TaskContext conformer
  P::TaskContext tctx;
  auto conf_ptr = std::make_shared<RDKit::Conformer>(top->conformer());
  tctx.set_object<RDKit::Conformer>(P::CTX_CONFORMER_KEY, conf_ptr);
  auto tf_frame = C::snapshot_of(*top, &tctx);
  ContactSet b  = engine.compute(tf_frame);

  EXPECT_EQ(a.size(), b.size());
  EXPECT_TRUE((a ^ b).empty());
}

} // namespace
