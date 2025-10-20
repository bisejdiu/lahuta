#include <gtest/gtest.h>

#include <filesystem>

#include "lahuta.hpp"
#include "logging.hpp"

#include "contacts/engine.hpp"
#include "contacts/molstar/provider.hpp"
#include "entities/find_contacts.hpp"
#include "entities/records.hpp"
#include "entities/search/config.hpp"

#include "compute/topology_snapshot.hpp"
#include "pipeline/dynamic/keys.hpp"
#include "pipeline/dynamic/types.hpp"

// clang-format off
namespace {

using namespace lahuta;

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
  namespace fs = std::filesystem;
  fs::path data = core_data_dir() / "fubi.cif";

  Luni sys(data.string());
  ASSERT_TRUE(sys.build_topology());
  auto top = sys.get_topology();
  ASSERT_TRUE(top != nullptr);

  // Predicates and tester (simple hydrophobic-like within 4.5A)
  auto pred_atom = +[](const AtomRec &rec) { return (rec.type & AtomType::Hydrophobic) == AtomType::Hydrophobic; };
  auto tester = +[](uint32_t, uint32_t, float d2, const ContactContext&) {
    return d2 <= 4.5f * 4.5f ? InteractionType::Hydrophobic : InteractionType::None;
  };
  search::SearchOptions opts; opts.distance_max = 4.5;

  _EmptyParams params;

  // Path A: explicit snapshot from topology default conformer
  auto tf_exp = compute::snapshot_of(*top, top->conformer());
  ContactContext ctx_exp(tf_exp, params);
  ContactSet a = find_contacts(ctx_exp, pred_atom, opts, tester);

  // Path B: per-frame-style snapshot via TaskContext carrying a conformer
  pipeline::dynamic::TaskContext tctx;
  auto conf_ptr = std::make_shared<RDKit::Conformer>(top->conformer());
  tctx.set_object<RDKit::Conformer>(pipeline::CTX_CONFORMER_KEY, conf_ptr);
  auto tf_frame = compute::snapshot_of(*top, &tctx);
  ContactContext ctx_frame(tf_frame, params);
  ContactSet b = find_contacts(ctx_frame, pred_atom, opts, tester);

  EXPECT_EQ(a.size(), b.size());
  EXPECT_TRUE((a ^ b).empty());
}

TEST(ContactsEngineSnapshotEquivalence, MolStarProvider) {
  Logger::get_logger()->set_level(spdlog::level::info);
  namespace fs = std::filesystem;
  fs::path data = core_data_dir() / "fubi.cif";

  Luni sys(data.string());
  ASSERT_TRUE(sys.build_topology());
  auto top = sys.get_topology();
  ASSERT_TRUE(top != nullptr);

  InteractionEngine<MolStarContactProvider> engine;

  // Path A: explicit snapshot
  auto tf_exp = compute::snapshot_of(*top, top->conformer());
  ContactSet a = engine.compute(tf_exp);

  // Path B: per-frame-style snapshot via TaskContext conformer
  pipeline::dynamic::TaskContext tctx;
  auto conf_ptr = std::make_shared<RDKit::Conformer>(top->conformer());
  tctx.set_object<RDKit::Conformer>(pipeline::CTX_CONFORMER_KEY, conf_ptr);
  auto tf_frame = compute::snapshot_of(*top, &tctx);
  ContactSet b = engine.compute(tf_frame);

  EXPECT_EQ(a.size(), b.size());
  EXPECT_TRUE((a ^ b).empty());
}

} // namespace
