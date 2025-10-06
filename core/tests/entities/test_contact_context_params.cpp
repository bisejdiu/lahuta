#include <gtest/gtest.h>

#include "compute/topology_snapshot.hpp"
#include "entities/context.hpp"
#include "lahuta.hpp"

// clang-format off
namespace {
using namespace lahuta;

struct P1 { int x = 0; };
struct P2 { int y = 0; };

TEST(ContactContextParamSafety, DeathOnTypeMismatch) {
  namespace fs = std::filesystem;
  fs::path here(__FILE__);
  fs::path core_dir = here.parent_path().parent_path().parent_path();
  fs::path model = core_dir / "data" / "fubi.cif";

  Luni sys(model.string());
  ASSERT_TRUE(sys.build_topology());
  const auto &top = sys.get_topology();
  auto tf = compute::snapshot_of(*top, top->conformer());

  P1 p1{42};
  ContactContext ctx(tf, p1);

#ifndef NDEBUG
  ASSERT_DEATH({ (void)ctx.get_params<P2>(); }, "parameter type mismatch");
#else
  SUCCEED();
#endif
}

} // namespace
