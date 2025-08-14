#include <gtest/gtest.h>

#include "lahuta.hpp"
#include "logging.hpp"
#include "models/factory.hpp"

// Verifies that constructing a model-fast-path system followed by constructing
// a pipeline-style system (or vice-versa) in the same process does not crash
// due to pool re-initialization or lifetime issues.
//

TEST(PoolsReinitLifetime, DirectThenDirectModelFastPath) {
  using namespace lahuta;
  Logger::get_logger()->set_level(spdlog::level::warn);

  // A baseline init
  InfoPoolFactory::initialize(1);
  BondPoolFactory::initialize(1);
  AtomPoolFactory::initialize(1);

  namespace fs = std::filesystem;
  fs::path here(__FILE__);
  fs::path core_dir = here.parent_path().parent_path().parent_path();
  fs::path model_path = core_dir / "data" / "fubi.cif";

  // Direct fast-path system
  {
    lahuta::Luni sys = lahuta::Luni::from_model_file(model_path.string());
    ASSERT_TRUE(sys.has_topology_built());
    EXPECT_GT(sys.n_atoms(), 10);
  }

  // Reinitialize pools to simulate a second-stage setup (e.g., pipeline sizing)
  InfoPoolFactory::initialize(1);
  BondPoolFactory::initialize(1);
  AtomPoolFactory::initialize(1);

  // Another direct fast-path system in the same process
  {
    lahuta::Luni sys = lahuta::Luni::from_model_file(model_path.string());
    ASSERT_TRUE(sys.has_topology_built());
    EXPECT_GT(sys.n_atoms(), 10);
  }
}
