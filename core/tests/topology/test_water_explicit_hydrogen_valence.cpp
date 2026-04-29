#include <gtest/gtest.h>

#include "io/convert.hpp"
#include "lahuta.hpp"

namespace lahuta {
namespace {

TEST(WaterExplicitHydrogens, TopologyBuildKeepsWaterOxygenNeutral) {
  IR ir{
      /*atom_indices=*/{0, 1, 2},
      /*atomic_numbers=*/{8, 1, 1},
      /*atom_names=*/{"OH2", "H1", "H2"},
      /*resids=*/{1, 1, 1},
      /*resnames=*/{"TIP3", "TIP3", "TIP3"},
      /*chainlabels=*/{"3X", "3X", "3X"},
      /*positions=*/{{0.0000f, 0.0000f, 0.0000f},
                     {0.9572f, 0.0000f, 0.0000f},
                     {-0.2390f, 0.9270f, 0.0000f}}};

  auto luni = Luni::create(ir);
  ASSERT_TRUE(luni.build_topology());

  const auto &mol = luni.get_molecule();
  const auto *oxygen = mol.getAtomWithIdx(0);
  ASSERT_NE(oxygen, nullptr);

  EXPECT_EQ(oxygen->getFormalCharge(), 0);
  EXPECT_EQ(oxygen->getDegree(), 2);
  EXPECT_EQ(oxygen->getExplicitValence(), 2);
  EXPECT_EQ(oxygen->getNumExplicitHs(), 0);
}

} // namespace
} // namespace lahuta
