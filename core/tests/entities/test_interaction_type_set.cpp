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

#include <gtest/gtest.h>

#include "entities/interaction_types.hpp"

// clang-format off
namespace {
using namespace lahuta;

TEST(InteractionTypeSetTest, DefaultStateIsEmptyAndNotAll) {
  InteractionTypeSet set;
  EXPECT_TRUE(set.empty());
  EXPECT_FALSE(set.is_all());
  EXPECT_EQ(set.count(), 0u);
  EXPECT_FALSE(set.contains(InteractionType::HydrogenBond));
  EXPECT_FALSE(set.contains(InteractionType::All));
}

TEST(InteractionTypeSetTest, AddingTypesTracksMembershipAndCount) {
  InteractionTypeSet set;
  set.add(InteractionType::HydrogenBond);
  EXPECT_FALSE(set.empty());
  EXPECT_TRUE(set.contains(InteractionType::HydrogenBond));
  EXPECT_FALSE(set.contains(InteractionType::WeakHydrogenBond));
  EXPECT_EQ(set.count(), 1u);

  auto members = set.members();
  ASSERT_EQ(members.size(), 1u);
  EXPECT_EQ(members.front(), InteractionType::HydrogenBond);
}

TEST(InteractionTypeSetTest, BitwiseOperatorsCombineUniqueMembers) {
  auto combined = InteractionType::HydrogenBond | InteractionType::WeakHydrogenBond;
  EXPECT_TRUE(combined.contains(InteractionType::HydrogenBond));
  EXPECT_TRUE(combined.contains(InteractionType::WeakHydrogenBond));
  EXPECT_FALSE(combined.contains(InteractionType::PiStacking));
  EXPECT_EQ(combined.count(), 2u);

  combined |= InteractionType::HydrogenBond; // duplicate should be ignored
  EXPECT_EQ(combined.count(), 2u);

  InteractionTypeSet via_set(InteractionType::HydrogenBond);
  via_set |= InteractionType::WeakHydrogenBond;
  EXPECT_EQ(combined, via_set);
}

TEST(InteractionTypeSetTest, AllSentinelDominatesAndExpandsMembers) {
  InteractionTypeSet set;
  set.add(InteractionType::HydrogenBond);
  set.add(InteractionType::All);

  EXPECT_TRUE(set.is_all());
  EXPECT_FALSE(set.empty());
  EXPECT_TRUE(set.contains(InteractionType::HydrogenBond));
  EXPECT_TRUE(set.contains(InteractionType::All));
  EXPECT_EQ(set.count(), all_interaction_types().size());

  auto names = set.members(false);
  ASSERT_EQ(names.size(), 1u);
  EXPECT_EQ(names.front(), InteractionType::All);

  auto expanded = set.members();
  EXPECT_EQ(expanded.size(), all_interaction_types().size());
}

TEST(InteractionTypeSetTest, UnionWithAllProducesAll) {
  InteractionTypeSet lhs(InteractionType::HydrogenBond);
  InteractionTypeSet rhs({InteractionType::WeakHydrogenBond, InteractionType::MetalCoordination});

  auto merged = lhs | rhs;
  EXPECT_FALSE(merged.is_all());
  EXPECT_EQ(merged.count(), 3u);
  EXPECT_TRUE(merged.contains(InteractionType::HydrogenBond));
  EXPECT_TRUE(merged.contains(InteractionType::WeakHydrogenBond));
  EXPECT_TRUE(merged.contains(InteractionType::MetalCoordination));

  lhs |= InteractionTypeSet::all();
  EXPECT_TRUE(lhs.is_all());

  auto with_all = rhs | InteractionTypeSet::all();
  EXPECT_TRUE(with_all.is_all());
  EXPECT_TRUE(with_all.contains(InteractionType::PiStacking));
}

} // namespace
