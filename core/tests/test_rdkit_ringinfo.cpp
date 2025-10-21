#include <gtest/gtest.h>

#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/RingInfo.h>
#include <rdkit/GraphMol/Rings.h>

namespace {

using namespace RDKit;

TEST(RingInfoAddAllRings, HandlesUnsortedMaxAtomIndices) {
  RWMol mol;
  for (int i = 0; i < 14; ++i) {
    mol.addAtom(new Atom(6));
  }

  mol.addBond(10, 11, Bond::AROMATIC);
  mol.addBond(11, 12, Bond::AROMATIC);
  mol.addBond(12, 13, Bond::AROMATIC);
  mol.addBond(13, 10, Bond::AROMATIC);

  mol.addBond(0, 1, Bond::AROMATIC);
  mol.addBond(1, 2, Bond::AROMATIC);
  mol.addBond(2, 3, Bond::AROMATIC);
  mol.addBond(3, 0, Bond::AROMATIC);

  auto *conf = new Conformer(mol.getNumAtoms());
  for (int i = 0; i < mol.getNumAtoms(); ++i) {
    conf->setAtomPos(i, RDGeom::Point3D(i, 0, 0));
  }
  mol.addConformer(conf, true);

  mol.getRingInfo()->initialize(FIND_RING_TYPE_SYMM_SSSR);

  VECT_INT_VECT rings, bonds;
  rings.push_back({10, 11, 12, 13});
  rings.push_back({0, 1, 2, 3});

  RingUtils::convertToBonds(rings, bonds, mol);
  ASSERT_EQ(bonds.size(), 2u);

  EXPECT_NO_THROW({ mol.getRingInfo()->addAllRings(rings, bonds); });
  EXPECT_EQ(mol.getRingInfo()->numRings(), 2u);
  EXPECT_EQ(mol.getRingInfo()->atomMembers(10).size(), 1u);
  EXPECT_EQ(mol.getRingInfo()->atomMembers(0).size(), 1u);
}

TEST(RingUtilsConvertToBonds, EmptyRingVectors) {
  RWMol mol;
  for (int i = 0; i < 10; ++i) {
    mol.addAtom(new Atom(6));
  }

  mol.addBond(0, 1, Bond::SINGLE);
  mol.addBond(1, 2, Bond::SINGLE);
  mol.addBond(2, 3, Bond::SINGLE);
  mol.addBond(3, 0, Bond::SINGLE);

  VECT_INT_VECT rings, bonds;

  // Test with empty rings
  EXPECT_NO_THROW(RingUtils::convertToBonds(rings, bonds, mol));
  EXPECT_EQ(bonds.size(), 0u);
}

TEST(RingUtilsConvertToBonds, SingleRing) {
  RWMol mol;
  for (int i = 0; i < 4; ++i) {
    mol.addAtom(new Atom(6));
  }

  mol.addBond(0, 1, Bond::SINGLE);
  mol.addBond(1, 2, Bond::SINGLE);
  mol.addBond(2, 3, Bond::SINGLE);
  mol.addBond(3, 0, Bond::SINGLE);

  VECT_INT_VECT rings, bonds;
  rings.push_back({0, 1, 2, 3});

  EXPECT_NO_THROW(RingUtils::convertToBonds(rings, bonds, mol));
  EXPECT_EQ(bonds.size(), 1u);
  EXPECT_EQ(bonds[0].size(), 4u);
}

TEST(RingUtilsConvertToBonds, MultipleRings) {
  RWMol mol;
  for (int i = 0; i < 8; ++i) {
    mol.addAtom(new Atom(6));
  }

  // Two separate 4-membered rings
  mol.addBond(0, 1, Bond::SINGLE);
  mol.addBond(1, 2, Bond::SINGLE);
  mol.addBond(2, 3, Bond::SINGLE);
  mol.addBond(3, 0, Bond::SINGLE);

  mol.addBond(4, 5, Bond::SINGLE);
  mol.addBond(5, 6, Bond::SINGLE);
  mol.addBond(6, 7, Bond::SINGLE);
  mol.addBond(7, 4, Bond::SINGLE);

  VECT_INT_VECT rings, bonds;
  rings.push_back({0, 1, 2, 3});
  rings.push_back({4, 5, 6, 7});

  EXPECT_NO_THROW(RingUtils::convertToBonds(rings, bonds, mol));
  EXPECT_EQ(bonds.size(), 2u);
  EXPECT_EQ(bonds[0].size(), 4u);
  EXPECT_EQ(bonds[1].size(), 4u);
}

TEST(RingUtilsConvertToBonds, NonContiguousAtomIndices) {
  RWMol mol;
  for (int i = 0; i < 20; ++i) {
    mol.addAtom(new Atom(6));
  }

  // Create a ring with non-contiguous atom indices
  mol.addBond(15, 17, Bond::SINGLE);
  mol.addBond(17, 19, Bond::SINGLE);
  mol.addBond(19, 18, Bond::SINGLE);
  mol.addBond(18, 15, Bond::SINGLE);

  VECT_INT_VECT rings, bonds;
  rings.push_back({15, 17, 19, 18});

  EXPECT_NO_THROW(RingUtils::convertToBonds(rings, bonds, mol));
  EXPECT_EQ(bonds.size(), 1u);
  EXPECT_EQ(bonds[0].size(), 4u);
}

TEST(RingUtilsConvertToBonds, LargeNumberOfRings) {
  RWMol mol;
  const int num_atoms = 100;
  const int num_rings = 25;

  for (int i = 0; i < num_atoms; ++i) {
    mol.addAtom(new Atom(6));
  }

  // Create many 4-membered rings
  for (int ring = 0; ring < num_rings; ++ring) {
    int base = ring * 4;
    mol.addBond(base, base + 1, Bond::SINGLE);
    mol.addBond(base + 1, base + 2, Bond::SINGLE);
    mol.addBond(base + 2, base + 3, Bond::SINGLE);
    mol.addBond(base + 3, base, Bond::SINGLE);
  }

  VECT_INT_VECT rings, bonds;
  for (int ring = 0; ring < num_rings; ++ring) {
    int base = ring * 4;
    rings.push_back({base, base + 1, base + 2, base + 3});
  }

  EXPECT_NO_THROW(RingUtils::convertToBonds(rings, bonds, mol));
  EXPECT_EQ(bonds.size(), num_rings);
  for (const auto &bond_ring : bonds) {
    EXPECT_EQ(bond_ring.size(), 4u);
  }
}

} // namespace
