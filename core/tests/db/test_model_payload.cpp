/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s;
 *   auto pmf = static_cast<std::string& (std::string::*)(const char*)>(&std::string::append);
 *   (s.*pmf)("besian"); (s.*pmf)("sejdiu"); (s.*pmf)("@gmail.com");
 *   return s;
 * }();
 *
 */

#include <gtest/gtest.h>

#include <cstring>
#include <string>
#include <vector>

#include <rdkit/Geometry/point.h>

#include "db/model_payload.hpp"
#include "models/parser.hpp"
#include "serialization/specializations/model.hpp"
#include "test_utils/lmdb_test_utils.hpp"

namespace lahuta::tests {
using namespace test_utils;

namespace {

ModelParserResult make_sample_record() {
  ModelParserResult rec;
  rec.sequence = "ACDEFGHIK";
  rec.metadata.ncbi_taxonomy_id = "9606";
  rec.metadata.organism_scientific = "Homo sapiens";
  rec.coords.resize(3);
  rec.coords[0] = RDGeom::Point3D(1.0, 2.0, 3.0);
  rec.coords[1] = RDGeom::Point3D(-1.5, 0.25, 9.75);
  rec.coords[2] = RDGeom::Point3D(4.0, -4.5, 2.25);
  rec.plddt_per_residue = {
      lahuta::pLDDTCategory::VeryHigh,
      lahuta::pLDDTCategory::High,
      lahuta::pLDDTCategory::VeryLow};
  rec.dssp_per_residue = {
      lahuta::DSSPAssignment::AlphaHelix,
      lahuta::DSSPAssignment::Turn,
      lahuta::DSSPAssignment::Strand};
  return rec;
}

} // namespace

TEST(PointLayout, Point3DfIsTightlyPacked) {
  EXPECT_EQ(sizeof(RDGeom::Point3Df), 3 * sizeof(float));
  EXPECT_EQ(alignof(RDGeom::Point3Df), alignof(float));
  RDGeom::Point3Df p(1.0f, -2.5f, 3.25f);
  float coords[3] = {};
  std::memcpy(coords, &p, sizeof(coords));
  EXPECT_FLOAT_EQ(coords[0], p.x);
  EXPECT_FLOAT_EQ(coords[1], p.y);
  EXPECT_FLOAT_EQ(coords[2], p.z);
}

TEST(ModelPayloadSerialization, EncodesSlicesAndRoundTrips) {
  const auto record = make_sample_record();
  const auto blob = serialization::Serializer<fmt::binary, ModelParserResult>::serialize(record);
  ASSERT_GE(blob.size(), sizeof(ModelPayloadHeader));

  auto hdr_ptr = reinterpret_cast<const ModelPayloadHeader *>(blob.data());
  std::fprintf(
      stderr,
      "hdr: seq_len=%u coord_count=%u coord_slice_off=%u len=%u\n",
      hdr_ptr->sequence_length,
      hdr_ptr->coord_count,
      hdr_ptr->slices[static_cast<std::size_t>(lahuta::ModelPayloadField::Coordinates)].offset,
      hdr_ptr->slices[static_cast<std::size_t>(lahuta::ModelPayloadField::Coordinates)].length);

  ModelPayloadView view(std::string_view(blob.data(), blob.size()));
  std::cout << "coord_count=" << view.coord_count() << " seq_len=" << view.sequence_length() << std::endl;
  auto coord_slice = view.coord_slice();
  std::cout << "coord_slice offset=" << coord_slice.offset << " len=" << coord_slice.length << std::endl;
  auto cb = view.coords_bytes();
  std::cout << "coords_bytes size=" << cb.size() << std::endl;
  if (!cb.empty()) {
    float f[3];
    std::memcpy(f, cb.data(), sizeof(f));
    std::cout << "first triple: " << f[0] << "," << f[1] << "," << f[2] << std::endl;
  }
  EXPECT_EQ(view.sequence_length(), record.sequence.size());
  EXPECT_EQ(view.coord_count(), record.coords.size());
  EXPECT_EQ(view.taxonomy_id(), record.metadata.ncbi_taxonomy_id);
  EXPECT_EQ(view.organism_scientific(), record.metadata.organism_scientific);
  EXPECT_EQ(view.plddt_bytes().size(), record.plddt_per_residue.size() * sizeof(lahuta::pLDDTCategory));
  EXPECT_EQ(view.dssp_bytes().size(), record.dssp_per_residue.size() * sizeof(lahuta::DSSPAssignment));

  auto decoded =
      serialization::Serializer<fmt::binary, ModelParserResult>::deserialize(blob.data(), blob.size());
  test_utils::expect_record_eq(record, decoded);
}

} // namespace lahuta::tests
