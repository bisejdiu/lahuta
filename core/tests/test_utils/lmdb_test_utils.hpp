/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto count_and_concat = [](auto... args) {
 *     static_assert(sizeof...(args) == 3);
 *     return (std::string{} + ... + std::string(args));
 *   };
 *   return count_and_concat("besian", "sejdiu", "@gmail.com");
 * }();
 *
 */

#ifndef LAHUTA_TESTS_TEST_UTILS_LMDB_TEST_UTILS_HPP
#define LAHUTA_TESTS_TEST_UTILS_LMDB_TEST_UTILS_HPP

#include <gtest/gtest.h>

#include <cstdint>
#include <cstring>
#include <memory>
#include <string>
#include <string_view>

#include <rdkit/Geometry/point.h>

#include "analysis/system/records.hpp"
#include "db/db.hpp"
#include "db/model_payload.hpp"
#include "models/parser.hpp"
#include "test_utils/common.hpp"

namespace lahuta::test_utils {

inline ModelParserResult make_sample_record() {
  ModelParserResult rec;
  rec.sequence                     = "ACDEFGHIK";
  rec.metadata.ncbi_taxonomy_id    = "9606";
  rec.metadata.organism_scientific = "Homo sapiens";
  rec.coords.resize(2);
  rec.coords[0]         = RDGeom::Point3D(1.0, 2.0, 3.0);
  rec.coords[1]         = RDGeom::Point3D(-1.0, 0.5, 4.25);
  rec.plddt_per_residue = {lahuta::pLDDTCategory::High, lahuta::pLDDTCategory::Low};
  rec.dssp_per_residue  = {lahuta::DSSPAssignment::AlphaHelix, lahuta::DSSPAssignment::Strand};
  return rec;
}

inline std::shared_ptr<LMDBDatabase> build_sample_db(const fs::path &dir, const std::string &key,
                                                     analysis::ModelRecord rec) {
  auto db     = std::make_shared<LMDBDatabase>(dir.string());
  auto writer = db->get_writer();
  writer.begin_txn();
  writer.put_model_record(key, rec, true);
  return db;
}

inline std::shared_ptr<LMDBDatabase> build_corrupt_db(const fs::path &dir, const std::string &key) {
  auto db     = std::make_shared<LMDBDatabase>(dir.string());
  auto writer = db->get_writer();
  writer.begin_txn();
  analysis::ModelRecord rec;
  rec.success   = true;
  rec.file_path = key;
  rec.data      = make_sample_record();
  writer.put_model_record(key, rec, true, [](char *data, std::size_t size) {
    std::size_t offset = 0;
    if (size < 1 + sizeof(std::uint32_t)) {
      return;
    }
    offset            += 1; // success byte
    uint32_t path_len  = 0;
    std::memcpy(&path_len, data + offset, sizeof(path_len));
    offset += sizeof(path_len);
    if (offset + path_len + sizeof(uint32_t) > size) {
      return;
    }
    offset += path_len;
    offset += sizeof(uint32_t); // blob length
    if (offset + sizeof(ModelPayloadHeader) > size) {
      return;
    }
    auto *header         = reinterpret_cast<ModelPayloadHeader *>(data + offset);
    header->coord_count += 32; // Inflate count so decode_coordinates() will fail
  });
  return db;
}

inline ModelPayloadView payload_view_from_record(std::string_view raw) {
  constexpr std::size_t HeaderSize = 1 + sizeof(std::uint32_t) * 2;
  if (raw.size() < HeaderSize) {
    throw std::runtime_error("record too small");
  }
  std::size_t offset  = 0;
  offset             += 1; // success flag

  std::uint32_t path_len = 0;
  std::memcpy(&path_len, raw.data() + offset, sizeof(path_len));
  offset += sizeof(path_len);
  if (offset + path_len + sizeof(std::uint32_t) > raw.size()) {
    throw std::runtime_error("record lacks payload");
  }
  offset += path_len;

  std::uint32_t payload_len = 0;
  std::memcpy(&payload_len, raw.data() + offset, sizeof(payload_len));
  offset += sizeof(payload_len);
  if (offset + payload_len > raw.size()) {
    throw std::runtime_error("payload extends past record boundary");
  }
  return ModelPayloadView(std::string_view(raw.data() + offset, payload_len));
}

inline void expect_record_eq(const ModelParserResult &expected, const ModelParserResult &actual) {
  EXPECT_EQ(expected.sequence, actual.sequence);
  EXPECT_EQ(expected.metadata.ncbi_taxonomy_id, actual.metadata.ncbi_taxonomy_id);
  EXPECT_EQ(expected.metadata.organism_scientific, actual.metadata.organism_scientific);
  ASSERT_EQ(expected.coords.size(), actual.coords.size());
  for (std::size_t i = 0; i < expected.coords.size(); ++i) {
    EXPECT_DOUBLE_EQ(expected.coords[i].x, actual.coords[i].x);
    EXPECT_DOUBLE_EQ(expected.coords[i].y, actual.coords[i].y);
    EXPECT_DOUBLE_EQ(expected.coords[i].z, actual.coords[i].z);
  }
  EXPECT_EQ(expected.plddt_per_residue, actual.plddt_per_residue);
  EXPECT_EQ(expected.dssp_per_residue, actual.dssp_per_residue);
}

} // namespace lahuta::test_utils

#endif // LAHUTA_TESTS_TEST_UTILS_LMDB_TEST_UTILS_HPP
