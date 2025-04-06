#ifndef LAHUTA_DB_READER_HPP
#define LAHUTA_DB_READER_HPP

#include "db/data_type.hpp"
#include "lmdb/lmdb++.h"
#include "models/parser.hpp"

// clang-format off
namespace lahuta {

class LMDBReader {
public:
  LMDBReader(lmdb::env &env, lmdb::dbi &dbi) : m_env(env), m_dbi(dbi) {}

  /// Retrieve a ModelParserResult based on a key. Returns false if key not found.
  bool retrieve(const std::string &key, ModelParserResult &result) {
    auto txn = lmdb::txn::begin(m_env.handle(), nullptr, MDB_RDONLY);
    std::string_view value;
    if (!m_dbi.get(txn.handle(), key, value)) { return false; }

    const SerializedModelData* serialized = reinterpret_cast<const SerializedModelData*>(value.data());

    // basic sanity check on the size.
    size_t expected_size = sizeof(SerializedModelData) +
                           serialized->sequence_length +
                           serialized->num_points * 3  *
                           sizeof(float); // 3 floats per point

    if (value.size() < expected_size) {
      throw std::runtime_error("Corrupted data for key: " + key);
    }

    result.sequence.assign(serialized->sequence_data(), serialized->sequence_length);

    // float coords --> Point3D (double)
    result.coords.resize(serialized->num_points);
    const float* float_coords = serialized->coords_data_float();
    for (uint32_t i = 0; i < serialized->num_points; ++i) {
      result.coords[i].x = static_cast<double>(float_coords[i * 3]);
      result.coords[i].y = static_cast<double>(float_coords[i * 3 + 1]);
      result.coords[i].z = static_cast<double>(float_coords[i * 3 + 2]);
    }

    return true; // Read-only transactions are automatically rolled back.
  }

private:
  lmdb::env &m_env;
  lmdb::dbi &m_dbi;
};

} // namespace lahuta

#endif // LAHUTA_DB_READER_HPP

