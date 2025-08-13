#ifndef LAHUTA_DB_DATA_TYPE_HPP
#define LAHUTA_DB_DATA_TYPE_HPP

#include <cstdint>

namespace lahuta {

// A simple binary header for the serialized data.
struct SerializedModelData {
  std::uint32_t sequence_length;
  std::uint32_t num_points;

  const char *sequence_data() const {
    return reinterpret_cast<const char *>(this) + sizeof(SerializedModelData);
  }

  const float *coords_data_float() const {
    return reinterpret_cast<const float *>(
        reinterpret_cast<const char *>(this) + sizeof(SerializedModelData) + sequence_length);
  }
};

} // namespace lahuta

#endif // LAHUTA_DB_DATA_TYPE_HPP
