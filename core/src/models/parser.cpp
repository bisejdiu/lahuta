#include <cstddef>
#include <cstdlib>
#include <cstring>

#include "models/fast_lookup.hpp"
#include "models/parser.hpp"
#include "simde/x86/avx2.h"

// clang-format off
namespace lahuta {

//
// Returns the offset where the required_count-th '#' appears,
// or 'size' if fewer than required_count '#' characters exist.
// It is not that efficient (e.g., '#' always appear at the start of the line). Nevertheless,
// this is already too deep in the microbenchmark region to be worth spending more time on.
//
inline size_t skip_hashes_avx2(const char* data, size_t size, int required_count) {

  const size_t simd_width = 32;
  size_t offset = 0;
  int hash_count = 0;

  // Create a 32-byte vector filled with '#' (0x23).
  simde__m256i hashVec = simde_mm256_set1_epi8('#');

  while ((offset + simd_width) <= size) {
    simde__m256i chunk = simde_mm256_loadu_si256(reinterpret_cast<const simde__m256i*>(data + offset));

    simde__m256i cmp = simde_mm256_cmpeq_epi8(chunk, hashVec);
    uint32_t mask = static_cast<uint32_t>(simde_mm256_movemask_epi8(cmp));

    if (mask == 0) {
      offset += simd_width;
      continue;
    }

    // There is at least one '#' in this chunk.
    // Iterate over each byte in the chunk (from least significant bit to most).
    for (int i = 0; i < 32; i++) {
      if (mask & 1) {
        hash_count++;
        if (hash_count == required_count) {
          // Found the required '#' at position (offset + i).
          return offset + i;
        }
      }
      mask >>= 1; // Move to the next byte's bit.
    }
    offset += simd_width;
  }

  while (offset < size) {
    if (data[offset] == '#') {
      hash_count++;
      if (hash_count == required_count) {
        return offset;
      }
    }
    offset++;
  }

  return size;
}


inline auto extract_marker_value(const char* data, size_t size, const char* marker) {
  std::string result;

  const char* marker_pos = std::strstr(data, marker);
  if (!marker_pos) return result;

  // move pointer to after the marker.
  const char* p = marker_pos + std::strlen(marker);

  const char* q = p;
  while (q < data + size && (*q == ' ' || *q == '\t')) { // first non-space character.
    ++q;
  }

  if (q < data + size && *q != '\n') { // value is inline
    const char* end_line = static_cast<const char*>(std::memchr(q, '\n', (data + size) - q));
    if (!end_line) {
      end_line = data + size;
    }
    size_t length = end_line - q;
    return std::string(q, length);
  }

  const char* newline = static_cast<const char*>(std::memchr(p, '\n', (data + size) - p));
  if (!newline) return result;

  p = newline + 1;
  if (p >= data + size || *p != ';') {
    return result;
  }
  p++;

  // ';' enclosed multi-line value.
  result.reserve(345);
  while (p < data + size) {
    const char* line_end = static_cast<const char*>(std::memchr(p, '\n', (data + size) - p));
    if (!line_end) {
      result.append(p, (data + size) - p);
      break;
    }

    result.append(p, line_end - p);
    p = line_end + 1;

    if (p < data + size && *p == ';') break;
  }

  return result;
}

inline std::string strip_cif_quotes(std::string value) {
  if (value.size() >= 2) {
    const char first = value.front();
    const char last  = value.back();
    if ((first == '"' && last == '"') || (first == '\'' && last == '\'')) {
      value.erase(value.size() - 1);
      value.erase(value.begin());
    }
  }
  return value;
}

ModelParserResult parse_model(const char *data, size_t size) {
  ModelParserResult output;
  const char *const end = data + size;
  const char *p = data;

  // difficult to benchmark
  // 34 is much more precise, but is inconsistent and some files fail.
  // p = data + skip_hashes_avx2(data, size, 1);
  output.sequence = extract_marker_value(data, size, "_struct_ref.pdbx_seq_one_letter_code");
  output.ncbi_taxonomy_id = strip_cif_quotes(
      extract_marker_value(data, size, "_ma_target_ref_db_details.ncbi_taxonomy_id"));
  output.organism_scientific = strip_cif_quotes(
      extract_marker_value(data, size, "_ma_target_ref_db_details.organism_scientific"));

  // we got the sequence, now we skip to the first ATOM record
  while (p + 4 < end) {
    if (p[0] == 'A' && p[1] == 'T' && p[2] == 'O' && p[3] == 'M' && p[4] == ' ') break;
    p = static_cast<const char *>(std::memchr(p, '\n', end - p));
    if (!p) return output;
    p++;
  }

  if (p + 4 >= end) return output;

  // we need to detect column positions from first ATOM line
  const char* first_atom = p;

  // find the '?' character which precedes the coordinates
  const char* q_mark = static_cast<const char*>(std::memchr(first_atom, '?', end - first_atom));
  if (!q_mark) return {};

  // skip '?' and any spaces after it
  const char* coord_start = q_mark + 1;
  while (coord_start < end && *coord_start == ' ') coord_start++;

  // find x, y, and z column offsets
  int x_offset = coord_start - first_atom;

  // find y
  const char* y_pos = coord_start;
  while (y_pos < end && ((*y_pos >= '0' && *y_pos <= '9') || *y_pos == '-' || *y_pos == '.')) y_pos++; // xkip x value (digits, minus sign, period)
  while (y_pos < end && *y_pos == ' ') y_pos++; // skip whitespace to get to y start
  int y_offset = y_pos - first_atom;

  // find z
  const char* z_pos = y_pos;
  while (z_pos < end && ((*z_pos >= '0' && *z_pos <= '9') || *z_pos == '-' || *z_pos == '.')) z_pos++;
  while (z_pos < end && *z_pos == ' ') z_pos++;
  int z_offset = z_pos - first_atom;

  // find the end of first line to determine line length
  const char* line_end = static_cast<const char*>(std::memchr(first_atom, '\n', end - first_atom));
  if (!line_end) line_end = end;

  size_t line_length = line_end - first_atom + 1; // +1 to include newline

  // count ATOM records to reserve space for coordinates
  size_t atom_count = 0;
  const char* counter = first_atom;

  for (const auto &c : output.sequence) {
    atom_count += StandardAminoAcidAtomSizeTable[c];
  }
  atom_count += 1; // for Frodo

  output.coords.reserve(atom_count);

  // 1. all atom recods are the same length
  // 2. atom records are consecutive
  // 3. positions have consistent offsets within the file
  const char* atom_ptr = first_atom;
  for (size_t i = 0; i < atom_count; i++) {
      if (atom_ptr + 4 >= end ||  atom_ptr[0] != 'A' || atom_ptr[1] != 'T' || atom_ptr[2] != 'O' || atom_ptr[3] != 'M') {
          break;
      }

      // --- x position ---
      double x = 0.0;
      bool x_negative = false;
      const char* x_start = atom_ptr + x_offset;

      // check sign
      if (*x_start == '-') {
          x_negative = true;
          x_start++;
      }

      // parse integer part
      while (*x_start >= '0' && *x_start <= '9') {
          x = x * 10.0 + (*x_start - '0');
          x_start++;
      }

      // skip '.'
      if (*x_start == '.') x_start++;

      // parse decimal part
      double fraction = 0.1;
      while (*x_start >= '0' && *x_start <= '9') {
          x += (*x_start - '0') * fraction;
          fraction *= 0.1;
          x_start++;
      }
      if (x_negative) x = -x;

      // --- y position ---
      double y = 0.0;
      bool y_negative = false;
      const char* y_start = atom_ptr + y_offset;

      if (*y_start == '-') {
          y_negative = true;
          y_start++;
      }

      while (*y_start >= '0' && *y_start <= '9') {
          y = y * 10.0 + (*y_start - '0');
          y_start++;
      }

      if (*y_start == '.') y_start++;

      fraction = 0.1;
      while (*y_start >= '0' && *y_start <= '9') {
          y += (*y_start - '0') * fraction;
          fraction *= 0.1;
          y_start++;
      }

      if (y_negative) y = -y;

      // --- z position ---
      double z = 0.0;
      bool z_negative = false;
      const char* z_start = atom_ptr + z_offset;

      if (*z_start == '-') {
          z_negative = true;
          z_start++;
      }

      while (*z_start >= '0' && *z_start <= '9') {
          z = z * 10.0 + (*z_start - '0');
          z_start++;
      }

      if (*z_start == '.') z_start++;

      fraction = 0.1;
      while (*z_start >= '0' && *z_start <= '9') {
          z += (*z_start - '0') * fraction;
          fraction *= 0.1;
          z_start++;
      }

      if (z_negative) z = -z;

      output.coords.push_back({x, y, z});

      atom_ptr += line_length;
  }

  return output;
}

} // namespace lahuta
