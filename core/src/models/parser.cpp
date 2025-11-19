#include <array>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <string_view>

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

inline bool is_space_char(char c) {
  return c == ' ' || c == '\t' || c == '\r';
}

template<std::size_t N>
inline std::size_t tokenize_cif_line(const char* begin, const char* end, std::array<std::string_view, N>& tokens) {
  std::size_t count = 0;
  const char* p = begin;
  while (p < end && count < N) {
    while (p < end && is_space_char(*p)) ++p;
    if (p >= end) break;
    if (*p == '#') break;
    const char* token_start = p;
    if (*p == '\'' || *p == '"') {
      const char quote = *p++;
      token_start = p;
      while (p < end && *p != quote) ++p;
      tokens[count++] = std::string_view(token_start, static_cast<std::size_t>(p - token_start));
      if (p < end) ++p;
    } else {
      while (p < end && !is_space_char(*p) && *p != '\n') ++p;
      tokens[count++] = std::string_view(token_start, static_cast<std::size_t>(p - token_start));
    }
  }
  return count;
}

inline int parse_int_sv(std::string_view sv) {
  if (sv.empty()) return 0;
  if (sv[0] == '?' || sv[0] == '.') return 0;
  bool negative = false;
  std::size_t idx = 0;
  if (sv[idx] == '-') {
    negative = true;
    ++idx;
  }
  int value = 0;
  for (; idx < sv.size(); ++idx) {
    char c = sv[idx];
    if (c < '0' || c > '9') break;
    value = value * 10 + (c - '0');
  }
  return negative ? -value : value;
}

inline void assign_dssp_range(std::vector<DSSPAssignment>& assignments, int start_seq_id, int end_seq_id, DSSPAssignment type) {
  if (assignments.empty() || start_seq_id <= 0 || end_seq_id < start_seq_id) return;
  const std::size_t seq_len = assignments.size();
  std::size_t start = static_cast<std::size_t>(start_seq_id - 1);
  if (start >= seq_len) return;
  std::size_t end = static_cast<std::size_t>(end_seq_id - 1);
  if (end >= seq_len) {
    end = seq_len - 1;
  }
  if (end < start) return;
  for (std::size_t i = start; i <= end; ++i) {
    assignments[i] = type;
  }
}

inline DSSPAssignment classify_struct_conf_type(std::string_view type) {
  if (type.empty()) return DSSPAssignment::Coil;
  if (type == "STRN") return DSSPAssignment::Strand;
  if (type == "BEND") return DSSPAssignment::Bend;
  if (type.rfind("TURN", 0) == 0) return DSSPAssignment::Turn;
  if (type.rfind("HELX_LH_PP_P", 0) == 0) return DSSPAssignment::PolyProlineHelix;
  if (type.rfind("HELX_RH_PI_P", 0) == 0) return DSSPAssignment::HelixPi;
  if (type.rfind("HELX_RH_3T_P", 0) == 0) return DSSPAssignment::Helix3_10;
  if (type.rfind("HELX_RH_AL_P", 0) == 0) return DSSPAssignment::AlphaHelix;
  if (type.rfind("HELX", 0) == 0) return DSSPAssignment::AlphaHelix;
  return DSSPAssignment::Coil;
}

inline const char* find_loop_data_start(const char* data, size_t size, const char* first_marker, const char* last_marker) {
  const char* start = std::strstr(data, first_marker);
  if (!start) return nullptr;
  const char* last = std::strstr(start, last_marker);
  if (!last) return nullptr;
  const char* end = data + size;
  const char* newline = static_cast<const char*>(std::memchr(last, '\n', end - last));
  if (!newline) return nullptr;
  ++newline;
  if (newline >= end) return nullptr;
  return newline;
}

inline void parse_struct_conf_entries(const char* data, size_t size, std::vector<DSSPAssignment>& assignments) {
  if (assignments.empty()) return;
  const char* end = data + size;
  const char* line = find_loop_data_start(data, size, "_struct_conf.beg_auth_asym_id", "_struct_conf.pdbx_end_PDB_ins_code");
  if (!line) return;
  while (line < end) {
    const char* trimmed = line;
    while (trimmed < end && is_space_char(*trimmed)) ++trimmed;
    if (trimmed >= end) break;
    if (*trimmed == '#') break;
    if (*trimmed == '\n') {
      line = trimmed + 1;
      continue;
    }
    if ((end - trimmed) >= 5 && std::memcmp(trimmed, "loop_", 5) == 0) break;
    const char* line_end = static_cast<const char*>(std::memchr(trimmed, '\n', end - trimmed));
    if (!line_end) line_end = end;
    std::array<std::string_view, 20> tokens;
    const std::size_t count = tokenize_cif_line(trimmed, line_end, tokens);
    if (count >= 13) {
      const auto type = classify_struct_conf_type(tokens[6]);
      if (type != DSSPAssignment::Coil) {
        const int start_seq = parse_int_sv(tokens[5]);
        const int end_seq   = parse_int_sv(tokens[12]);
        assign_dssp_range(assignments, start_seq, end_seq, type);
      }
    }
    if (line_end >= end) break;
    line = line_end + 1;
  }
}

inline void parse_struct_sheet_range_entries(const char* data, size_t size, std::vector<DSSPAssignment>& assignments) {
  if (assignments.empty()) return;
  const char* end = data + size;
  const char* line = find_loop_data_start(data, size,
                                          "_struct_sheet_range.sheet_id",
                                          "_struct_sheet_range.end_auth_seq_id");
  if (!line) return;
  while (line < end) {
    const char* trimmed = line;
    while (trimmed < end && is_space_char(*trimmed)) ++trimmed;
    if (trimmed >= end) break;
    if (*trimmed == '#') break;
    if (*trimmed == '\n') {
      line = trimmed + 1;
      continue;
    }
    if ((end - trimmed) >= 5 && std::memcmp(trimmed, "loop_", 5) == 0) break;
    const char* line_end = static_cast<const char*>(std::memchr(trimmed, '\n', end - trimmed));
    if (!line_end) line_end = end;
    std::array<std::string_view, 20> tokens;
    const std::size_t count = tokenize_cif_line(trimmed, line_end, tokens);
    if (count >= 9) {
      const int start_seq = parse_int_sv(tokens[4]);
      const int end_seq   = parse_int_sv(tokens[8]);
      assign_dssp_range(assignments, start_seq, end_seq, DSSPAssignment::Strand);
    }
    if (line_end >= end) break;
    line = line_end + 1;
  }
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

inline double parse_fixed_float(const char* start) {
  bool negative = false;
  const char* p = start;
  if (*p == '-') {
    negative = true;
    ++p;
  }
  double value = 0.0;
  while (*p >= '0' && *p <= '9') {
    value = value * 10.0 + static_cast<double>(*p - '0');
    ++p;
  }
  if (*p == '.') {
    ++p;
    double factor = 0.1;
    while (*p >= '0' && *p <= '9') {
      value += static_cast<double>(*p - '0') * factor;
      factor *= 0.1;
      ++p;
    }
  }
  return negative ? -value : value;
}

ModelParserResult parse_model(const char *data, size_t size) {
  ModelParserResult output;
  const char *const end = data + size;
  const char *p = data;

  // difficult to benchmark
  // 34 is much more precise, but is inconsistent and some files fail.
  // p = data + skip_hashes_avx2(data, size, 1);
  output.sequence = extract_marker_value(data, size, "_struct_ref.pdbx_seq_one_letter_code");
  output.metadata.ncbi_taxonomy_id = strip_cif_quotes(
      extract_marker_value(data, size, "_ma_target_ref_db_details.ncbi_taxonomy_id"));
  output.metadata.organism_scientific = strip_cif_quotes(
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

  // find B-factor (pLDDT)
  const char* bfactor_pos = z_pos; // z_pos points to the start of the z value
  while (bfactor_pos < end && *bfactor_pos != ' ') bfactor_pos++; // Skip z value   -> spaces
  while (bfactor_pos < end && *bfactor_pos == ' ') bfactor_pos++; // Skip spaces    -> start of occupancy
  while (bfactor_pos < end && *bfactor_pos != ' ') bfactor_pos++; // Skip occupancy -> spaces
  while (bfactor_pos < end && *bfactor_pos == ' ') bfactor_pos++; // Skip spaces    -> start of B-factor
  int bfactor_offset = bfactor_pos - first_atom;

  // find the end of first line to determine line length
  const char* line_end = static_cast<const char*>(std::memchr(first_atom, '\n', end - first_atom));
  if (!line_end) line_end = end;

  size_t line_length = line_end - first_atom + 1; // +1 to include newline
  if (bfactor_pos >= line_end) bfactor_offset = -1;

  // count ATOM records to reserve space for coordinates
  size_t atom_count = 0;
  const char* counter = first_atom;

  for (const auto &c : output.sequence) {
    atom_count += StandardAminoAcidAtomSizeTable[static_cast<unsigned char>(c)];
  }
  atom_count += 1; // for Frodo

  output.coords.reserve(atom_count);
  output.plddt_per_residue.assign(output.sequence.size(), pLDDTCategory::VeryLow);
  output.dssp_per_residue .assign(output.sequence.size(), DSSPAssignment::Coil);
  parse_struct_conf_entries(data, size, output.dssp_per_residue);
  parse_struct_sheet_range_entries(data, size, output.dssp_per_residue);

  size_t residue_index = 0;
  size_t atoms_in_residue = 0;
  if (!output.sequence.empty()) {
    atoms_in_residue = StandardAminoAcidAtomSizeTable[static_cast<unsigned char>(output.sequence[0])];
    if (atoms_in_residue == 0) atoms_in_residue = 1;
  }
  size_t atoms_seen_in_residue = 0;
  bool recorded_plddt = false;

  // 1. all atom recods are the same length
  // 2. atom records are consecutive
  // 3. positions have consistent offsets within the file
  const char *atom_ptr = first_atom;
  for (size_t i = 0; i < atom_count; i++) {
    if (atom_ptr + 4 >= end || atom_ptr[0] != 'A' || atom_ptr[1] != 'T' || atom_ptr[2] != 'O' || atom_ptr[3] != 'M') {
      break;
    }

    // x position
    double x = 0.0;
    bool x_negative = false;
    const char *x_start = atom_ptr + x_offset;

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

    // y position
    double y = 0.0;
    bool y_negative = false;
    const char *y_start = atom_ptr + y_offset;

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

    // z position
    double z = 0.0;
    bool z_negative = false;
    const char *z_start = atom_ptr + z_offset;

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

    if (!recorded_plddt && bfactor_offset >= 0 && residue_index < output.plddt_per_residue.size()) {
      const char *b_start = atom_ptr + bfactor_offset;
      if (b_start < atom_ptr + line_length) {
        double bfactor = parse_fixed_float(b_start);
        auto category = categorize_plddt(bfactor);
        output.plddt_per_residue[residue_index] = category;
        recorded_plddt = true;
      }
    }

    if (residue_index < output.plddt_per_residue.size() && atoms_in_residue > 0) {
      atoms_seen_in_residue++;
      if (atoms_seen_in_residue >= atoms_in_residue) {
        residue_index++;
        atoms_seen_in_residue = 0;
        recorded_plddt = false;
        if (residue_index < output.sequence.size()) {
          atoms_in_residue =
              StandardAminoAcidAtomSizeTable[static_cast<unsigned char>(output.sequence[residue_index])];
          if (atoms_in_residue == 0) atoms_in_residue = 1;
        } else {
          atoms_in_residue = 0;
        }
      }
    }

    atom_ptr += line_length;
  }

  return output;
}

} // namespace lahuta
