/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto f = [](auto&&... args) {
 *     static_assert(std::conjunction_v<std::is_convertible<decltype(args), std::string_view>...>);
 *     return (std::string{} + ... + std::string(args));
 *   };
 *   return f("besian", "sejdiu", "@gmail.com");
 * }();
 *
 */

#include <cctype>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/Conformer.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <gemmi/third_party/fast_float.h>

#include "chemistry/elements.hpp"
#include "md/element_utils.hpp"
#include "md/gro_reader.hpp"
#include "mmap/MemoryMapped.h"
#include "residues/definitions.hpp"

// clang-format off
namespace lahuta::md {
namespace {

[[nodiscard]] inline bool next_line(const char *data, std::size_t size, std::size_t &pos, std::string_view &out) {
  if (pos >= size) return false;
  const char *begin = data + pos;
  const char *end   = data + size;
  const void *nl_ptr = std::memchr(begin, '\n', static_cast<std::size_t>(end - begin));
  if (nl_ptr) {
    const char *nl  = static_cast<const char *>(nl_ptr);
    std::size_t len = static_cast<std::size_t>(nl - begin);
    // Handle optional CR before LF
    if (len > 0 && begin[len - 1] == '\r') --len;
    out = std::string_view(begin, len);
    pos = static_cast<std::size_t>(nl - data) + 1; // move past '\n'
    return true;
  } else {
    // Last line without trailing newline
    std::size_t len = static_cast<std::size_t>(end - begin);
    if (len > 0 && begin[len - 1] == '\r') --len;
    out = std::string_view(begin, len);
    pos = size;
    return true;
  }
}

inline std::string trim_copy(std::string_view s) {
  const std::size_t first = s.find_first_not_of(" \t\r\n");
  if (first == std::string_view::npos) return std::string();
  const std::size_t last = s.find_last_not_of(" \t\r\n");
  return std::string{s.substr(first, last - first + 1)};
}

[[nodiscard]] inline bool is_ascii_space(char c) noexcept {
  // Matches the set used in ws: " \t\r\n\f\v"
  return c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\f' || c == '\v';
}

inline std::string trim_span(std::string_view str, std::size_t start, std::size_t length) {
  const std::size_t end = std::min(start + length, str.size());
  if (start >= end) return std::string();
  return trim_copy(str.substr(start, end - start));
}

[[nodiscard]] inline std::optional<int> parse_int(std::string_view line, std::size_t start, std::size_t width) {
  if (start >= line.size() || width == 0) return std::nullopt;
  const std::size_t end = std::min(start + width, line.size());
  const char *begin = line.data() + start;
  const char *last  = line.data() + end;

  // Trim leading whitespace
  while (begin < last && is_ascii_space(*begin)) ++begin;
  if (begin >= last) return std::nullopt;

  // Optional sign
  bool negative = false;
  if (*begin == '+' || *begin == '-') {
    negative = (*begin == '-');
    ++begin;
  }
  if (begin >= last) return std::nullopt;

  // Parse digits
  int value = 0;
  const char *p = begin;
  if (*p < '0' || *p > '9') return std::nullopt;
  while (p < last && *p >= '0' && *p <= '9') {
    value = value * 10 + (*p - '0');
    ++p;
  }
  // Remaining must be whitespace
  while (p < last && is_ascii_space(*p)) ++p;
  if (p != last) return std::nullopt;

  return negative ? -value : value;
}

[[nodiscard]] inline std::optional<unsigned> detect_field_width_ddist(std::string_view line, std::string *error_msg = nullptr) {
  //
  // Follows the GROMACS rule: Determine spacing between the first three decimal points in the line.
  // ddist = distance between consecutive '.' characters must be consistent.
  //
  const std::size_t p1 = line.find('.');
  if (p1 == std::string_view::npos) {
    if (error_msg) *error_msg = "gro_reader: a coordinate line does not contain a '.'";
    return std::nullopt;
  }
  const std::size_t p2 = line.find('.', p1 + 1);
  if (p2 == std::string_view::npos) {
    if (error_msg) *error_msg = "gro_reader: a coordinate line does not contain a '.'";
    return std::nullopt;
  }
  const std::size_t p3 = line.find('.', p2 + 1);
  if (p3 == std::string_view::npos) {
    if (error_msg) *error_msg = "gro_reader: a coordinate line does not contain a '.'";
    return std::nullopt;
  }
  const std::ptrdiff_t ddist = static_cast<std::ptrdiff_t>(p2) - static_cast<std::ptrdiff_t>(p1);
  if (static_cast<std::ptrdiff_t>(p3) - static_cast<std::ptrdiff_t>(p2) != ddist) {
    if (error_msg) *error_msg = "gro_reader: inconsistent spacing of decimal points for x, y, z";
    return std::nullopt;
  }
  if (ddist <= 0) {
    if (error_msg) *error_msg = "gro_reader: invalid decimal point spacing";
    return std::nullopt;
  }
  return static_cast<unsigned>(ddist);
}

[[nodiscard]] inline std::optional<double>
parse_fixed_width_double(std::string_view line, std::size_t start, std::size_t width, std::string *error_msg = nullptr) {
  if (start >= line.size() || width == 0) return std::nullopt;

  const std::size_t end = std::min(start + width, line.size());
  const char *begin = line.data() + start;
  const char *last  = line.data() + end;

  // Trim leading ASCII whitespace
  while (begin < last && is_ascii_space(*begin)) ++begin;
  if (begin >= last) return std::nullopt; // all whitespace

  double value{};
  auto [ptr, ec] = fast_float::from_chars(begin, last, value);
  if (ec == std::errc::invalid_argument) return std::nullopt; // nothing parsed
  if (ec == std::errc::result_out_of_range) {
    if (error_msg) *error_msg = "gro_reader: value out of range";
    return std::nullopt;
  }

  // Only allow trailing ASCII whitespace
  while (ptr < last && is_ascii_space(*ptr)) ++ptr;
  if (ptr != last) {
    if (error_msg) *error_msg = "gro_reader: malformed fixed-width coordinate field (extra characters)";
    return std::nullopt;
  }

  return value;
}

struct AtomFixedFields {
  int resid;
  std::string resname;
  std::string atom_name;
  int atom_id;
};

inline void add_rdkit_atom_to_mol(RDKit::RWMol &mol, Element el, const std::string &atom_name, const std::string &resname, int atom_id, int resid) {
  auto *rd_atom = new RDKit::Atom(static_cast<int>(atomic_number_from_element(el)));
  if (el == Element::D) rd_atom->setIsotope(2);

  auto *info = new RDKit::AtomPDBResidueInfo(atom_name, atom_id, "", resname, resid, "");

  const bool is_protein = definitions::is_protein_extended(resname);
  info->setIsHeteroAtom(!is_protein);
  info->setMonomerType(RDKit::AtomMonomerInfo::PDBRESIDUE);
  rd_atom->setMonomerInfo(info);
  mol.addAtom(rd_atom, true, true);
}

inline void attach_conformer(RDKit::RWMol &mol, std::vector<RDGeom::Point3D> &&coords, int natoms) {
  auto *conformer = new RDKit::Conformer(static_cast<unsigned int>(natoms));
  conformer->setAllAtomPositions(std::move(coords));
  conformer->setId(0);
  mol.addConformer(conformer, true);
  mol.updatePropertyCache(false);
}

} // namespace

static std::shared_ptr<RDKit::RWMol> read_gro_mmap_to_rwmol(const MemoryMapped &mm) {
  const unsigned char *raw = mm.getData();
  if (raw == nullptr || mm.size() == 0) {
    throw std::runtime_error("gro_reader: cannot open file: mmap empty or invalid");
  }

  const char *data = reinterpret_cast<const char *>(raw);
  const std::size_t size = static_cast<std::size_t>(mm.size());
  std::size_t pos = 0;

  std::string_view title;
  if (!next_line(data, size, pos, title)) {
    throw std::runtime_error("gro_reader: missing title line");
  }

  std::string_view natoms_sv;
  if (!next_line(data, size, pos, natoms_sv)) {
    throw std::runtime_error("gro_reader: missing atom count line");
  }
  int natoms = 0;
  {
    std::string natoms_line = trim_copy(natoms_sv);
    std::istringstream iss(natoms_line);
    if (!(iss >> natoms) || natoms < 0) {
      throw std::runtime_error("gro_reader: invalid atom count: " + natoms_line);
    }
  }

  auto mol = std::make_shared<RDKit::RWMol>();
  std::vector<RDGeom::Point3D> coords;
  coords.reserve(static_cast<std::size_t>(natoms));

  unsigned ddist = 0;
  for (int i = 0; i < natoms; ++i) {
    std::string_view line;
    if (!next_line(data, size, pos, line)) {
      throw std::runtime_error("gro_reader: unexpected EOF while reading atoms");
    }

    if (i == 0) {
      if (line.size() < 39) throw std::runtime_error("gro_reader: invalid atom line (too short)");
      std::string ddist_err;
      auto ddist_opt = detect_field_width_ddist(line, &ddist_err);
      if (!ddist_opt.has_value()) throw std::runtime_error(ddist_err);
      ddist = *ddist_opt;
    }

    if (line.size() < 39) throw std::runtime_error(std::string("gro_reader: invalid atom line at atom ") + std::to_string(i + 1));

    AtomFixedFields fields{
        parse_int(line, 0,  5).value_or(0),
        trim_span(line, 5,  5),
        trim_span(line, 10, 5),
        parse_int(line, 15, 5).value_or(0)};

    std::string coord_err;
    auto x_opt = parse_fixed_width_double(line, 20 + static_cast<std::size_t>(0) * ddist, ddist, &coord_err);
    auto y_opt = parse_fixed_width_double(line, 20 + static_cast<std::size_t>(1) * ddist, ddist, &coord_err);
    auto z_opt = parse_fixed_width_double(line, 20 + static_cast<std::size_t>(2) * ddist, ddist, &coord_err);

    if (!x_opt.has_value() || !y_opt.has_value() || !z_opt.has_value()) {
      if (!coord_err.empty()) throw std::runtime_error(coord_err);
      throw std::runtime_error("gro_reader: malformed coordinate fields in atom line");
    }

    constexpr double angstrom_per_nm = 10.0;
    RDGeom::Point3D pt(*x_opt * angstrom_per_nm, *y_opt * angstrom_per_nm, *z_opt * angstrom_per_nm);
    coords.push_back(pt);

    const Element el = guess_element_from_gro(fields.atom_name, fields.resname);
    add_rdkit_atom_to_mol(*mol, el, fields.atom_name, fields.resname, fields.atom_id, fields.resid);
  }

  std::string_view box_line;
  if (!next_line(data, size, pos, box_line)) {} // No box

  if (static_cast<int>(coords.size()) != natoms || static_cast<int>(mol->getNumAtoms()) != natoms) {
    throw std::runtime_error("gro_reader: atom count mismatch after parsing");
  }

  attach_conformer(*mol, std::move(coords), natoms);
  return mol;
}

std::shared_ptr<RDKit::RWMol> read_gro_stream_to_rwmol(std::istream &stream) {
  std::string title;
  if (!std::getline(stream, title)) {
    throw std::runtime_error("gro_reader: missing title line");
  }

  std::string natoms_line;
  if (!std::getline(stream, natoms_line)) {
    throw std::runtime_error("gro_reader: missing atom count line");
  }
  int natoms = 0;
  {
    std::istringstream iss(natoms_line);
    if (!(iss >> natoms) || natoms < 0) {
      throw std::runtime_error("gro_reader: invalid atom count: " + natoms_line);
    }
  }

  auto mol = std::make_shared<RDKit::RWMol>();
  std::vector<RDGeom::Point3D> coords;
  coords.reserve(static_cast<std::size_t>(natoms));

  unsigned ddist = 0;
  for (int i = 0; i < natoms; ++i) {
    std::string line;
    if (!std::getline(stream, line)) throw std::runtime_error("gro_reader: unexpected EOF while reading atoms");

    if (i == 0) {
      if (line.size() < 39) throw std::runtime_error("gro_reader: invalid atom line (too short)");

      std::string ddist_err;
      auto ddist_opt = detect_field_width_ddist(line, &ddist_err);
      if (!ddist_opt.has_value()) throw std::runtime_error(ddist_err);
      ddist = *ddist_opt;
    }

    if (line.size() < 39) throw std::runtime_error(std::string("gro_reader: invalid atom line at atom ") + std::to_string(i + 1));

    AtomFixedFields fields{
        parse_int(line, 0,  5).value_or(0),
        trim_span(line, 5,  5),
        trim_span(line, 10, 5),
        parse_int(line, 15, 5).value_or(0)};

    std::string coord_err;
    auto x_opt = parse_fixed_width_double(line, 20 + static_cast<std::size_t>(0) * ddist, ddist, &coord_err);
    auto y_opt = parse_fixed_width_double(line, 20 + static_cast<std::size_t>(1) * ddist, ddist, &coord_err);
    auto z_opt = parse_fixed_width_double(line, 20 + static_cast<std::size_t>(2) * ddist, ddist, &coord_err);

    if (!x_opt.has_value() || !y_opt.has_value() || !z_opt.has_value()) {
      if (!coord_err.empty()) throw std::runtime_error(coord_err);
      throw std::runtime_error("gro_reader: malformed coordinate fields in atom line");
    }

    constexpr double angstrom_per_nm = 10.0;
    RDGeom::Point3D pt(*x_opt * angstrom_per_nm, *y_opt * angstrom_per_nm, *z_opt * angstrom_per_nm);
    coords.push_back(pt);

    const Element el = guess_element_from_gro(fields.atom_name, fields.resname);
    add_rdkit_atom_to_mol(*mol, el, fields.atom_name, fields.resname, fields.atom_id, fields.resid);
  }

  // Try to read box line (3 or 9 floats).
  std::string box_line; std::getline(stream, box_line);

  if (static_cast<int>(coords.size()) != natoms || static_cast<int>(mol->getNumAtoms()) != natoms) {
    throw std::runtime_error("gro_reader: atom count mismatch after parsing");
  }

  attach_conformer(*mol, std::move(coords), natoms);
  return mol;
}

std::shared_ptr<RDKit::RWMol> read_gro_to_rwmol(std::string_view filename) {
  {
    MemoryMapped mm;
    if (mm.open(std::string(filename), MemoryMapped::WholeFile, MemoryMapped::SequentialScan) && mm.isValid()) {
      try {
        return read_gro_mmap_to_rwmol(mm);
      } catch (...) {}
    }
  }

  std::ifstream in{std::string(filename)};
  if (!in.good()) throw std::runtime_error(std::string("gro_reader: cannot open file: ") + std::string(filename));

  return read_gro_stream_to_rwmol(in);
}

} // namespace lahuta::md
