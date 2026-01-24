#ifndef LAHUTA_FSEEK_UTILS_HPP
#define LAHUTA_FSEEK_UTILS_HPP

#include <array>
#include <fstream>
#include <iostream>
#include <string>

#include <gemmi/resinfo.hpp>
#include <rdkit/GraphMol/Conformer.h>
#include <rdkit/GraphMol/RWMol.h>

#include "fseek/seq_aligner.hpp"
#include "lahuta.hpp"
#include "matcher.hpp"

namespace lahuta {

inline std::string format_seq(const std::string &seq) {
  std::string formatted_seq;
  int seq_size_to_print = 50;
  formatted_seq += seq.substr(0, seq_size_to_print);
  formatted_seq += "...";
  formatted_seq += seq.substr(seq.size() - seq_size_to_print, seq_size_to_print);
  return formatted_seq;
}

class Transform3D {
public:
  static constexpr std::array<float, 3> IdentityTranslation = {0.0f, 0.0f, 0.0f};
  static constexpr std::array<std::array<float, 3>, 3> IdentityRotation = {
      {{1.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f}, {0.0f, 0.0f, 1.0f}}};

  constexpr Transform3D() noexcept : translation_(IdentityTranslation), rotation_(IdentityRotation) {}

  constexpr Transform3D(
      const std::array<float, 3> &translation, const std::array<std::array<float, 3>, 3> &rotation) noexcept
      : translation_(translation), rotation_(rotation) {}

  Transform3D(const float translation[3], const float rotation[3][3]) noexcept
      : translation_({translation[0], translation[1], translation[2]}),
        rotation_(
            {{{rotation[0][0], rotation[0][1], rotation[0][2]},
              {rotation[1][0], rotation[1][1], rotation[1][2]},
              {rotation[2][0], rotation[2][1], rotation[2][2]}}}) {}

  constexpr std::array<float, 3> operator()(float x, float y, float z) const noexcept {
    return {
        translation_[0] + x * rotation_[0][0] + y * rotation_[0][1] + z * rotation_[0][2],
        translation_[1] + x * rotation_[1][0] + y * rotation_[1][1] + z * rotation_[1][2],
        translation_[2] + x * rotation_[2][0] + y * rotation_[2][1] + z * rotation_[2][2]};
  }

private:
  std::array<float, 3> translation_;
  std::array<std::array<float, 3>, 3> rotation_;
};

// FIX: SeqData does not store residue IDs
// FIX: Save the entire structure, not just CA
inline void
write_aligned_pdb(const std::string &filename, const SeqData &sd, const Transform3D *transform = nullptr) {
  std::ofstream out_file(filename);
  if (!out_file.is_open()) {
    std::cerr << "Failed to open file: " << filename << '\n';
    return;
  }

  std::ostringstream result;
  result << "MODEL\n";

  if (!transform) transform = new Transform3D();
  for (size_t tpos = 0; tpos < sd.SeqAA.size(); ++tpos) {
    float tx = sd.CaData[tpos];
    float ty = sd.CaData[sd.SeqAA.size() + tpos];
    float tz = sd.CaData[2 * sd.SeqAA.size() + tpos];

    auto [x, y, z] = (*transform)(tx, ty, tz);

    // clang-format off
    result << std::setw(6) << std::left << "ATOM"
           << std::setw(5) << std::right << tpos + 1 << ' '
           << std::setw(4) << std::right << "CA"
           << ' ' << std::setw(3) << std::left << gemmi::expand_protein_one_letter(sd.SeqAA[tpos]) << ' '
           << 'A' << std::setw(4) << std::right << tpos + 1 << "    "
           << std::fixed << std::setprecision(3)
           << std::setw(8) << x
           << std::setw(8) << y
           << std::setw(8) << z
           << std::setw(6) << std::fixed << std::setprecision(2) << 1.0
           << std::setw(6) << std::fixed << std::setprecision(2) << 0.0
           << '\n';
    // clang-format on
  }

  result << "ENDMDL\n";
  out_file << result.str();
}

inline void write_aligned_pdb_w(
    const std::string &filename, const SeqData &sd, const Matcher::result_t &ar,
    const Transform3D *transform = nullptr) {
  std::ofstream out_file(filename);
  if (!out_file.is_open()) {
    std::cerr << "Failed to open file: " << filename << '\n';
    return;
  }

  auto luni = Luni::create(sd.st);
  const RDKit::RWMol &mol = luni.get_molecule();
  const RDKit::Conformer &conf = luni.get_conformer();

  Transform3D local_transform = transform ? *transform : Transform3D();

  std::ostringstream result;
  result << "MODEL\n";

  std::cout << "DEBUG: " << ar.dbStartPos << " " << ar.dbEndPos << std::endl;
  for (unsigned int atom_id = 0; atom_id < mol.getNumAtoms(); ++atom_id) {
    const auto atom = mol.getAtomWithIdx(atom_id);
    auto *ri = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
    if (!ri) continue;

    std::string atom_name = ri->getName();
    std::string residue_name = ri->getResidueName();
    auto chain_id = ri->getChainId();
    int residue_num = ri->getResidueNumber();

    if (chain_id != sd.chain_name) continue;
    /*if (residue_num < ar.dbStartPos || residue_num > ar.dbEndPos) continue;*/

    const RDGeom::Point3D &pos = conf.getAtomPos(atom_id);

    auto [x, y, z] =
        local_transform(static_cast<float>(pos.x), static_cast<float>(pos.y), static_cast<float>(pos.z));

    // clang-format off
    result << std::setw(6) << std::left << "ATOM"
           << std::setw(5) << std::right << atom_id + 1 << ' '
           << std::setw(4) << std::right << atom_name
           << ' ' << std::setw(3) << std::left << residue_name << ' '
           << 'A' << std::setw(4) << std::right << residue_num << "    "
           << std::fixed << std::setprecision(3)
           << std::setw(8) << x
           << std::setw(8) << y
           << std::setw(8) << z
           << std::setw(6) << std::fixed << std::setprecision(2) << 1.0
           << std::setw(6) << std::fixed << std::setprecision(2) << 0.0
           << '\n';
    // clang-format on
  }

  result << "ENDMDL\n";
  out_file << result.str();
}

inline void log_sequence(const SeqData &sd, const std::string &prefix = "") {
  std::cout << prefix << " Sequence: " << sd.file_name << " " << sd.chain_name //
            << " " << sd.Seq3Di.size() << " residues: " << format_seq(sd.SeqAA) << std::endl;
}

inline void save_file(const std::string &filename, SeqData &sd, const AlignmentResult &ar) {
  Transform3D transform(ar.tmres.t, ar.tmres.u);
  write_aligned_pdb(filename, sd, &transform);
};

inline void save_file_w(const std::string &filename, SeqData &sd, const AlignmentResult &ar) {
  Transform3D transform(ar.tmres.t, ar.tmres.u);
  write_aligned_pdb_w(filename, sd, ar.ar.front(), &transform);
};

} // namespace lahuta

#endif // LAHUTA_FSEEK_UTILS_HPP
