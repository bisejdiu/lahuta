#include "CalcProbTP.h"
#include "align.hpp"
#include "file_system.hpp"
#include "fseek/ops.hpp"
#include "lahuta.hpp"

#include "extract.hpp"
#include "gemmi/resinfo.hpp"
#include "prefilter.hpp"
#include "seq_aligner_builder.hpp"

using namespace lahuta;

std::string format_seq(const std::string &seq) {
  std::string formatted_seq;
  formatted_seq += seq.substr(0, 10);
  formatted_seq += "...";
  formatted_seq += seq.substr(seq.size() - 10, 10);
  return formatted_seq;
}

std::string alignment_from_cigar(const char *seq, const unsigned int offset, const std::string &bt) {

  std::string out{};
  unsigned int seq_pos{0};
  for (const auto &symbol : bt) {
    switch (symbol) {
      case 'M':
      case 'D':
        out.push_back(seq[offset + seq_pos]);
        seq_pos++;
        break;
      case 'I':
        out.push_back('-');
        break;
    }
  }
  return out;
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

void log_sequence(const SeqData &sd, const std::string &prefix = "") {
  std::cout << prefix << " Sequence: " << sd.file_name << " " << sd.chain_name << " " << sd.Seq3Di.size()
            << " residues: " << format_seq(sd.SeqAA) << std::endl;
}

int main(int argc, char const *argv[]) {

  FoldSeekOps ops;
  PrefilterOptions pf_ops;

  std::string source_file_name = "2rh1.cif";
  std::string target_file_name = "3sn6.cif";

  std::vector<std::string> fileNames = {source_file_name, target_file_name};

  SeqCollection queries = extract_all(ops, {fileNames[0]});
  SeqCollection targets = extract_all(ops, {fileNames[1]});

  std::cout << "Queries: " << queries.size() << std::endl;
  std::cout << "Targets: " << targets.size() << std::endl;

  SeqFilter seq_filter(queries, targets, pf_ops); // TODO: make pf_opts optional?
  seq_filter.build_index();

  SeqAlignerBuilder builder(ops);

  int query_count = 0;
  std::unique_ptr<SeqAligner> aligner = builder.build(queries, targets);

  for (auto &query : queries) {
    log_sequence(query, "Q:");
    Hits hits = seq_filter.filter(query);

    for (const auto &hit_id : hits) {
      SeqData target = targets[hit_id];
      log_sequence(target, "T:");
      auto AR = aligner->align(query, target, ops);
      if (!AR.success) continue;

      std::string filename = "3sn6_aligned" + std::to_string(query_count++) + ".pdb";
      // std::string filename = "aligned_" + query.file_name + "_" + target.file_name + ".pdb";

      Transform3D transform(AR.tmres.t, AR.tmres.u);
      write_aligned_pdb(filename, target, &transform);
      auto out = alignment_from_cigar(
          target.SeqAA.c_str(),
          AR.res.dbStartPos,
          Matcher::uncompressAlignment(AR.res.backtrace));
      std::cout << "Alignment: " << out << std::endl;

      auto lddt_result = AR.scores->get_lddt_score();
      auto alig_result = AR.scores->get_alignment_score();
      auto prob = CalcProbTP::calculate(AR.res.score);

      std::cout << "----------------------------------------" << std::endl;
      std::cout << "RMSD: " << alig_result.rmsd << std::endl;
      std::cout << "LDDT: " << lddt_result.avgLddtScore << std::endl;
      std::cout << "TMscore: " << alig_result.tmscore << std::endl;
      std::cout << "QTMscore: " << AR.scores->get_query_score().tmscore << std::endl;
      std::cout << "TTMscore: " << AR.scores->get_target_score().tmscore << std::endl;
      std::cout << "Prob: " << prob << std::endl;
      std::cout << "----------------------------------------" << std::endl;
    }
  }
  return 0;
};
