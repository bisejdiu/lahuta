#ifndef LAHUTA_OPS_HPP
#define LAHUTA_OPS_HPP

#include "blosum62.out.h"
#include "mat3di.out.h"
#include "TMaligner.h"

enum class AlignType { AA, _3Di, AA_3Di };

struct SubMatrixData {
  const char *name;
  const unsigned char *data;
  const unsigned int len;
  const unsigned int nameLen;
};

constexpr SubMatrixData subMatrices[] = {
    {"3di.out", mat3di_out, mat3di_out_len, 7}, //
    {"blosum62.out", blosum62_out, blosum62_out_len, 12}};

class CompileTimeMatrix {
public:
  constexpr CompileTimeMatrix(const SubMatrixData &matrix)
      : name(matrix.name), data(matrix.data), size(matrix.len), name_size(matrix.nameLen) {}

  const char *get_data_as_char() const { return reinterpret_cast<const char *>(data); }
  const unsigned char *get_data_as_uchar() const { return data; }

public:
  const char *name;
  const unsigned char *data;
  const unsigned int size;
  const unsigned int name_size;
};

namespace lahuta {

struct FoldSeekOps {
  AlignType alignType = AlignType::AA_3Di;

  // Accept alignments with a TMscore > threshold [0.0,1.0]
  float tmScoreThr = 0.0;

  // TMscore threshold mode: 0 - alignment, 1 - query, 2 - target length, 3 - min(query, target)
  TMScoreThrMode tmScoreThrMode = TMScoreThrMode::alignment;

  // Turn on fast exact TMscore (slower), default is approximate
  bool exactTMscore = false;

  // Accept alignments with a lddt > threshold [0.0,1.0]
  float lddtThr = 0.0;

  // Sort by structure bit score: bits * sqrt(alnlddt * alntmscore)
  bool sortByStructureBits = true;

  // Alignment computation mode:
  // 0 - automatic, 1 - only score and end_pos, 2 - also start_pos and cov,
  // 3 - also seq.id, 4 - only ungapped alignment
  int alignmentMode = 3;

  // Alignment output mode:
  // 0 - automatic, 1 - only score and end_pos, 2 - also start_pos and cov,
  // 3 - also seq.id, 4 - only ungapped alignment, 5 - score only (output) cluster format
  int alignmentOutputMode = 0;

  // Allow wrapped scoring by doubling the (nucleotide) query sequence for scoring
  bool wrappedScoring = false;

  // Maximum sequence length
  int maxSeqLen = 65535;

  // Correct for locally biased amino acid composition (range 0-1)
  bool compBiasCorrection = true;

  // Scale for composition bias correction (range 0-1)
  float compBiasCorrectionScale = 0.5;

  // Score bias for SW alignment (in bits)
  float scoreBias = 0.0;

  // Compute more conservative, shorter alignments (scores and E-values unchanged)
  bool realign = false;

  // Weight of backtrace correlation score added to the alignment score
  float correlationScoreWeight = 0.0;

  // Add backtrace string, used for converting to alignments with mmseqs convertalis
  int addBacktrace = 1;

  // Coverage threshold: List matches above this fraction of aligned residues
  float covThr = 0.0;

  // Coverage mode:
  // 0 - query and target, 1 - target, 2 - query,
  // 3 - target length >= x% of query length,
  // 4 - query length >= x% of target length,
  // 5 - short seq. >= x% of other seq. length
  int covMode = 0;

  // List matches below this E-value (range 0.0-inf)
  int evalThr = 10;

  // Sequence identity threshold for clustering (range 0.0-1.0)
  float seqIdThr = 0.0;

  // Sequence identity mode:
  // 0 - alignment length, 1 - shorter, 2 - longer sequence
  bool seqIdMode = 0;

  // Minimum alignment length (range 0-INT_MAX)
  float alnLenThr = 0.0;

  // Chain name mode: 0 - auto, 1 - always add
  int chainNameMode = 0;

  // Mask b-factor threshold (mask residues for seeding if b-factor < threshold [0,100])
  float maskBfactorThreshold = 0.0;

  // number of alternative alignments
  int altAlignment{0};

  // Format of input structures:
  // 0: Auto-detect by extension
  // 1: PDB
  // 2: mmCIF
  int inputFormat = 0;

  // Gap opening and extension penalties
  int gapOpen = 10;
  int gapExtend = 1;
};

constexpr const CompileTimeMatrix SubMatrix3Di(subMatrices[0]);
constexpr const CompileTimeMatrix SubMatrixBlosum62(subMatrices[1]);

namespace {
constexpr bool compare_strings(const char *str1, const char *str2) {
  while (*str1 && *str2) {
    if (*str1 != *str2) return false;
    ++str1;
    ++str2;
  }
  return *str1 == *str2;
}
static_assert(compare_strings(SubMatrix3Di.name, "3di.out"), "Matrix name is not correct");
} // namespace

} // namespace lahuta

#endif
