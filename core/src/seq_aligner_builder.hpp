#ifndef LAHUTA_SEQ_ALIGNER_BUILDER_HPP
#define LAHUTA_SEQ_ALIGNER_BUILDER_HPP

#include "fseek/ops.hpp"
#include "seq_aligner.hpp"
#include <memory>

namespace lahuta {

class SeqAlignerBuilder {
public:
  explicit SeqAlignerBuilder(FoldSeekOps &ops) : ops(ops) {}

  std::unique_ptr<SeqAligner> build(const SeqCollection &Q, const SeqCollection &T) {
    auto matrix3Di = get_3Di_matrix();
    auto matrixAA = get_aa_matrix();

    MatrixContainer M{
        std::move(matrix3Di),
        std::move(matrixAA),
        init_sw(*M.subMatAA, *M.subMat3Di),
        init_sw(*M.subMatAA, *M.subMat3Di),
        init_tiny_matrix(*M.subMatAA),
        init_tiny_matrix(*M.subMat3Di)};

    auto seq_aligner = std::unique_ptr<SeqAligner>(new SeqAligner(ops, M, T.total_length));
    /*auto seq_aligner = std::unique_ptr<SeqAligner>(ops, M, T.total_length);*/
    seq_aligner->initialize(Q.max_length, T.max_length);

    return seq_aligner;
  }

private:
  std::unique_ptr<Sequence> build_sequence(SubstitutionMatrix &matrix) {
    return std::make_unique<Sequence>(ops.maxSeqLen, &matrix, ops.compBiasCorrection);
  }

  std::unique_ptr<SubstitutionMatrix> get_3Di_matrix() {
    return std::make_unique<SubstitutionMatrix>(SubMatrix3Di, 2.1, ops.scoreBias);
  }

  std::unique_ptr<SubstitutionMatrix> get_aa_matrix() {
    float aa_factor = (ops.alignType == AlignType::AA_3Di) ? 1.4f : 0.0f;
    return std::make_unique<SubstitutionMatrix>(SubMatrixBlosum62, aa_factor, ops.scoreBias);
  }

  std::vector<int8_t> init_tiny_matrix(const SubstitutionMatrix &matrix) {
    std::vector<int8_t> tinyMatrix(matrix.alphabetSize * 32);

    for (int i = 0; i < matrix.alphabetSize; ++i) {
      for (int j = 0; j < matrix.alphabetSize; ++j) {
        tinyMatrix[i * matrix.alphabetSize + j] = matrix.subMatrix[i][j];
      }
    }

    return std::move(tinyMatrix);
  }

  std::unique_ptr<StructureSmithWaterman>
  init_sw(SubstitutionMatrix &matrixAA, SubstitutionMatrix &matrix3Di) {
    return std::make_unique<StructureSmithWaterman>(
        ops.maxSeqLen,
        matrix3Di.alphabetSize,
        ops.compBiasCorrection,
        ops.compBiasCorrectionScale,
        &matrixAA,
        &matrix3Di);
  }

  FoldSeekOps &ops;
};

} // namespace lahuta

#endif // LAHUTA_SEQ_ALIGNER_BUILDER_HPP
