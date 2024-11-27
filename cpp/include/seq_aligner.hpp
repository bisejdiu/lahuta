#ifndef LAHUTA_SEQ_ALIGNER_HPP
#define LAHUTA_SEQ_ALIGNER_HPP

#include "align.hpp"
#include "extract.hpp"
#include "fseek/ops.hpp"
#include <StructureSmithWaterman.h>
#include <SubstitutionMatrix.h>
#include <TMaligner.h>

namespace lahuta {

class AlignmentScores;

// NOTE: uniqe_ptr because the underlying objects do not provide move semantics
struct MatrixContainer {
  std::unique_ptr<SubstitutionMatrix> subMat3Di;
  std::unique_ptr<SubstitutionMatrix> subMatAA;
  std::unique_ptr<StructureSmithWaterman> sw;
  std::unique_ptr<StructureSmithWaterman> rsw;
  std::vector<int8_t> tinySubMatAA;
  std::vector<int8_t> tinySubMat3Di;
};

struct AlignmentResult {
  std::vector<Matcher::result_t> alignmentResult;
  std::string resultBuffer;
  TMaligner::TMscoreResult tmres;
  Matcher::result_t res;
  std::unique_ptr<AlignmentScores> scores;
  bool success{false};
};

class AlignmentScores {
  struct TMscoreFunctor {
    SeqData &target;
    const Matcher::result_t &res;
    const std::string backtrace;

    TMscoreFunctor(SeqData &target_, const Matcher::result_t &res_)
        : target(target_), res(res_), backtrace(Matcher::uncompressAlignment(res.backtrace)) {}

    TMaligner::TMscoreResult operator()(TMaligner *aligner, int len) const {
      return aligner->computeTMscore(
          target.x(),
          target.y(),
          target.z(),
          res.dbLen,
          res.qStartPos,
          res.dbStartPos,
          backtrace.c_str(),
          len);
    }

    LDDTCalculator::LDDTScoreResult operator()(LDDTCalculator *lddt_calculator, int len) const {
      return lddt_calculator->computeLDDTScore(
          res.dbLen,
          res.qStartPos,
          res.dbStartPos,
          backtrace,
          target.x(),
          target.y(),
          target.z());
    }
  };

public:
  AlignmentScores(SeqData &query, SeqData &target, const Matcher::result_t &ar)
      : calculator(target, ar), //
        tm_calculator(std::max(query.size(), target.size()) + 1, false, true, false),
        lddt_calculator(query.size() + 1, target.size() + 1),
        alignment_len(std::min(ar.qEndPos - ar.qStartPos, ar.dbEndPos - ar.dbStartPos)), query_len(ar.qLen),
        target_len(ar.dbLen) {

    // FIX: move out of constructor and make optional
    tm_calculator.initQuery(query.x(), query.y(), query.z(), nullptr, ar.qLen);
    lddt_calculator.initQuery(query.size(), query.x(), query.y(), query.z());
  }

  auto get_alignment_score() { return calculator(&tm_calculator, alignment_len); }
  auto get_query_score() { return calculator(&tm_calculator, query_len); }
  auto get_target_score() { return calculator(&tm_calculator, target_len); }
  auto get_lddt_score() { return calculator(&lddt_calculator, alignment_len); }

private:
  TMscoreFunctor calculator;
  TMaligner tm_calculator;
  LDDTCalculator lddt_calculator;
  const int alignment_len;
  const int query_len;
  const int target_len;
};

class SeqAligner {
public:
  SeqAligner(const SeqAligner &) = delete;
  SeqAligner &operator=(const SeqAligner &) = delete;
  SeqAligner(SeqAligner &&) = delete;
  SeqAligner &operator=(SeqAligner &&) = delete;
  ~SeqAligner() = default;

  void initialize(int max_q_len, int max_t_len) {
    auto max_len = std::max(max_q_len, max_t_len);
    tmaligner = std::make_unique<TMaligner>(max_len + 1, false, true, ops.exactTMscore);
    lddt_calculator = std::make_unique<LDDTCalculator>(max_q_len + 1, max_t_len + 1);
    is_initialized = true;
  }

  AlignmentResult align(SeqData &Q, SeqData &T, FoldSeekOps &ops) {

    if (!is_initialized) throw std::runtime_error("SeqAligner is not initialized");

    std::vector<Matcher::result_t> alignment_result;

    // FIX: not really needed
    std::string backtrace;
    char buffer[1024 + 32768];
    std::string result_buffer;

    // FIX: not really needed
    TMaligner::TMscoreResult tmres;
    LDDTCalculator::LDDTScoreResult lddtres;

    // FIX: can be computed only when query or target changes
    std::unique_ptr<Sequence> q3 = Q.map_3di(*M.subMat3Di, ops);
    std::unique_ptr<Sequence> qa = Q.map_aa(*M.subMatAA, ops);
    std::unique_ptr<Sequence> t3 = T.map_3di(*M.subMat3Di, ops);
    std::unique_ptr<Sequence> ta = T.map_aa(*M.subMatAA, ops);

    // FIX: this is a setup parameter. It should not be necessary to check every time
    if (need_tmaligner) {
      tmaligner->initQuery(Q.x(), Q.y(), Q.z(), NULL, Q.size());
    }
    if (need_lddt) {
      lddt_calculator->initQuery(Q.size(), Q.x(), Q.y(), Q.z());
    }

    std::pair<double, double> mu_lambda = evaluer.predictMuLambda(q3->numSequence, q3->L);
    std::cout << "mu: " << mu_lambda.first << " lambda: " << mu_lambda.second << std::endl;

    init_sw(*M.sw, *q3, *qa);
    q3->reverse();
    qa->reverse();
    init_sw(*M.rsw, *q3, *qa);

    if (!Util::canBeCovered(ops.covThr, ops.covMode, Q.Seq3Di.size(), T.SeqAA.size())) return {};

    Matcher::result_t res;

    bool success = align_structure(*ta, *t3, q3->L, res, mu_lambda, backtrace);
    if (!success) return {};

    // NOTE: this part of the code must be run in parallel
    if (!alignment_criteria_valid(res)) return {};

    if (need_tmaligner || need_lddt) {

      if (need_tmaligner) {
        tmres = compute_tm(T, res);

        if (tmres.tmscore < ops.tmScoreThr) {
          std::cout << "TMscore is lower than threshold" << std::endl;
          return {};
        }
      }
      if (need_lddt) {
        lddtres = compute_lddt(res, T);

        if (lddtres.avgLddtScore < ops.lddtThr) {
          std::cout << "LDDT score is lower than threshold" << std::endl;
          return {};
        }
        res.dbcov = lddtres.avgLddtScore;
      }
      if (ops.sortByStructureBits && need_tmaligner && need_lddt) {
        res.score = res.score * sqrt(lddtres.avgLddtScore * tmres.tmscore);
      }
    }

    alignment_result.emplace_back(res);
    for (int alt_ali = ops.altAlignment; alt_ali > 0; --alt_ali) {
      Matcher::result_t alt_res;

      bool success = alt_align_structure(*ta, *t3, q3->L, res, alt_res, mu_lambda, backtrace);
      if (!success) break;

      alignment_result.push_back(alt_res);
      res = alt_res;
    }

    if (alignment_result.size() > 1) {
      if (ops.sortByStructureBits) {
        std::sort(alignment_result.begin(), alignment_result.end(), foldseek::compareHitsByStructureBits);
      } else {
        std::sort(alignment_result.begin(), alignment_result.end(), Matcher::compareHits);
      }
    }

    for (size_t result = 0; result < alignment_result.size(); result++) {
      auto resutl = alignment_result[result];
      size_t len = Matcher::resultToBuffer(buffer, alignment_result[result], ops.addBacktrace);
      result_buffer.append(buffer, len);
    }

    std::cout << "Result: " << result_buffer << std::endl;
    return {alignment_result, result_buffer, tmres, res, std::make_unique<AlignmentScores>(Q, T, res), true};
  }

private:
  /// target_length is the combined length of all sequences in the target collection
  SeqAligner(FoldSeekOps &ops_, MatrixContainer &c, int target_length)
      : ops(ops_), M(std::move(c)), evaluer(target_length, M.subMat3Di.get()) {}

  FoldSeekOps &ops;
  MatrixContainer M;

  std::unique_ptr<TMaligner> tmaligner;
  std::unique_ptr<LDDTCalculator> lddt_calculator;
  EvalueNeuralNet evaluer;

  bool need_tmaligner = ops.sortByStructureBits;
  bool need_lddt = ops.sortByStructureBits ? true : ops.lddtThr > 0;
  bool is_initialized = false;

  friend class SeqAlignerBuilder;

  void init_sw(StructureSmithWaterman &sw, Sequence &s3, Sequence &sa) {
    sw.ssw_init(&sa, &s3, M.tinySubMatAA.data(), M.tinySubMat3Di.data(), M.subMatAA.get());
  }

  bool alignment_criteria_valid(Matcher::result_t &res) {
    bool is_identity = false;
    if (alignmentCheckCriteria(
            res,
            is_identity,
            ops.evalThr,
            ops.seqIdThr,
            ops.alnLenThr,
            ops.covMode,
            ops.covThr)) {
      return true;
    }
    return false;
  }

  bool align_structure(
      Sequence &sa, Sequence &s3, int q_len, Matcher::result_t &res, std::pair<double, double> &mu_lambda,
      std::string &backtrace) {
    auto result =
        foldseek::alignStructure(*M.sw, *M.rsw, sa, s3, q_len, sa.L, evaluer, mu_lambda, res, backtrace, ops);
    return result != -1;
  }

  bool alt_align_structure(
      Sequence &sa, Sequence &s3, int q_len, Matcher::result_t &res, Matcher::result_t &alt_res,
      std::pair<double, double> &mu_lambda, std::string &backtrace) {
    auto result = foldseek::computeAlternativeAlignment(
        *M.sw,
        *M.rsw,
        sa,
        s3,
        q_len,
        sa.L,
        evaluer,
        mu_lambda,
        res,
        alt_res,
        backtrace,
        ops);
    return result != -1;
  }

  LDDTCalculator::LDDTScoreResult compute_lddt(Matcher::result_t &res, SeqData &s) {
    return lddt_calculator
        ->computeLDDTScore(res.dbLen, res.qStartPos, res.dbStartPos, res.backtrace, s.x(), s.y(), s.z());
  }

  TMaligner::TMscoreResult compute_tm(SeqData &s, Matcher::result_t &res) {
    unsigned int norm_mode = TMaligner::normalization(
        ops.tmScoreThrMode,
        std::min(res.qEndPos - res.qStartPos, res.dbEndPos - res.dbStartPos),
        res.qLen,
        res.dbLen);

    return tmaligner->computeTMscore(
        s.x(),
        s.y(),
        s.z(),
        res.dbLen,
        res.qStartPos,
        res.dbStartPos,
        res.backtrace,
        norm_mode);
  }
};

} // namespace lahuta

#endif // SEQ_ALIGNER_HPP
