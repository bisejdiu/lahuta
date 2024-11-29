#ifndef LAHUTA_SEQ_ALIGNER_HPP
#define LAHUTA_SEQ_ALIGNER_HPP

#include "CalcProbTP.h"
#include "align.hpp"
#include "fseek/ops.hpp"
#include "matcher.hpp"
#include "seq.hpp"
#include <StructureSmithWaterman.h>
#include <SubstitutionMatrix.h>
#include <TMaligner.h>
#include <array>

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

struct Scores {
  double rmsd;
  double tmscore;
  double avgLddtScore;
  double prob;
};


struct AlignmentResult {
  // FIX: why do we need both Matcher::result_t and a vector of Matcher::result_t?
  std::vector<Matcher::result_t> ar;
  Scores scores;
  /*std::string resultBuffer;*/
  TMaligner::TMscoreResult tmres;
  /*Matcher::result_t res;*/
  /*std::shared_ptr<AlignmentScores> scores;*/
  bool success{false};
  std::unique_ptr<SeqData> query;
  std::unique_ptr<SeqData> target;

  std::string query_alignment() const {
    return alignment_from_cigar(query->SeqAA.c_str(), ar[0].qStartPos, SeqType::Query);
  }

  std::string target_alignment() const {
    return alignment_from_cigar(target->SeqAA.c_str(), ar[0].dbStartPos, SeqType::Target);
  }

private:
  enum class SeqType { Query = 0, Target = 1 };

  // constexpr not possible with std::function<bool(char)>
  constexpr static std::array<bool (*)(char), 2> behaviors = {{
      [](char symbol) { return symbol == 'M' || symbol == 'I'; }, // SeqType::Query
      [](char symbol) { return symbol == 'M' || symbol == 'D'; }, // SeqType::Target
  }};

  std::string alignment_from_cigar(const char *seq, unsigned int offset, SeqType type) const {
    const auto &advance_seq_pos = behaviors[static_cast<size_t>(type)];
    std::string out{};
    unsigned int seq_pos{0};
    /*std::string backtrace = Matcher::uncompressAlignment(res.backtrace);*/
    for (const auto &symbol : ar[0].backtrace) {
      if (advance_seq_pos(symbol)) {
        out.push_back(seq[offset + seq_pos]);
        seq_pos++;
      } else {
        out.push_back('-');
      }
    }
    return out;
  }
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

    Scores scores;
    TMaligner::TMscoreResult tmres;
    LDDTCalculator::LDDTScoreResult lddtres;
    std::vector<Matcher::result_t> result;

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

    init_sw(*M.sw, *q3, *qa);
    q3->reverse();
    qa->reverse();
    init_sw(*M.rsw, *q3, *qa);

    if (!Util::canBeCovered(ops.covThr, ops.covMode, Q.Seq3Di.size(), T.SeqAA.size())) return {};

    Matcher::result_t res;
    std::string backtrace;
    bool success = align_structure(*ta, *t3, q3->L, res, mu_lambda, backtrace);
    if (!success) return {};

    // NOTE: this part of the code must be run in parallel
    if (!alignment_criteria_valid(res)) return {};

    if (need_tmaligner || need_lddt) {

      if (need_tmaligner) {
        tmres = compute_tm(T, res);
        std::cout << "TMscore: " << tmres.tmscore << std::endl;
        std::cout << "RMSD: " << tmres.rmsd << std::endl;

        if (tmres.tmscore < ops.tmScoreThr) {
          std::cout << "TMscore is lower than threshold" << std::endl;
          return {};
        }

        scores.tmscore = tmres.tmscore;
        scores.rmsd = tmres.rmsd;
      }
      if (need_lddt) {
        lddtres = compute_lddt(res, T);
        std::cout << "lddtres: " << lddtres.avgLddtScore << std::endl;

        if (lddtres.avgLddtScore < ops.lddtThr) {
          std::cout << "LDDT score is lower than threshold" << std::endl;
          return {};
        }
        res.dbcov = lddtres.avgLddtScore;
        scores.avgLddtScore = lddtres.avgLddtScore;
      }
      if (ops.sortByStructureBits && need_tmaligner && need_lddt) {
        res.score = res.score * sqrt(lddtres.avgLddtScore * tmres.tmscore);
      }
    }

    result.push_back(res);
    for (int alt_ali = ops.altAlignment; alt_ali > 0; --alt_ali) {
      Matcher::result_t alt_res;

      bool success = alt_align_structure(*ta, *t3, q3->L, res, alt_res, mu_lambda, backtrace);
      if (!success) break;

      if (!alignment_criteria_valid(alt_res)) break;
      if (alt_res.backtrace == result.back().backtrace) break; // in-lieu of convergence check
      result.push_back(alt_res);
    }

    auto sorter = ops.sortByStructureBits ? foldseek::compareHitsByStructureBits : Matcher::compareHits;
    std::sort(result.begin(), result.end(), sorter);

    scores.prob = CalcProbTP::calculate(result.front().score);
    AlignmentResult ro = {result, scores, tmres, true};
    ro.query = std::make_unique<SeqData>(Q);
    ro.target = std::make_unique<SeqData>(T);
    return ro;
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
