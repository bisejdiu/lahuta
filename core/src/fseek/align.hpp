#ifndef LAHUTA_ALIGN_HPP
#define LAHUTA_ALIGN_HPP

#include "Debug.h"
#include "EvalueNeuralNet.h"
#include "LDDT.h"
#include "StructureSmithWaterman.h"
#include "fseek/ops.hpp"
#include "matcher.hpp"

namespace lahuta {

namespace foldseek {

static bool compareHitsByStructureBits(const Matcher::result_t &first, const Matcher::result_t &second) {
  if (first.score != second.score) {
    return first.score > second.score;
  }
  if (first.dbLen != second.dbLen) {
    return first.dbLen < second.dbLen;
  }
  return first.dbKey < second.dbKey;
}

inline int alignStructure(
    StructureSmithWaterman &ssm, StructureSmithWaterman &r_ssm, Sequence &tSeqAA, Sequence &tSeq3Di,
    unsigned int querySeqLen, unsigned int targetSeqLen, EvalueNeuralNet &evaluer,
    std::pair<double, double> muLambda, Matcher::result_t &res, std::string &backtrace, FoldSeekOps &ops) {

  using SSW = StructureSmithWaterman;

  float seqId = 0.0;
  backtrace.clear();
  // align only score and end pos
  auto align = ssm.alignScoreEndPos<SSW::PROFILE>(
      tSeqAA.numSequence,
      tSeq3Di.numSequence,
      targetSeqLen,
      ops.gapOpen,
      ops.gapExtend,
      querySeqLen / 2);

  if (!Util::hasCoverage(ops.covThr, ops.covMode, align.qCov, align.tCov)) return -1;

  // we can already stop if this e-value isn't good enough, it wont be any better in the next step
  align.evalue = evaluer.computeEvalueCorr(align.score1, muLambda.first, muLambda.second);
  if (align.evalue > ops.evalThr) return -1;

  StructureSmithWaterman::s_align rev_align;
  if (ssm.isProfileSearch()) {
    rev_align.score1 = 0;
  } else {
    rev_align = r_ssm.alignScoreEndPos<SSW::PROFILE>(
        tSeqAA.numSequence,
        tSeq3Di.numSequence,
        targetSeqLen,
        ops.gapOpen,
        ops.gapExtend,
        querySeqLen / 2);
  }

  int32_t score = static_cast<int32_t>(align.score1) - static_cast<int32_t>(rev_align.score1);
  align.evalue = evaluer.computeEvalueCorr(score, muLambda.first, muLambda.second);
  if (align.evalue > ops.evalThr) return -1;

  bool blockAlignFailed = false;
  if (ssm.isProfileSearch() == false) {
    StructureSmithWaterman::s_align alignTmp = ssm.alignStartPosBacktraceBlock(
        tSeqAA.numSequence,
        tSeq3Di.numSequence,
        targetSeqLen,
        ops.gapOpen,
        ops.gapExtend,
        backtrace,
        align);

    if (align.score1 == UINT32_MAX) {
      Debug(Debug::WARNING) << "block-align failed, falling back to normal alignment\n";
      blockAlignFailed = true;
    } else {
      align = alignTmp;
    }
  }

  if (blockAlignFailed || ssm.isProfileSearch()) {
    align = ssm.alignStartPosBacktrace<StructureSmithWaterman::PROFILE>(
        tSeqAA.numSequence,
        tSeq3Di.numSequence,
        targetSeqLen,
        ops.gapOpen,
        ops.gapExtend,
        ops.alignmentMode,
        backtrace,
        align,
        ops.covMode,
        ops.covThr,
        querySeqLen / 2);
  }

  unsigned int alnLength =
      Matcher::computeAlnLength(align.qStartPos1, align.qEndPos1, align.dbStartPos1, align.dbEndPos1);
  if (backtrace.size() > 0) {
    alnLength = backtrace.size();
    seqId = Util::computeSeqId(ops.seqIdMode, align.identicalAACnt, querySeqLen, targetSeqLen, alnLength);
  }
  align.score1 = score;
  res = Matcher::result_t(
      tSeqAA.getDbKey(),
      align.score1,
      align.qCov,
      align.tCov,
      seqId,
      align.evalue,
      alnLength,
      align.qStartPos1,
      align.qEndPos1,
      querySeqLen,
      align.dbStartPos1,
      align.dbEndPos1,
      targetSeqLen,
      backtrace);
  return 0;
}

inline int computeAlternativeAlignment(
    StructureSmithWaterman &structureSmithWaterman, StructureSmithWaterman &reverseStructureSmithWaterman,
    Sequence &tSeqAA, Sequence &tSeq3Di, unsigned int querySeqLen, unsigned int targetSeqLen,
    EvalueNeuralNet &evaluer, std::pair<double, double> muLambda, Matcher::result_t &result,
    Matcher::result_t &altRes, std::string &backtrace, FoldSeekOps &ops) {

  const unsigned char xAAIndex = tSeqAA.subMat->aa2num[static_cast<int>('X')];
  const unsigned char x3DiIndex = tSeq3Di.subMat->aa2num[static_cast<int>('X')];
  for (int pos = result.dbStartPos; pos < result.dbEndPos; ++pos) {
    tSeqAA.numSequence[pos] = xAAIndex;
    tSeq3Di.numSequence[pos] = x3DiIndex;
  }
  if (alignStructure(
          structureSmithWaterman,
          reverseStructureSmithWaterman,
          tSeqAA,
          tSeq3Di,
          querySeqLen,
          targetSeqLen,
          evaluer,
          muLambda,
          altRes,
          backtrace,
          ops)
      == -1) {
    return -1;
  }
  if (alignmentCheckCriteria(
          altRes,
          false,
          ops.evalThr,
          ops.seqIdThr,
          ops.alnLenThr,
          ops.covMode,
          ops.covThr)) {
    return 0;
  } else {
    return -1;
  }
}

} // namespace foldseek

} // namespace lahuta

#endif // NEW_LOCAL_ALIGN_HPP
