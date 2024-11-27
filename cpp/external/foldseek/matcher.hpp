#ifndef CUSTOM_UTILS_HPP
#define CUSTOM_UTILS_HPP

#include "itoa.h"
#include <EvalueComputation.h>
#include <string>

class Matcher {

public:
  struct result_t {
    unsigned int dbKey;
    int score;
    float qcov;
    float dbcov;
    float seqId;
    double eval;
    unsigned int alnLength;
    int qStartPos;
    int qEndPos;
    unsigned int qLen;
    int dbStartPos;
    int dbEndPos;
    unsigned int dbLen;
    int queryOrfStartPos;
    int queryOrfEndPos;
    int dbOrfStartPos;
    int dbOrfEndPos;
    std::string backtrace;
    result_t(unsigned int dbkey, int score, float qcov, float dbcov,
             float seqId, double eval, unsigned int alnLength, int qStartPos,
             int qEndPos, unsigned int qLen, int dbStartPos, int dbEndPos,
             unsigned int dbLen, int queryOrfStartPos, int queryOrfEndPos,
             int dbOrfStartPos, int dbOrfEndPos, std::string backtrace)
        : dbKey(dbkey), score(score), qcov(qcov), dbcov(dbcov), seqId(seqId),
          eval(eval), alnLength(alnLength), qStartPos(qStartPos),
          qEndPos(qEndPos), qLen(qLen), dbStartPos(dbStartPos),
          dbEndPos(dbEndPos), dbLen(dbLen), queryOrfStartPos(queryOrfStartPos),
          queryOrfEndPos(queryOrfEndPos), dbOrfStartPos(dbOrfStartPos),
          dbOrfEndPos(dbOrfEndPos), backtrace(backtrace){};

    result_t(unsigned int dbkey, int score, float qcov, float dbcov,
             float seqId, double eval, unsigned int alnLength, int qStartPos,
             int qEndPos, unsigned int qLen, int dbStartPos, int dbEndPos,
             unsigned int dbLen, std::string backtrace)
        : dbKey(dbkey), score(score), qcov(qcov), dbcov(dbcov), seqId(seqId),
          eval(eval), alnLength(alnLength), qStartPos(qStartPos),
          qEndPos(qEndPos), qLen(qLen), dbStartPos(dbStartPos),
          dbEndPos(dbEndPos), dbLen(dbLen), queryOrfStartPos(-1),
          queryOrfEndPos(-1), dbOrfStartPos(-1), dbOrfEndPos(-1),
          backtrace(backtrace){};

    result_t(){};

    static void swapResult(result_t &res, EvalueComputation &evaluer,
                           bool hasBacktrace) {
      double rawScore = evaluer.computeRawScoreFromBitScore(res.score);
      res.eval = evaluer.computeEvalue(rawScore, res.dbLen);

      unsigned int qstart = res.qStartPos;
      unsigned int qend = res.qEndPos;
      unsigned int qLen = res.qLen;
      res.qStartPos = res.dbStartPos;
      res.qEndPos = res.dbEndPos;
      res.qLen = res.dbLen;
      res.dbStartPos = qstart;
      res.dbEndPos = qend;
      res.dbLen = qLen;
      if (hasBacktrace) {
        for (size_t j = 0; j < res.backtrace.size(); j++) {
          if (res.backtrace.at(j) == 'I') {
            res.backtrace.at(j) = 'D';
          } else if (res.backtrace.at(j) == 'D') {
            res.backtrace.at(j) = 'I';
          }
        }
      }
    }

    static void protein2nucl(std::string &backtrace,
                             std::string &newBacktrace) {
      char buffer[256];
      for (size_t pos = 0; pos < backtrace.size(); pos++) {
        int cnt = 0;
        if (isdigit(backtrace[pos])) {
          cnt += Util::fast_atoi<int>(backtrace.c_str() + pos);
          while (isdigit(backtrace[pos])) {
            pos++;
          }
        }
        bool update = false;
        switch (backtrace[pos]) {
        case 'M':
        case 'D':
        case 'I':
          update = true;
          break;
        }
        if (update) {
          char *buffNext = Itoa::i32toa_sse2(cnt * 3, buffer);
          size_t len = buffNext - buffer;
          newBacktrace.append(buffer, len - 1);
          newBacktrace.push_back(backtrace[pos]);
        }
      }
    }
  };

  static bool compareHits(const result_t &first, const result_t &second) {
    if (first.eval != second.eval) {
      return first.eval < second.eval;
    }
    if (first.score != second.score) {
      return first.score > second.score;
    }
    if (first.dbLen != second.dbLen) {
      return first.dbLen < second.dbLen;
    }
    return first.dbKey < second.dbKey;
  }

  static int computeAlnLength(int qStart, int qEnd, int dbStart, int dbEnd) {
    return std::max(abs(qEnd - qStart), abs(dbEnd - dbStart)) + 1;
  }

  static std::string compressAlignment(const std::string &bt) {
    std::string ret;
    char state = 'M';
    size_t counter = 0;
    for (size_t i = 0; i < bt.size(); ++i) {
      if (bt[i] != state) {
        // we could leave this out if counter == 1
        // to save a few byte (~5% of total cigar strings)
        ret.append(SSTR(counter));
        ret.push_back(state);
        state = bt[i];
        counter = 1;
      } else {
        counter++;
      }
    }
    ret.append(SSTR(counter));
    ret.push_back(state);
    return ret;
  }

  static std::string uncompressAlignment(const std::string &cbt) {
    std::string bt;
    bt.reserve(cbt.size());
    size_t count = 0;
    for (size_t i = 0; i < cbt.size(); ++i) {
        char c = cbt[i];
        if (c >= '0' && c <= '9') {
            count = count * 10 + c - '0';
        } else {
            bt.append(count == 0 ? 1 : count, c);
            count = 0;
        }
    }
    return bt;
}

  static size_t resultToBuffer(char *buff1, const Matcher::result_t &result,
                               bool addBacktrace, bool compress = true,
                               bool addOrfPosition = false) {
    char *basePos = buff1;
    char *tmpBuff = Itoa::u32toa_sse2((uint32_t)result.dbKey, buff1);
    *(tmpBuff - 1) = '\t';
    tmpBuff = Itoa::i32toa_sse2(result.score, tmpBuff);
    *(tmpBuff - 1) = '\t';
    tmpBuff = Util::fastSeqIdToBuffer(result.seqId, tmpBuff);
    *(tmpBuff - 1) = '\t';
    tmpBuff += snprintf(tmpBuff, 32, "%.3E", result.eval);
    tmpBuff++;
    *(tmpBuff - 1) = '\t';
    tmpBuff = Itoa::i32toa_sse2(result.qStartPos, tmpBuff);
    *(tmpBuff - 1) = '\t';
    tmpBuff = Itoa::i32toa_sse2(result.qEndPos, tmpBuff);
    *(tmpBuff - 1) = '\t';
    tmpBuff = Itoa::i32toa_sse2(result.qLen, tmpBuff);
    *(tmpBuff - 1) = '\t';
    tmpBuff = Itoa::i32toa_sse2(result.dbStartPos, tmpBuff);
    *(tmpBuff - 1) = '\t';
    tmpBuff = Itoa::i32toa_sse2(result.dbEndPos, tmpBuff);
    *(tmpBuff - 1) = '\t';
    tmpBuff = Itoa::i32toa_sse2(result.dbLen, tmpBuff);
    if (addOrfPosition) {
      *(tmpBuff - 1) = '\t';
      tmpBuff = Itoa::i32toa_sse2(result.queryOrfStartPos, tmpBuff);
      *(tmpBuff - 1) = '\t';
      tmpBuff = Itoa::i32toa_sse2(result.queryOrfEndPos, tmpBuff);
      *(tmpBuff - 1) = '\t';
      tmpBuff = Itoa::i32toa_sse2(result.dbOrfStartPos, tmpBuff);
      *(tmpBuff - 1) = '\t';
      tmpBuff = Itoa::i32toa_sse2(result.dbOrfEndPos, tmpBuff);
    }
    if (addBacktrace == true) {
      if (compress) {
        *(tmpBuff - 1) = '\t';
        std::string compressedCigar = compressAlignment(result.backtrace);
        tmpBuff =
            strncpy(tmpBuff, compressedCigar.c_str(), compressedCigar.length());
        tmpBuff += compressedCigar.length() + 1;
      } else {
        *(tmpBuff - 1) = '\t';
        tmpBuff = strncpy(tmpBuff, result.backtrace.c_str(),
                          result.backtrace.length());
        tmpBuff += result.backtrace.length() + 1;
      }
    }
    *(tmpBuff - 1) = '\n';
    *(tmpBuff) = '\0';
    return tmpBuff - basePos;
  }
};

static bool alignmentCheckCriteria(Matcher::result_t &res, bool isIdentity,
                          double evalThr, double seqIdThr, int alnLenThr,
                          int covMode, float covThr) {
  const bool evalOk = (res.eval <= evalThr);    // -e
  const bool seqIdOK = (res.seqId >= seqIdThr); // --min-seq-id
  const bool covOK = Util::hasCoverage(covThr, covMode, res.qcov, res.dbcov);
  const bool alnLenOK = Util::hasAlignmentLength(alnLenThr, res.alnLength);

  // check first if it is identity
  if (isIdentity ||
      // general accaptance criteria
      (evalOk && seqIdOK && covOK && alnLenOK)) {
    return true;
  } else {
    return false;
  }
};

#endif // CUSTOM_UTILS_HPP
