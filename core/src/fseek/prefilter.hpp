#ifndef LAHUTA_PREFILTER_HPP
#define LAHUTA_PREFILTER_HPP

#include <climits>
#include <optional>

#include <IndexTable.h>
#include <QueryMatcher.h>
#include <Sequence.h>
#include <substitution_matrix.hpp>
#include <tantan.h>

#include "fseek/seq.hpp"
#include "score_matrix.hpp"

namespace lahuta {

// clang-format off
enum class SeqType { AminoAcid, Nucleotide, HMM };

struct PrefilterOptions {
  bool use_prefilter{true};

  BaseMatrix *kmerSubMat;
  BaseMatrix *ungappedSubMat;
  ScoreMatrix _2merSubMatrix;
  ScoreMatrix _3merSubMatrix;
  SequenceLookup *sequenceLookup;

  int alphabetSize;

  bool maskMode{};
  bool maskLowerCaseMode{true};
  float maskProb{0.99995};

  int kmerSize{6};
  int kmerThr{};
  bool spacedKmer{true};
  std::string spacedKmerPattern{};
  bool takeOnlyBestKmer{false};

  SeqType querySeqType = SeqType::AminoAcid;
  SeqType targetSeqType = SeqType::AminoAcid;
  int targetSearchMode{};

  float sensitivity{9.5};
  unsigned short maxSeqLen = std::numeric_limits<unsigned short>::max();
  unsigned int diagonalScoring{1};
  unsigned int minDiagScoreThr{30};
  bool aaBiasCorrection{true};
  float aaBiasCorrectionScale{0.15};
  float covThr{};
  int covMode{};

  size_t maxResListLen;
  // int exactKmerMatching;
  // bool includeIdentical;
};

namespace foldseek {

class DbInfo {
public:
  DbInfo(size_t dbFrom, size_t dbTo, unsigned int effectiveKmerSize, SeqCollection &seqData) {
    tableSize = 0;
    aaDbSize = 0;
    auto dbSize = dbTo - dbFrom;
    sequenceOffsets = new size_t[dbSize];
    sequenceOffsets[0] = 0;
    std::cout << "DbInfo: " << dbFrom << " " << dbTo << " " << effectiveKmerSize << " " << seqData.size() << "\n";
    for (size_t id = dbFrom; id < dbTo; id++) {
      const int seqLen = seqData[id].Seq3Di.size();
      aaDbSize += seqLen;
      size_t idFromNull = (id - dbFrom);
      if (id < dbSize - 1) {
        sequenceOffsets[idFromNull + 1] = sequenceOffsets[idFromNull] + seqLen;
      }
      if (Util::overlappingKmers(seqLen, effectiveKmerSize > 0)) {
        tableSize += 1;
      }
    }
  }

  ~DbInfo() { delete[] sequenceOffsets; }

  size_t tableSize;
  size_t aaDbSize;
  size_t *sequenceOffsets;
};

inline char *getScoreLookup(BaseMatrix &matrix) {
  char *idScoreLookup = NULL;
  idScoreLookup = new char[matrix.alphabetSize];
  for (int aa = 0; aa < matrix.alphabetSize; aa++) {
    short score = matrix.subMatrix[aa][aa];
    if (score > SCHAR_MAX || score < SCHAR_MIN) {
      Debug(Debug::WARNING) << "Truncating substitution matrix diagonal score!";
    }
    idScoreLookup[aa] = (char)score;
  }
  return idScoreLookup;
}

inline void fill_database(
    IndexTable *indexTable, SequenceLookup **maskedLookup, SequenceLookup **unmaskedLookup,
    BaseMatrix &subMat, ScoreMatrix &three, ScoreMatrix &two, Sequence *seq, SeqCollection *dbr,
    size_t dbFrom, size_t dbTo, int kmerThr, bool mask, bool maskLowerCaseMode, float maskProb,
    int targetSearchMode) {

  Debug(Debug::INFO) << "Index table: counting k-mers\n";

  const bool isProfile = Parameters::isEqualDbtype(seq->getSeqType(), Parameters::DBTYPE_HMM_PROFILE);
  const bool isTargetSimiliarKmerSearch = isProfile || targetSearchMode;
  dbTo = std::min(dbTo, dbr->size());
  size_t dbSize = dbTo - dbFrom;
  DbInfo *info = new DbInfo(dbFrom, dbTo, seq->getEffectiveKmerSize(), *dbr);

  SequenceLookup *sequenceLookup;
  if (unmaskedLookup != NULL && maskedLookup == NULL) {
    *unmaskedLookup = new SequenceLookup(dbSize, info->aaDbSize);
    sequenceLookup = *unmaskedLookup;
  } else if (unmaskedLookup == NULL && maskedLookup != NULL) {
    *maskedLookup = new SequenceLookup(dbSize, info->aaDbSize);
    sequenceLookup = *maskedLookup;
  } else if (unmaskedLookup != NULL && maskedLookup != NULL) {
    *unmaskedLookup = new SequenceLookup(dbSize, info->aaDbSize);
    *maskedLookup = new SequenceLookup(dbSize, info->aaDbSize);
    sequenceLookup = *maskedLookup;
  } else {
    Debug(Debug::ERROR) << "This should not happen\n";
    EXIT(EXIT_FAILURE);
  }

  // need to prune low scoring k-mers through masking
  ProbabilityMatrix *probMatrix = NULL;
  if (maskedLookup != NULL) {
    probMatrix = new ProbabilityMatrix(subMat);
  }

  // identical scores for memory reduction code
  char *idScoreLookup = getScoreLookup(subMat);

  size_t maskedResidues = 0;
  size_t totalKmerCount = 0;

  unsigned int thread_idx = 0;
  Indexer idxer(static_cast<unsigned int>(indexTable->getAlphabetSize()), seq->getKmerSize());
  Sequence s(
      seq->getMaxLen(),
      seq->getSeqType(),
      &subMat,
      seq->getKmerSize(),
      seq->isSpaced(),
      false,
      true,
      seq->getUserSpacedKmerPattern());

  KmerGenerator *generator = NULL;
  if (isTargetSimiliarKmerSearch) {
    generator = new KmerGenerator(seq->getKmerSize(), indexTable->getAlphabetSize(), kmerThr);
    if (isProfile) {
      generator->setDivideStrategy(s.profile_matrix);
    } else {
      generator->setDivideStrategy(&three, &two);
    }
  }

  unsigned int *buffer = static_cast<unsigned int *>(malloc(seq->getMaxLen() * sizeof(unsigned int)));
  unsigned int bufferSize = seq->getMaxLen();

  // FIX: Process using multiple threads
  for (size_t id = dbFrom; id < dbTo; id++) {
    s.resetCurrPos();

    SeqData seqData = dbr->get_data().at(id);
    s.mapSequence(id, id, seqData.Seq3Di.c_str(), seqData.Seq3Di.size());
    if (s.getMaxLen() >= bufferSize) {
      buffer = static_cast<unsigned int *>(realloc(buffer, s.getMaxLen() * sizeof(unsigned int)));
      bufferSize = seq->getMaxLen();
    }

    // count similar or exact k-mers based on sequence type
    if (isTargetSimiliarKmerSearch) {
      // Find out if we should also mask profiles
      totalKmerCount += indexTable->addSimilarKmerCount(&s, generator);
      unsigned char *seq = (isProfile) ? s.numConsensusSequence : s.numSequence;
      if (unmaskedLookup != NULL) {
        (*unmaskedLookup)->addSequence(seq, s.L, id - dbFrom, info->sequenceOffsets[id - dbFrom]);
      } else if (maskedLookup != NULL) {
        (*maskedLookup)->addSequence(seq, s.L, id - dbFrom, info->sequenceOffsets[id - dbFrom]);
      }
    } else {
      // Do not mask if column state sequences are used
      if (unmaskedLookup != NULL) {
        (*unmaskedLookup)->addSequence(s.numSequence, s.L, id - dbFrom, info->sequenceOffsets[id - dbFrom]);
      }
      if (mask == true) {
        // s.print();
        maskedResidues += tantan::maskSequences(
            (char *)s.numSequence,
            (char *)(s.numSequence + s.L),
            50 /*options.maxCycleLength*/,
            probMatrix->probMatrixPointers,
            0.005 /*options.repeatProb*/,
            0.05 /*options.repeatEndProb*/,
            0.9 /*options.repeatOffsetProbDecay*/,
            0,
            0,
            maskProb /*options.minMaskProb*/,
            probMatrix->hardMaskTable);
      }

      if (maskLowerCaseMode == true
          && (Parameters::isEqualDbtype(s.getSequenceType(), Parameters::DBTYPE_AMINO_ACIDS)
              || Parameters::isEqualDbtype(s.getSequenceType(), Parameters::DBTYPE_NUCLEOTIDES))) {

        const char *charSeq = s.getSeqData();
        unsigned char maskLetter = subMat.aa2num[static_cast<int>('X')];
        for (int i = 0; i < s.L; i++) {
          bool isLowerCase = (islower(charSeq[i]));
          maskedResidues += isLowerCase;
          s.numSequence[i] = isLowerCase ? maskLetter : s.numSequence[i];
        }
      }
      if (maskedLookup != NULL) {
        (*maskedLookup)->addSequence(s.numSequence, s.L, id - dbFrom, info->sequenceOffsets[id - dbFrom]);
      }

      totalKmerCount += indexTable->addKmerCount(&s, &idxer, buffer, kmerThr, idScoreLookup);
    }
  }

  free(buffer);

  if (generator != NULL) {
    delete generator;
  }

  if (probMatrix != NULL) {
    delete probMatrix;
  }

  Debug(Debug::INFO) << "Index table: Masked residues: " << maskedResidues << "\n";
  if (totalKmerCount == 0) {
    Debug(Debug::ERROR) << "No k-mer could be extracted.\n"
                        << "Maybe the sequences length is less than 14 residues.\n";
    if (maskedResidues == true) {
      Debug(Debug::ERROR) << " or contains only low complexity regions.";
      Debug(Debug::ERROR) << "Use --mask 0 to deactivate the low complexity filter.\n";
    }
    EXIT(EXIT_FAILURE);
  }
  // dbr->remapData();

  // TODO find smart way to remove extrem k-mers without harming huge protein families
  //    size_t lowSelectiveResidues = 0;
  //    const float dbSize = static_cast<float>(dbTo - dbFrom);
  //    for(size_t kmerIdx = 0; kmerIdx < indexTable->getTableSize(); kmerIdx++){
  //        size_t res = (size_t) indexTable->getOffset(kmerIdx);
  //        float selectivityOfKmer = (static_cast<float>(res)/dbSize);
  //        if(selectivityOfKmer > 0.005){
  //            indexTable->getOffset()[kmerIdx] = 0;
  //            lowSelectiveResidues += res;
  //        }
  //    }
  //    Debug(Debug::INFO) << "Index table: Remove "<< lowSelectiveResidues <<" none selective residues\n";
  //    Debug(Debug::INFO) << "Index table: init... from "<< dbFrom << " to "<< dbTo << "\n";

  ctpl::thread_pool pool(std::thread::hardware_concurrency());
  std::cout << "Index table size: " << info->tableSize << "\n";
  indexTable->initMemory(info->tableSize, pool);
  indexTable->init(pool);

  delete info;

  Debug(Debug::INFO) << "Index table: fill\n";

  // FIX: Process using multiple threads
  {
    unsigned int thread_idx = 0;
    Sequence s(
        seq->getMaxLen(),
        seq->getSeqType(),
        &subMat,
        seq->getKmerSize(),
        seq->isSpaced(),
        false,
        true,
        seq->getUserSpacedKmerPattern());
    Indexer idxer(static_cast<unsigned int>(indexTable->getAlphabetSize()), seq->getKmerSize());
    IndexEntryLocalTmp *buffer = static_cast<IndexEntryLocalTmp *>(malloc(seq->getMaxLen() * sizeof(IndexEntryLocalTmp)));
    size_t bufferSize = seq->getMaxLen();
    KmerGenerator *generator = NULL;
    if (isTargetSimiliarKmerSearch) {
      generator = new KmerGenerator(seq->getKmerSize(), indexTable->getAlphabetSize(), kmerThr);
      if (isProfile) {
        generator->setDivideStrategy(s.profile_matrix);
      } else {
        generator->setDivideStrategy(&three, &two);
      }
    }

    for (size_t id = dbFrom; id < dbTo; id++) {
      s.resetCurrPos();

      if (isTargetSimiliarKmerSearch) {
        SeqData seqData = dbr->get_data().at(id);
        s.mapSequence(id - dbFrom, id, seqData.Seq3Di.c_str(), dbr->get_data().at(id).Seq3Di.size());
        indexTable->addSimilarSequence(&s, generator, &buffer, bufferSize, &idxer);
      } else {
        const auto &[seqData, seqLen] = sequenceLookup->getSequence(id - dbFrom);
        auto pair = sequenceLookup->getSequence(id - dbFrom);
        s.mapSequence(id - dbFrom, id, pair);
        indexTable->addSequence(&s, &idxer, &buffer, bufferSize, kmerThr, idScoreLookup);
      }
    }

    if (generator != NULL) {
      delete generator;
    }

    free(buffer);
  }
  if (idScoreLookup != NULL) {
    delete[] idScoreLookup;
  }

  indexTable->revertPointer();
  indexTable->sortDBSeqLists(pool);

}

} // namespace foldseek

inline void fill_database(
    IndexTable *indexTable, PrefilterOptions &fo, Sequence *seq, SeqCollection *dbr, size_t dbFrom,
    size_t dbTo, int kmerThr) {

  SequenceLookup **unmaskedLookup = fo.maskMode == 0 && fo.maskLowerCaseMode == 0 ? &fo.sequenceLookup : NULL;
  SequenceLookup **maskedLookup = fo.maskMode == 1 || fo.maskLowerCaseMode == 1 ? &fo.sequenceLookup : NULL;

  // clang-format off
  foldseek::fill_database(
      indexTable, maskedLookup, unmaskedLookup,
      *fo.kmerSubMat, fo._3merSubMatrix, fo._2merSubMatrix,
      seq, dbr, dbFrom, dbTo, kmerThr,
      fo.maskMode, fo.maskLowerCaseMode, fo.maskProb, fo.targetSearchMode
  );
  // clang-format on
}

struct KmerThreshold {
  SeqType sequenceType;
  int kmerSize;
  float base;
  float sensPerStep;
};

inline std::vector<KmerThreshold> external_thr = {{SeqType::AminoAcid, 7, 197.0, 11.22}};
inline int getKmerThreshold(const float sensitivity, const int kmerSize) {
  float kmerThrBest = FLT_MAX;
  for (size_t i = 0; i < external_thr.size(); i++) {
    if (kmerSize == external_thr[i].kmerSize && external_thr[i].sequenceType == SeqType::AminoAcid) {
      return external_thr[i].base - (external_thr[i].sensPerStep * sensitivity);
    }
  }
  if (kmerSize == 5) {
    float base = 160.75;
    kmerThrBest = base - (sensitivity * 12.75);
  } else if (kmerSize == 6) {
    float base = 163.2;
    kmerThrBest = base - (sensitivity * 8.917);
  } else if (kmerSize == 7) {
    float base = 186.15;
    kmerThrBest = base - (sensitivity * 11.22);
  } else {
    Debug(Debug::ERROR) << "The k-mer size " << kmerSize << " is not valid\n";
    EXIT(EXIT_FAILURE);
  }
  return static_cast<int>(kmerThrBest);
}

using Hits = std::vector<int>;

class SeqFilter {
public:
  SeqFilter(SeqCollection &queries_, SeqCollection &targets_, PrefilterOptions &ops_)
      : targets(targets_), ops(ops_) {
    max_seq_len = queries_.max_length > targets_.max_length ? queries_.max_length : targets_.max_length;
    ops.maxResListLen = targets.size();
  }

  std::unique_ptr<IndexTable>
  get_index_table(std::optional<int> from = std::nullopt, std::optional<int> to = std::nullopt) {
    return build_index(from.value_or(0), to.value_or(targets.size()));
  }

  void build_index(std::optional<int> from = std::nullopt, std::optional<int> to = std::nullopt) {
    std::cout << "Building index table\n";
    indexTable = build_index(from.value_or(0), to.value_or(targets.size()));
    matcher = create_query_matcher();
  }

  void print_statistics(IndexTable *it) { it->printStatistics(ops.kmerSubMat->num2aa); }
  int get_max_seq_len() { return max_seq_len; }

  Hits filter(SeqData &Q) {

    Sequence seq = map_sequence(Q);
    set_seq_matrix(seq);

    size_t target_seq_id = UINT_MAX;
    std::pair<hit_t *, size_t> prefilter_result = matcher->matchQuery(&seq, target_seq_id, false);

    Hits hits;
    for (size_t i = 0; i < prefilter_result.second; i++) {
      hit_t *res = prefilter_result.first + i;

      // COV_MODE_BIDIRECTIONAL:0, COV_MODE_QUERY:2, COV_MODE_LENGTH_SHORTER:5
      if (ops.covThr > 0.0 && (ops.covMode == 0 || ops.covMode == 2 || ops.covMode == 5)) {
        if (!Util::canBeCovered(ops.covThr, ops.covMode, Q.size(), targets[res->seqId].size())) continue;
      }

      hits.push_back(res->seqId);
    }

    return hits;
  }

private:
  std::unique_ptr<QueryMatcher> create_query_matcher() {
    return std::make_unique<QueryMatcher>(
        indexTable.get(),
        ops.sequenceLookup,
        ops.kmerSubMat,
        ops.ungappedSubMat,
        ops.kmerThr,
        ops.kmerSize,
        targets.size(), // dbSize
        max_seq_len,
        ops.maxResListLen,
        ops.aaBiasCorrection,
        ops.aaBiasCorrectionScale,
        ops.diagonalScoring,
        ops.minDiagScoreThr,
        ops.takeOnlyBestKmer,
        ops.targetSeqType == SeqType::Nucleotide);
  }

  void set_seq_matrix(Sequence &seq) {
    if (seq.profile_matrix != NULL) {
      matcher->setProfileMatrix(seq.profile_matrix);
    } else if (ops._3merSubMatrix.isValid() && ops._2merSubMatrix.isValid()) {
      matcher->setSubstitutionMatrix(&ops._3merSubMatrix, &ops._2merSubMatrix);
    } else {
      matcher->setSubstitutionMatrix(NULL, NULL);
    }
  }

  Sequence map_sequence(const SeqData &Q, bool map = true) {
    Sequence seq = make_sequence(Q.size(), static_cast<int>(ops.querySeqType));
    seq.mapSequence(0, 0, Q.Seq3Di.c_str(), Q.Seq3Di.size());
    return seq;
  }

  Sequence make_sequence(int size, int seq_type) {
    return Sequence(
        size,
        seq_type,
        ops.kmerSubMat,
        ops.kmerSize,
        ops.spacedKmer,
        ops.aaBiasCorrection,
        true,
        ops.spacedKmerPattern);
  }

  [[nodiscard]] std::unique_ptr<IndexTable> build_index(int from, int to) {

    ops.kmerSubMat = create_base_matrix(CTSM8);
    ops.ungappedSubMat = create_base_matrix(CTSM2);
    ops.alphabetSize = ops.kmerSubMat->alphabetSize;

    ops.kmerThr = getKmerThreshold(ops.sensitivity, ops.kmerSize);
    ops.kmerSubMat->alphabetSize -= 1;
    ops._2merSubMatrix = calculate_score_matrix<2, 20>(CTSM8);
    ops._3merSubMatrix = calculate_score_matrix<3, 20>(CTSM8);
    ops.kmerSubMat->alphabetSize = ops.alphabetSize;
    auto it = std::make_unique<IndexTable>(ops.alphabetSize - 1, ops.kmerSize, false);
    fill_database(it.get(), from, to);

    return std::move(it);
  }

  void fill_database(IndexTable *it, size_t from, size_t to) {

    SequenceLookup **unmasked = !ops.maskMode && !ops.maskLowerCaseMode ? &ops.sequenceLookup : nullptr;
    SequenceLookup **masked = ops.maskMode || ops.maskLowerCaseMode ? &ops.sequenceLookup : nullptr;

    Sequence seq = make_sequence(ops.maxSeqLen, static_cast<int>(ops.targetSeqType));

    foldseek::fill_database(
        it,
        masked,
        unmasked,
        *ops.kmerSubMat,
        ops._3merSubMatrix,
        ops._2merSubMatrix,
        &seq,
        &targets,
        from,
        to,
        ops.kmerThr,
        ops.maskMode,
        ops.maskLowerCaseMode,
        ops.maskProb,
        ops.targetSearchMode);
  }

  SeqCollection &targets;
  PrefilterOptions &ops;
  std::unique_ptr<IndexTable> indexTable;
  /// The maximum possible sequence length
  int max_seq_len;
  std::unique_ptr<QueryMatcher> matcher;
};

} // namespace lahuta

#endif // LAHUTA_PREFILTER_HPP
