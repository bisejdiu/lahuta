#ifndef LAHUTA_SEQ_HPP
#define LAHUTA_SEQ_HPP

#include "GraphMol/RWMol.h"
#include "fseek/ops.hpp"
#include <GemmiWrapper.h>
#include <Sequence.h>
#include <SubstitutionMatrix.h>

namespace lahuta {

struct SeqData {
  /// Returns the number of residues in the sequence
  int size() const { return SeqAA.size(); }
  int size() { return SeqAA.size(); }
  float *x() { return CaData.data(); }
  float *y() { return CaData.data() + SeqAA.size(); }
  float *z() { return CaData.data() + 2 * SeqAA.size(); }

  bool operator<(const SeqData &rhs) const { return SeqAA < rhs.SeqAA; }

  std::unique_ptr<Sequence> map_3di(SubstitutionMatrix &matrix, FoldSeekOps &ops) {
    std::unique_ptr<Sequence> qSeq = build_sequence(matrix, ops);
    qSeq->mapSequence(0, 0, Seq3Di.c_str(), Seq3Di.size());
    return qSeq;
  }

  std::unique_ptr<Sequence> map_aa(SubstitutionMatrix &matrix, FoldSeekOps &ops) {
    std::unique_ptr<Sequence> qSeq = build_sequence(matrix, ops);
    qSeq->mapSequence(0, 0, SeqAA.c_str(), SeqAA.size());
    return qSeq;
  }

  std::unique_ptr<Sequence> build_sequence(SubstitutionMatrix &matrix, FoldSeekOps &ops) {
    return std::make_unique<Sequence>(ops.maxSeqLen, &matrix, ops.compBiasCorrection);
  }

  std::string Seq3Di;        // 3Di sequence
  std::string SeqAA;         // Amino acid sequence
  std::vector<float> CaData; // C-alpha coordinates
  std::string file_name;     // File name
  std::string chain_name;    // Chain name

  std::shared_ptr<RDKit::RWMol> mol = std::make_shared<RDKit::RWMol>();
};

struct SeqCollection {
  void add_data(const SeqData entity) { data.push_back(entity); }

  const std::vector<SeqData> &get_data() const { return data; }
  std::vector<SeqData> &get_data() { return data; }

  SeqData &operator[](size_t index) { return data[index]; };
  const SeqData &operator[](size_t index) const { return data[index]; };

  /// Returns the number of sequences in the collection
  size_t size() const { return data.size(); }

  typename std::vector<SeqData>::iterator begin() { return data.begin(); }
  typename std::vector<SeqData>::iterator end() { return data.end(); }
  typename std::vector<SeqData>::const_iterator cbegin() const { return data.cbegin(); }
  typename std::vector<SeqData>::const_iterator cend() const { return data.cend(); }

private:
  std::vector<SeqData> data;

public:
  /// Maximum sequence length
  int max_length{};
  /// Sum of all sequence lengths
  int total_length{};
  /// Indicates if the file was read successfully
  bool success{};
};

} // namespace lahuta

#endif // LAHUTA_SEQ_HPP
