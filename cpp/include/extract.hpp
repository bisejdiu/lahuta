#ifndef LAHUTA_EXTRACT_HPP
#define LAHUTA_EXTRACT_HPP

#include "fseek/ops.hpp"
#include <GemmiWrapper.h>
#include <LDDT.h>
#include <Sequence.h>
#include <SubstitutionMatrix.h>

#include <cassert>
#include <iostream>
#include <matcher.hpp>
#include <vector>

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
  /// Indicates if the file was read successfully
  bool success{};
  /// Maximum sequence length
  int max_length{};
  /// Sum of all sequence lengths
  int total_length{};
};

inline SeqCollection
process_file(const FoldSeekOps &ops, const SubstitutionMatrix &mat, const std::string &file_name) {
  SeqCollection seqDataList;

  Alphabet3Di::StructureTo3Di structureTo3Di;
  GemmiWrapper readStructure;

  bool readSuccess = readStructure.load(file_name, static_cast<GemmiWrapper::Format>(ops.inputFormat));
  if (!readSuccess) {
    std::cerr << "Error: Failed to read structure file " << file_name << ".\n";
    return seqDataList;
  }

  for (size_t ch = 0; ch < readStructure.chain.size(); ch++) {
    size_t chainStart = readStructure.chain[ch].first;
    size_t chainEnd = readStructure.chain[ch].second;
    size_t chainLen = chainEnd - chainStart;

    SeqData seqData;

    if (chainLen <= 3) {
      std::cout << "Skipping chain " << readStructure.names[ch] << " reason: TOO SHORT\n";
      continue;
    }
    bool allX = true;
    for (size_t pos = 0; pos < chainLen; pos++) {
      const char aa = readStructure.ami[chainStart + pos];
      if (aa != 'X' && aa != 'x') {
        allX = false;
        break;
      }
    }
    if (allX) {
      std::cout << "Skipping chain " << readStructure.names[ch] << " reason: NON-PROTEINS NOT SUPPORTED\n";
      continue;
    }

    std::string alphabet3di;
    std::string alphabetAA;

    char *states = structureTo3Di.structure2states(
        &readStructure.ca[chainStart],
        &readStructure.n[chainStart],
        &readStructure.c[chainStart],
        &readStructure.cb[chainStart],
        chainLen);

    for (size_t pos = 0; pos < chainLen; pos++) {
      // Handle masking based on B-factor threshold
      if (readStructure.ca_bfactor[chainStart + pos] < ops.maskBfactorThreshold) {
        alphabet3di += tolower(mat.num2aa[static_cast<int>(states[pos])]);
        alphabetAA += tolower(readStructure.ami[chainStart + pos]);
      } else {
        alphabet3di += mat.num2aa[static_cast<int>(states[pos])];
        alphabetAA += readStructure.ami[chainStart + pos];
      }

      // Collect C-alpha coordinates
      seqData.CaData.resize(3 * chainLen);
      auto &ca_pos = readStructure.ca[chainStart + pos];
      seqData.CaData[pos] = ca_pos.x;
      seqData.CaData[chainLen + pos] = ca_pos.y;
      seqData.CaData[2 * chainLen + pos] = ca_pos.z;
    }

    seqData.Seq3Di = alphabet3di;
    seqData.SeqAA = alphabetAA;
    seqData.file_name = file_name;
    seqData.chain_name = readStructure.chainNames[ch];

    if (seqData.Seq3Di.size() != seqData.SeqAA.size()) {
      std::cerr << "Error: Sequence length mismatch for chain " << readStructure.names[ch] << ".\n";
      continue;
    }

    seqDataList.add_data(seqData);
    seqDataList.total_length += seqData.SeqAA.size();
    if (seqData.SeqAA.size() > seqDataList.max_length) seqDataList.max_length = seqData.SeqAA.size();
  }

  if (!seqDataList.get_data().empty()) seqDataList.success = true;

  return seqDataList;
}

inline SeqCollection extract_all(FoldSeekOps &ops, const std::vector<std::string> &file_names) {
  SubstitutionMatrix mat(SubMatrix3Di, 2.0, ops.scoreBias);

  SeqCollection results;
  for (auto &file_name : file_names) {
    auto reas = process_file(ops, mat, file_name);
    for (auto &seqData : reas) {
      results.add_data(seqData);
      results.total_length += seqData.SeqAA.size();
      if (seqData.SeqAA.size() > results.max_length) results.max_length = seqData.SeqAA.size();
    }
  }

  if (!results.get_data().empty()) results.success = true;

  return results;
}

} // namespace lahuta

#endif // LAHUTA_EXTRACT_HPP
