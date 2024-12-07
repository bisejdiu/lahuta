#ifndef LAHUTA_EXTRACT_HPP
#define LAHUTA_EXTRACT_HPP

#include "convert.hpp"
#include "gemmi/model.hpp"
#include "seq.hpp"

#include "ctpl/ctpl_stl.h"
#include <cassert>
#include <iostream>
#include <matcher.hpp>
#include <vector>

namespace lahuta {

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

  // std::shared_ptr<RDKit::RWMol> mol = std::make_shared<RDKit::RWMol>();
  // RDKit::Conformer *conformer = new RDKit::Conformer();
  // create_RDKit_repr(*mol, readStructure.st, *conformer, false);
  // mol->updatePropertyCache(false);
  // mol->addConformer(conformer, true);

  for (size_t ch = 0; ch < readStructure.chain.size(); ch++) {
    size_t chainStart = readStructure.chain[ch].first;
    size_t chainEnd = readStructure.chain[ch].second;
    size_t chainLen = chainEnd - chainStart;

    SeqData seqData;

    if (chainLen <= 3) {
      /*std::cout << "Skipping chain " << readStructure.names[ch] << " reason: TOO SHORT\n";*/
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
      /*std::cout << "Skipping chain " << readStructure.names[ch] << " reason: NON-PROTEINS NOT
       * SUPPORTED\n";*/
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
    seqData.file_name =
        file_name.rfind('/') != std::string::npos ? file_name.substr(file_name.rfind('/') + 1) : file_name;
    seqData.chain_name = readStructure.chainNames[ch];

    if (seqData.Seq3Di.size() != seqData.SeqAA.size()) {
      std::cerr << "Error: Sequence length mismatch for chain " << readStructure.names[ch] << ".\n";
      continue;
    }

    seqData.st = readStructure.st;
    /*seqData.mol = mol;*/
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
    SeqCollection reas = process_file(ops, mat, file_name);
    for (auto &seqData : reas) {
      results.total_length += seqData.SeqAA.size();
      if (seqData.SeqAA.size() > results.max_length) results.max_length = seqData.SeqAA.size();
      results.add_data(std::move(seqData));
    }
  }

  if (!results.get_data().empty()) results.success = true;

  return results;
}

inline SeqCollection extract_all_parallel(FoldSeekOps &ops, const std::vector<std::string> &file_names) {
  SubstitutionMatrix mat(SubMatrix3Di, 2.0, ops.scoreBias);

  SeqCollection results;
  int num_threads = std::thread::hardware_concurrency();
  ctpl::thread_pool pool(num_threads);
  std::vector<std::future<SeqCollection>> futures;

  for (auto &file_name : file_names) {
    futures.emplace_back(
        pool.push([file_name, &ops, &mat](int id) { return process_file(ops, mat, file_name); }));
  }

  for (auto &f : futures) {
    auto reas = f.get();
    for (auto &seqData : reas) {
      results.total_length += seqData.SeqAA.size();
      if (seqData.SeqAA.size() > results.max_length) results.max_length = seqData.SeqAA.size();
      results.add_data(std::move(seqData));
    }
  }

  if (!results.get_data().empty()) results.success = true;

  return results;
}

} // namespace lahuta

#endif // LAHUTA_EXTRACT_HPP
