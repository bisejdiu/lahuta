#include "CalcProbTP.h"
#include "align.hpp"
#include "fseek/ops.hpp"

#include "prefilter.hpp"
#include "seq_aligner_builder.hpp"
#include "mapper.hpp"

using namespace lahuta;

// TODO: 1. Make prefilter optional if the number of alignments is small (benchmark)
int main(int argc, char const *argv[]) {

  FoldSeekOps ops;
  PrefilterOptions pf_ops;

  std::string source_file_name = "2rh1.cif";
  std::string target_file_name = "3sn6.cif";

  std::vector<std::string> fileNames = {source_file_name, target_file_name};

  SeqCollection queries = extract_all(ops, {fileNames[0]});
  SeqCollection targets = extract_all(ops, {fileNames[1]});

  std::cout << "Queries: " << queries.size() << std::endl;
  std::cout << "Targets: " << targets.size() << std::endl;

  SeqFilter seq_filter(queries, targets, pf_ops); // TODO: make pf_opts optional?
  seq_filter.build_index();

  SeqAlignerBuilder builder(ops);

  ops.altAlignment = 10;
  int query_count = 0;
  std::unique_ptr<SeqAligner> aligner = builder.build(queries, targets);

  for (auto &query : queries) {
    std::cout << "Q: " << query.SeqAA << std::endl;
    Hits hits = seq_filter.filter(query);

    for (const auto &hit_id : hits) {
      SeqData target = targets[hit_id];
      std::cout << "T: " << target.SeqAA << std::endl;
      auto AR = aligner->align(query, target, ops);
      if (!AR.success) continue;

      Mapper mapper(query, target, AR.ar.front());
      mapper.map();
      mapper.print();

      auto Scores = std::make_shared<AlignmentScores>(query, target, AR.ar[0]);

      std::cout << "New result: \n" << Matcher::results_to_string(AR.ar) << std::endl;
      std::cout << "Q data: " << AR.ar[0].qStartPos << " " << AR.ar[0].qEndPos << std::endl;
      std::cout << "T data: " << AR.ar[0].dbStartPos << " " << AR.ar[0].dbEndPos << std::endl;

      std::cout << "Q Alignment: " << AR.query_alignment() << std::endl;
      std::cout << "T Alignment: " << AR.target_alignment() << std::endl;

      std::cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << std::endl;
      std::cout << "RMSD: " << AR.scores.rmsd << std::endl;
      std::cout << "LDDT: " << AR.scores.avgLddtScore << std::endl;
      std::cout << "TMscore: " << AR.scores.tmscore << std::endl;
      std::cout << "QTMscore: " << Scores->get_query_score().tmscore << std::endl;
      std::cout << "TTMscore: " << Scores->get_target_score().tmscore << std::endl;
      std::cout << "Prob: " << AR.scores.prob << std::endl;
      std::cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << std::endl;
    }
  }
  return 0;
};
