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


      std::cout << "----------------------------------------" << std::endl;
      Mapper mapper(query, target, AR.res);
      mapper.map();
      mapper.print();
      std::cout << "----------------------------------------" << std::endl;

      std::cout << "Q data: " << AR.res.qStartPos << " " << AR.res.qEndPos << " " << AR.res.qLen << std::endl;
      std::cout << "T data: " << AR.res.dbStartPos << " " << AR.res.dbEndPos << " " << AR.res.dbLen << std::endl;

      std::cout << "Q Alignment: " << AR.query_alignment() << std::endl;
      std::cout << "T Alignment: " << AR.target_alignment() << std::endl;

      auto non_d_indices = get_non_deletion_indices(AR.res.backtrace); // query
      auto non_i_indices = get_non_insertion_indices(AR.res.backtrace); // target

      auto lddt_result = AR.scores->get_lddt_score();
      auto alig_result = AR.scores->get_alignment_score();
      auto prob = CalcProbTP::calculate(AR.res.score);

      std::cout << "----------------------------------------" << std::endl;
      std::cout << "RMSD: " << alig_result.rmsd << std::endl;
      std::cout << "LDDT: " << lddt_result.avgLddtScore << std::endl;
      std::cout << "TMscore: " << alig_result.tmscore << std::endl;
      std::cout << "QTMscore: " << AR.scores->get_query_score().tmscore << std::endl;
      std::cout << "TTMscore: " << AR.scores->get_target_score().tmscore << std::endl;
      std::cout << "Prob: " << prob << std::endl;
      std::cout << "----------------------------------------" << std::endl;
    }
  }
  return 0;
};
