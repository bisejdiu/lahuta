#include "CalcProbTP.h"
#include "fseek/align.hpp"
#include "fseek/ops.hpp"

#include "lahuta.hpp"
#include "mapper.hpp"
#include "prefilter.hpp"
#include "seq_aligner_builder.hpp"

using namespace lahuta;

int main(int argc, char const *argv[]) {

  FoldSeekOps ops;
  PrefilterOptions pf_ops;

  // std::string source_file_name = "2rh1.cif";
  // std::string target_file_name = "3sn6.cif";

  // std::vector<std::string> file_names = {source_file_name, target_file_name};

  // SeqCollection queries = extract_all(ops, {file_names[0]});
  // SeqCollection targets = extract_all(ops, {file_names[1]});

  std::vector<std::string> q_files = {
      "/Users/bsejdiu/data/gpcrs/small/a/4nc3.cif",
      "/Users/bsejdiu/data/gpcrs/small/a/6drz.cif",
      "/Users/bsejdiu/data/gpcrs/small/a/7srr.cif",
      "/Users/bsejdiu/data/gpcrs/small/a/7um4.cif",
      "/Users/bsejdiu/data/gpcrs/small/a/7voe.cif",
      "/Users/bsejdiu/data/gpcrs/small/a/7wc4.cif",
      "/Users/bsejdiu/data/gpcrs/small/a/7xt8.cif",
      "/Users/bsejdiu/data/gpcrs/small/a/7xtb.cif",
      "/Users/bsejdiu/data/gpcrs/small/a/8dpf.cif",
      "/Users/bsejdiu/data/gpcrs/small/a/8w8b.cif"};

  std::vector<std::string> t_files = {
      "/Users/bsejdiu/data/gpcrs/small/b/2r4r.cif",
      "/Users/bsejdiu/data/gpcrs/small/b/2r4s.cif",
      "/Users/bsejdiu/data/gpcrs/small/b/3nya.cif",
      "/Users/bsejdiu/data/gpcrs/small/b/3sn6.cif",
      "/Users/bsejdiu/data/gpcrs/small/b/4qkx.cif",
      "/Users/bsejdiu/data/gpcrs/small/b/6e67.cif",
      "/Users/bsejdiu/data/gpcrs/small/b/6prz.cif",
      "/Users/bsejdiu/data/gpcrs/small/b/6ps2.cif",
      "/Users/bsejdiu/data/gpcrs/small/b/7dhi.cif",
      "/Users/bsejdiu/data/gpcrs/small/b/8jj8.cif"};

  SeqCollection queries = extract_all(ops, q_files);
  SeqCollection targets = extract_all(ops, t_files);

  std::cout << "Queries: " << queries.size() << std::endl;
  std::cout << "Targets: " << targets.size() << std::endl;

  if (!queries.success || !targets.success) {
    std::cerr << "Failed to process files" << std::endl;
    return 1;
  }

  SeqFilter seq_filter(queries, targets, pf_ops); // TODO: make pf_opts optional?
  seq_filter.build_index();

  SeqAlignerBuilder builder(ops);

  ops.altAlignment = 0;
  int query_count = 0;
  std::unique_ptr<SeqAligner> aligner = builder.build(queries, targets);

  int run_count = 0;
  for (auto &query : queries) {
    std::cout << "Q: " << query.SeqAA << std::endl;
    Hits hits = seq_filter.filter(query);

    for (const auto &hit_id : hits) {
      SeqData target = targets[hit_id];
      run_count++;
      std::cout << "T: " << target.SeqAA << std::endl;
      auto AR = aligner->align(query, target);
      if (!AR.success) continue;

      Mapper mapper(query, target); // , AR.ar.front());
      mapper.map(AR.ar.front());
      mapper.print();

      // FIX: building from mol changes it and crashes on the 2nd iteration bc bonds have been added
      /*Luni luni = Luni::from_RDKit(query.mol);*/
      Luni luni = Luni::build(query.st);
      /*Luni luni(query.file_name);*/
      std::cout << "Helo from Luni: " << luni.get_molecule().getNumAtoms() << std::endl;

      std::cout << "New result: \n" << Matcher::results_to_string(AR.ar) << std::endl;
      std::cout << "Q data: " << AR.ar.front().qStartPos << " " << AR.ar.front().qEndPos << std::endl;
      std::cout << "T data: " << AR.ar.front().dbStartPos << " " << AR.ar.front().dbEndPos << std::endl;

      /*std::cout << "Q Alignment: " << AR.query_alignment() << std::endl;*/
      /*std::cout << "T Alignment: " << AR.target_alignment() << std::endl;*/

      std::cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << std::endl;
      auto Scores = std::make_shared<AlignmentScores>(query, target, AR.ar[0]);
      std::cout << "RMSD: " << AR.scores.rmsd << std::endl;
      std::cout << "LDDT: " << AR.scores.avgLddtScore << std::endl;
      std::cout << "TMscore: " << AR.scores.tmscore << std::endl;
      std::cout << "QTMscore: " << Scores->get_query_score().tmscore << std::endl;
      std::cout << "TTMscore: " << Scores->get_target_score().tmscore << std::endl;
      std::cout << "Prob: " << AR.scores.prob << std::endl;
      std::cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << std::endl;
    }
  }

  std::cout << "Run count: " << run_count << std::endl;
  return 0;
};
