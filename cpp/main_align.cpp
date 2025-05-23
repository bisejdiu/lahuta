#include "GraphMol/RWMol.h"
#include "contacts/interactions.hpp"
#include "contacts/utils.hpp"
#include "file_system.hpp"
#include "lahuta.hpp"
#include "logging.hpp"
#include "mapper.hpp"
/*#include "processor.hpp"*/
#include "seq.hpp"
#include "seq_aligner.hpp"
/*#include "proc.hpp"*/
#include "procs/search_and_align.hpp"

// clang-format off

using namespace lahuta;

void mapping_processor(SeqData &query, SeqData &target, AlignmentResult &ar) {

  if (query.file_name == target.file_name && query.chain_name >= target.chain_name) return;

  if (ar.ar[0].eval > 10e-2) return;

  std::cout << "Mapping: " << query.file_name << "_" << query.chain_name << " : " //
            << target.file_name << "_" << target.chain_name << " : " << Matcher::results_to_string(ar.ar)
            << std::endl;

  // spdlog::warn("Mapping: {}_{} : {}_{} : {}", query.file_name, query.chain_name, target.file_name,
  // target.chain_name,
  //              Matcher::results_to_string(ar.ar));

  /*mapper.print();*/
  /*alignment_computers::print_result(query, target, ar);*/
}

void _mapping_processor_w(SeqData &query, SeqData &target, AlignmentResult &ar) {

  if (query.file_name == target.file_name && query.chain_name >= target.chain_name) return;
  if (ar.ar[0].eval > 1) return;

  alignment_computers::print_result(query, target, ar);

  auto luni = Luni::create(query.st);
  if (!luni.build_topology()) {
    Logger::get_logger()->warn("Failed building topology for {}!", luni.get_file_name());
  }
  auto mol = luni.get_molecule();
  std::cout << "Query Molecule: " << query.file_name << " " << mol.getNumAtoms() << std::endl;

  InteractionOptions opts{5.0};
  Interactions interactions(luni, opts);

  std::cout << "Ionic" << std::endl;
  auto _5 = interactions.ionic();
  _5.sort_interactions();
  _5.print_interactions();

  std::cout << "CationPi" << std::endl;
  auto _7 = interactions.cationpi();
  _7.sort_interactions();
  _7.print_interactions();

  std::cout << "PiStacking" << std::endl;
  auto _8 = interactions.pistacking();
  _8.sort_interactions();
  _8.print_interactions();

  auto luni_t = Luni::create(target.st);
  if (!luni_t.build_topology()) {
    Logger::get_logger()->warn("Failed building topology for {}!", luni_t.get_file_name());
  }
  auto mol_t = luni_t.get_molecule();
  std::cout << "Target Molecule: " << target.file_name << " " << mol_t.getNumAtoms() << std::endl;

  InteractionOptions opts_t{5.0};
  Interactions interactions_t(luni_t, opts_t);

  std::cout << "Ionic" << std::endl;
  auto _5_t = interactions_t.ionic();
  _5_t.sort_interactions();
  _5_t.print_interactions();

  std::cout << "CationPi" << std::endl;
  auto _7_t = interactions_t.cationpi();
  _7_t.sort_interactions();
  _7_t.print_interactions();

  std::cout << "PiStacking" << std::endl;
  auto _8_t = interactions_t.pistacking();
  _8_t.sort_interactions();
  _8_t.print_interactions();

  save_file("/Users/bsejdiu/projects/lahuta/cpp/build/aligned_q.pdb", query, ar);
  save_file("/Users/bsejdiu/projects/lahuta/cpp/build/aligned_t.pdb", target, ar);
}

Contacts compute_contacts(const Luni &luni, const InteractionOptions &opts) {
  Interactions interactions(luni, opts);
  auto ionic = interactions.ionic();
  auto cationpi = interactions.cationpi();
  auto pistacking = interactions.pistacking();
  auto hbond = interactions.hbond();
  auto weak_hbond = interactions.weak_hbond();
  auto hydrophobic = interactions.hydrophobic();
  auto halogen = interactions.halogen();
  auto metalic = interactions.metalic();

  Contacts contacts;
  contacts.add(ionic);
  contacts.add(cationpi);
  contacts.add(pistacking);
  contacts.add(hbond);
  contacts.add(weak_hbond);
  contacts.add(hydrophobic);
  contacts.add(halogen);
  contacts.add(metalic);

  contacts.sort_interactions();
  std::cout << "Total Contacts: " << contacts.size() << std::endl;
  return contacts;
}

void mapping_processor_w(SeqData &query, SeqData &target, AlignmentResult &ar) {

  using Mapping = TopologyMapper::MappingType;

  if (query.file_name == target.file_name && query.chain_name >= target.chain_name) return;
  if (ar.ar[0].eval > 1) return;

  std::cout << "Mapping: " << query.file_name << "_" << query.chain_name << " : " //
            << target.file_name << "_" << target.chain_name << " : "              //
            << Matcher::results_to_string(ar.ar) << std::endl;

  alignment_computers::print_result(query, target, ar);

  if (ar.ar.empty()) return;

  // 0. choose the alignment result
  Matcher::result_t res = ar.ar.front();

  auto mapper = LahutaMapper(query, target);
  mapper.map(res);

  InteractionOptions opts{5.0};

  Interactions ic_q(mapper.get_luni(Mapping::Query), opts);
  Interactions ic_t(mapper.get_luni(Mapping::Target), opts);

  auto contacts_q = ic_q.compute_contacts();
  auto contacts_t = ic_t.compute_contacts();

  auto count = mapper.evaluate(contacts_q, Mapping::Query, contacts_t, Mapping::Target);
  std::cout << "result: " << count << std::endl;
}

Contacts compute_neighbor_contacts(const Luni &luni, double cutoff) {
  Contacts contacts(&luni);
  auto grid = FastNS(luni.get_conformer().getPositions());
  auto ok = grid.build(cutoff);
  if (!ok) {
    throw std::runtime_error("Box dimension too small for the given cutoff.");
  }
  NSResults neighbors = grid.self_search();

  for (const auto &[pair, dist] : neighbors) {
    auto [i, j] = pair;
    const RDKit::Atom *a1 = luni.get_molecule().getAtomWithIdx(i);
    const RDKit::Atom *a2 = luni.get_molecule().getAtomWithIdx(j);
    if (are_residueids_close(luni.get_molecule(), *a1, *a2, 4)) continue;

    contacts.add(Contact(
        make_entity_id(EntityType::Atom, a1->getIdx()),
        make_entity_id(EntityType::Atom, a2->getIdx()),
        dist));

  }

  contacts.sort_interactions();
  std::cout << "Total Contacts: " << contacts.size() << std::endl;
  /*contacts.print_interactions();*/
  return contacts;
};

std::optional<Matcher::result_t> mapping_processor_w2(SeqData &query, SeqData &target, AlignmentResult &ar) {

  using Mapping = TopologyMapper::MappingType;

  if (query.file_name == target.file_name && query.chain_name >= target.chain_name) return std::nullopt;
  if (ar.ar[0].eval > 1) return std::nullopt;

  std::cout << "Mapping: " << query.file_name << "_" << query.chain_name << " : " //
            << target.file_name << "_" << target.chain_name << " : "              //
            << Matcher::results_to_string(ar.ar) << std::endl;

  alignment_computers::print_result(query, target, ar);

  if (ar.ar.empty()) return std::nullopt;

  Matcher::result_t res = ar.ar.front();

  auto mapper = LahutaMapper(query, target);
  mapper.map(res);

  InteractionOptions opts{5.0};
  /*opts.ionic = false;*/
  /*opts.hbond = false;*/
  /*opts.weak_hbond = true;*/
  /*opts.hydrophobic = false;*/
  /*opts.halogen = false;*/
  /*opts.metalic = false;*/
  /*opts.cationpi = false;*/
  /*opts.pistacking = false;*/

  Interactions ic_q(mapper.get_luni(Mapping::Query), opts);
  Interactions ic_t(mapper.get_luni(Mapping::Target), opts);

  auto contacts_q = ic_q.compute_contacts();
  auto contacts_t = ic_t.compute_contacts();
  contacts_q.print_interactions();
  contacts_t.print_interactions();

  std::cout << "Q contacts: " << contacts_q.size() << std::endl;
  std::cout << "T contacts: " << contacts_t.size() << std::endl;

  /*auto count = mapper.evaluate(contacts_q, Mapping::Query, contacts_t, Mapping::Target);*/

  std::cout << "returning result" << std::endl;
  auto qc = compute_neighbor_contacts(mapper.get_luni(Mapping::Query), 6.0);
  auto tc = compute_neighbor_contacts(mapper.get_luni(Mapping::Target), 6.0);

  std::cout << "--> Q contacts: " << qc.size() << std::endl;
  std::cout << "--> T contacts: " << tc.size() << std::endl;

  /*auto count = mapper.evaluate(contacts_q, Mapping::Query, contacts_t, Mapping::Target);*/
  auto count = mapper.evaluate(qc, Mapping::Query, tc, Mapping::Target);
  std::cout << "end: " << count << std::endl;

  return res;
}

std::optional<Matcher::result_t> mapping_processor_w3(SeqData &query, SeqData &target, AlignmentResult &ar) {

  using Mapping = TopologyMapper::MappingType;

  if (query.file_name == target.file_name && query.chain_name >= target.chain_name) return std::nullopt;
  if (ar.ar[0].eval > 1) return std::nullopt;

  std::cout << "Mapping: " << query.file_name << "_" << query.chain_name << " : "
            << target.file_name << "_" << target.chain_name << " : "
            << Matcher::results_to_string(ar.ar) << std::endl;

  std::cout << "xMapping: " << target.file_name << "_" << target.chain_name << " : " << ar.scores.rmsd << std::endl;
  alignment_computers::print_result(query, target, ar);

  if (ar.ar.empty()) return std::nullopt;

  Matcher::result_t res = ar.ar.front();

  std::cout << "Returning result" << std::endl;
  return res;
}


int main() {

  /*Logger::get_instance().set_log_level(Logger::LogLevel::Info);*/

  /*std::string dir = "/Users/bsejdiu/data/PDB_ARCHIVE_UNCOMP";*/
  /*std::string dir = "/Users/bsejdiu/data/gpcrs/small/a";*/
  /*std::string dir = "/Users/bsejdiu/data/mini_pdb_uncomp";*/

  // DirectoryHandler dir_handler(dir, ".cif", /* recursive= */ true);
  // std::vector<std::string> pdb_targets = dir_handler.get_all_files();

  std::vector<std::string> query_files_t = {
      "/Users/bsejdiu/progs/foldseek/build/src/test/4nc3.cif",
      /*"/Users/bsejdiu/progs/foldseek/build/src/test/4ami.cif",*/
      /*"/Users/bsejdiu/progs/foldseek/build/src/test/8w8b.cif",*/
  };

  std::vector<std::string> target_files = {
      "/Users/bsejdiu/data/gpcrs/small/b/2r4r.cif",
      "/Users/bsejdiu/progs/foldseek/build/src/test/4ami.cif",
      /*"/Users/bsejdiu/data/gpcrs/small/b/2r4s.cif",*/
      /*"/Users/bsejdiu/data/gpcrs/small/b/3nya.cif",*/
      "/Users/bsejdiu/data/gpcrs/small/b/8jj8.cif"
  };

  std::vector<std::string> target_files_t = {
      /*"/Users/bsejdiu/progs/foldseek/build/src/test2/4ami.cif",*/
      /*"/Users/bsejdiu/progs/foldseek/build/src/test2/8w8b.cif",*/
  };

  std::vector<std::string> qq = {"/Users/bsejdiu/projects/lahuta/cpp/build/2rh1.cif"};
  std::vector<std::string> tt = {"/Users/bsejdiu/projects/lahuta/cpp/build/3sn6_chainR.cif"};


  // {

  //   FoldSeekOps ops;
  //   PrefilterOptions pf_ops;
  //   pf_ops.use_prefilter = false;


  //   ProcessingConfig config{
  //       .query_chunk_size = 10,
  //       .target_chunk_size = 10000, // NOTE: the number of targets will be higher because each chain makes up a new target
  //       .allow_self_ops = false,
  //   };

  //   std::shared_ptr<LahutaAligner> aligner = std::make_shared<LahutaAligner>(ops, pf_ops);
  //   LahutaProcessor processor(aligner, config);
  //   processor.process(query_files_t, pdb_targets);

  //   auto & results = aligner->get_results();
  //   std::cout << "Results: " << results.size() << std::endl;

  //   for (auto &res : results) {
  //     auto &ars = res.results;
  //     for (auto &ar : ars) {
  //       std::cout << "Result: " << res.query->file_name << " : "
  //                 /*<< res.t_file << " : " << res.q_ch << " : " << res.t_ch //*/
  //                 << res.target->file_name << " : " << res.query->chain_name << " : " << res.target->chain_name
  //                 << " : " << Matcher::result_to_string(ar)
  //                 << std::endl;
  //     }
  //   }
  // }


  {

    FoldSeekOps ops;
    PrefilterOptions pf_ops;
    pf_ops.use_prefilter = false; // FIX: LDDT works only with prefilter true (I think this is due to some memory-related issue in the source code)
    ops.altAlignment = 8;


    ProcessingConfig config{
        .query_chunk_size = 10,
        .target_chunk_size = 20000, // NOTE: the number of targets will be higher because each chain makes up a new target
        .allow_self_ops = true,
    };

    std::shared_ptr<LahutaAligner> aligner = std::make_shared<LahutaAligner>(ops, pf_ops);
    LahutaProcessor processor(aligner, config);
    processor.process(query_files_t, target_files);

    auto & results = aligner->get_results();
    std::cout << "Results: " << results.size() << std::endl;

    for (auto &res : results) {
      auto &ars = res.results;
      for (auto &ar : ars) {
        std::cout << "Result: " << res.query->file_name << " : "
                  /*<< res.t_file << " : " << res.q_ch << " : " << res.t_ch //*/
                  << res.target->file_name << " : " << res.query->chain_name << " : " << res.target->chain_name
                  << " : " << Matcher::result_to_string(ar)
                  << std::endl;
      }
    }
  }

  {

    FoldSeekOps ops;
    PrefilterOptions pf_ops;
    pf_ops.use_prefilter = true;
    ops.altAlignment = 8;


    ProcessingConfig config{
        .query_chunk_size = 10,
        .target_chunk_size = 20000, // NOTE: the number of targets will be higher because each chain makes up a new target
        .allow_self_ops = false,
    };

    std::shared_ptr<LahutaAligner> aligner = std::make_shared<LahutaAligner>(ops, pf_ops);
    LahutaProcessor processor(aligner, config);
    processor.process(query_files_t, target_files);

    auto & results = aligner->get_results();
    std::cout << "Results: " << results.size() << std::endl;

    for (auto &res : results) {
      auto &ars = res.results;
      for (auto &ar : ars) {
        std::cout << "Result: " << res.query->file_name << " : "
                  << res.target->file_name << " : " << res.query->chain_name << " : " << res.target->chain_name
                  << " : " << Matcher::result_to_string(ar)
                  << std::endl;
      }
    }
  }

  {
    FoldSeekOps ops;
    PrefilterOptions pf_ops;
    pf_ops.use_prefilter = false;
    ops.altAlignment = 8;

    ProcessingConfig config{
        .query_chunk_size = 10,
        .target_chunk_size = 20000, // NOTE: the number of targets will be higher because each chain makes up a new target
        .allow_self_ops = true,
    };

    auto aligner = std::make_shared<LahutaAligner>(ops, pf_ops);
    LahutaProcessor processor(aligner, config);
    processor.process(target_files);

    auto & results = aligner->get_results();
    std::cout << "Results: " << results.size() << std::endl;

    for (auto &res : results) {
      auto &ars = res.results;
      for (auto &ar : ars) {
        std::cout << "Result: " << res.query->file_name << " : "
                  << res.target->file_name << " : " << res.query->chain_name << " : " << res.target->chain_name
                  << " : " << Matcher::result_to_string(ar)
                  << std::endl;
      }
    }
  }


}

