#include "GraphMol/RWMol.h"
#include "lahuta.hpp"
#include "mapper.hpp"
#include "processor.hpp"
#include "residues.hpp"
#include "seq.hpp"
#include "seq_aligner.hpp"
#include "spdlog/sinks/stdout_color_sinks.h"

using namespace lahuta;

/*void print_processor(SeqData &query, SeqData &target, AlignmentResult &ar) {*/
/*  std::cout << "Result: \n"*/
/*            << Matcher::results_to_string(ar.ar) << "\n"*/
/*            << "Q Alignment: " << ar.query_alignment() << "\n"*/
/*            << "T Alignment: " << ar.target_alignment() << "\n";*/
/*}*/

void mapping_processor(SeqData &query, SeqData &target, AlignmentResult &ar) {

  if (query.file_name == target.file_name && query.chain_name >= target.chain_name) return;
  /*std::cout << "Considering: " << query.file_name << "_" << query.chain_name << " : " << target.file_name*/
  /*          << "_" << target.chain_name << std::endl;*/

  /*if (ar.ar[0].eval > 10e-2) return;*/
  if (ar.ar[0].eval > 1) return;

  /*Mapper mapper(query, target, ar.ar.front());*/
  /*mapper.map();*/

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

  auto luni = Luni::build(query.st);
  auto mol = luni.get_molecule();
  std::cout << "Molecule: " << query.file_name << " " << mol.getNumAtoms() << std::endl;

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

  auto luni_t = Luni::build(target.st);
  auto mol_t = luni_t.get_molecule();
  std::cout << "Molecule: " << target.file_name << " " << mol_t.getNumAtoms() << std::endl;

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

std::vector<int> get_residue_numbers(const Residues &residues, const SeqData &sd) {
  std::vector<int> res_nums;
  /*res_nums.reserve(mol.getNumAtoms() / 10);*/
  for (const auto &residue : residues) {
    if (residue.chain_id != sd.chain_name) continue;
    res_nums.push_back(residue.number);
  }
  return res_nums;
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

  // 0. Choose the alignment result
  Matcher::result_t res = ar.ar.front();

  auto mapper = LahutaMapper(query, target);
  mapper.map(res);

  // 2. Map indices
  InteractionOptions opts{5.0};

  Interactions ic_q(mapper.get_luni(Mapping::Query), opts);
  Interactions ic_t(mapper.get_luni(Mapping::Target), opts);

  std::cout << "Q Ionic" << std::endl;
  auto _5 = ic_q.ionic();
  _5.sort_interactions();
  _5.print_interactions();

  std::cout << "T Ionic" << std::endl;
  Contacts _5_t = ic_t.ionic();
  _5_t.sort_interactions();
  _5_t.print_interactions();

  std::cout << "mapping check" << std::endl;
  /*nm.evaluate(_5, Mapper::MappingType::Query, _5_t, Mapper::MappingType::Target);*/
  mapper.evaluate(ic_q.hbond(), Mapping::Query, ic_t.hbond(), Mapping::Target);
  std::cout << "end" << std::endl;
}

struct LahutaOptions {
  bool use_prefilter = false;
  bool use_alt_alignment = false;
  bool allow_query_self_operations = false;
};

int main() {
  auto logger = spdlog::stdout_color_mt("console");
  spdlog::level::level_enum log_level = spdlog::level::info;
  spdlog::set_level(log_level);

  /*DirectoryHandler dir_handler("/Users/bsejdiu/data/PDB_ARCHIVE_UNCOMP", ".cif", true);*/
  /*DirectoryHandler dir_handler("/Users/bsejdiu/data/mini_pdb_uncomp", ".cif", true);*/

  /*auto pdb_targets = dir_handler.get_all_files();*/

  try {
    LahutaProcessor::FileList query_files = {
        /*"/Users/bsejdiu/data/gpcrs/small/a/8w8b.cif",*/
        "/Users/bsejdiu/projects/lahuta/cpp/build/4ami.cif",
        /*"/Users/bsejdiu/data/gpcrs/small/a/4nc3.cif",*/
        /*"/Users/bsejdiu/data/gpcrs/small/a/8w8b.cif",*/
    };

    LahutaProcessor::FileList query_files_t = {
        "/Users/bsejdiu/progs/foldseek/build/src/test/4ami.cif",
        "/Users/bsejdiu/progs/foldseek/build/src/test/4nc3.cif",
        "/Users/bsejdiu/progs/foldseek/build/src/test/8w8b.cif",
    };

    LahutaProcessor::FileList target_files = {
        "/Users/bsejdiu/data/gpcrs/small/b/2r4r.cif",
        /*"/Users/bsejdiu/data/gpcrs/small/b/2r4s.cif",*/
        /*"/Users/bsejdiu/data/gpcrs/small/b/3nya.cif",*/
        /*"/Users/bsejdiu/data/gpcrs/small/b/8jj8.cif"*/
    };

    LahutaProcessor::FileList target_files_t = {
        "/Users/bsejdiu/progs/foldseek/build/src/test2/4ami.cif",
        "/Users/bsejdiu/progs/foldseek/build/src/test2/8w8b.cif",
    };

    LahutaProcessor::ProcessingConfig config{
        .query_chunk_size = 10,
        .target_chunk_size = 20000,
        .allow_query_self_operations = true,
        /*.file_extension = ".cif.gz",*/
        .file_extension = ".cif",
        .recursive_search = true,
    };

    // FIX: test with using alternative alignments
    FoldSeekOps ops;
    PrefilterOptions pf_ops;
    pf_ops.use_prefilter = true;

    std::unique_ptr<LahutaAligner> aligner = std::make_unique<LahutaAligner>(ops, pf_ops);
    aligner->set_computer(mapping_processor_w);
    /*aligner->set_computer(mapping_processor);*/

    LahutaProcessor processor(config); // , ops, pf_ops);
    processor.set_runner(std::move(aligner));
    /*processor.process_files(query_files);*/
    processor.process_files(query_files, target_files);
    /*processor.process_files({"/Users/bsejdiu/projects/lahuta/cpp/build/4ami.cif"}, pdb_targets);*/

    /*processor.process_files(query_files, pdb_targets);*/
    /*processor.process_files(query_files_t, target_files_t);*/

  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
  }
}
