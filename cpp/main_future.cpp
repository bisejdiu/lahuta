#include "GraphMol/RWMol.h"
#include "definitions.hpp"
#include "lahuta.hpp"
#include "mapper.hpp"
#include "processor.hpp"
#include "residues.hpp"
#include "seq.hpp"
#include "seq_aligner.hpp"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "common.hpp"
#include <stdexcept>

using namespace lahuta;

void mapping_processor(SeqData &query, SeqData &target, AlignmentResult &ar) {

  if (query.file_name == target.file_name && query.chain_name >= target.chain_name) return;
  /*std::cout << "Considering: " << query.file_name << "_" << query.chain_name << " : " << target.file_name*/
  /*          << "_" << target.chain_name << std::endl;*/

  /*if (ar.ar[0].eval > 10e-2) return;*/
  if (ar.ar[0].eval > 1) return;

  /*Mapper mapper(query, target, ar.ar.front());*/
  /*mapper.map();*/

  std::cout << "Mapping: " << query.file_name << "_" << query.chain_name << " : " //
            << target.file_name << "_" << target.chain_name << " : " 
            << Matcher::results_to_string(ar.ar) << std::endl;
  /*mapper.print();*/
  /*alignment_computers::print_result(query, target, ar);*/
}

void print_atoms(const RDKit::RWMol &mol, std::vector<int> atom_ids) {
  for (const auto atom_id: atom_ids) {
    auto atom = mol.getAtomWithIdx(atom_id);
    auto ri = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
    std::cout << "Atom: " << atom_id << " " << ri->getResidueName()  //
              << "-" << ri->getResidueNumber() << "-" << ri->getName()  // 
              << " :: " << ri->getResidueIndex() << "\n";
  }
}

std::vector<int> get_residue_numbers(const Residues &residues, const SeqData &sd) {
  std::vector<int> res_nums;
  /*res_nums.reserve(mol.getNumAtoms() / 10);*/
  for (const auto &residue: residues) {
    if (residue.chain_id != sd.chain_name) continue;
    res_nums.push_back(residue.number);
  }
  return res_nums;
}

void find_and_log_residue(const Luni &luni, const LuniMapper &luni_mapper, int atom_id) {

  auto &mol = luni.get_molecule();
  auto atom = mol.getAtomWithIdx(atom_id);
  auto ri = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
  auto r_ix = ri->getResidueIndex();
  auto r_id = ri->getResidueNumber();

  auto mapped_resid = luni_mapper.get_mapped_resid(atom_id);
  if (mapped_resid) {
    std::cout << "Residue index: "
              << r_id << " " << r_ix << " -> "
              << *mapped_resid << " "
              << ri->getName() << "-" << ri->getResidueName()
              << "-" << ri->getResidueNumber() <<  std::endl;
  }
}

void mapping_processor_w(SeqData &query, SeqData &target, AlignmentResult &ar) {

  if (query.file_name == target.file_name && query.chain_name >= target.chain_name) return;
  if (ar.ar[0].eval > 1) return;

  std::cout << "Mapping: " << query.file_name << "_" << query.chain_name << " : " //
            << target.file_name << "_" << target.chain_name << " : " 
            << Matcher::results_to_string(ar.ar) << std::endl;

  alignment_computers::print_result(query, target, ar);

  if (ar.ar.empty()) return;

  // 0. Choose the alignment result
  Matcher::result_t res = ar.ar.front();

  // 1. Initialize the mapper
  LuniMapper luni_mapper_q(query,  LuniMapper::MappingType::Query);
  LuniMapper luni_mapper_t(target, LuniMapper::MappingType::Target);

  luni_mapper_q.map(res);
  luni_mapper_t.map(res);

  // 2. Map indices
  // 795 797 798  --- 786 787
  find_and_log_residue(luni_mapper_q.get_luni(), luni_mapper_q, 795);
  find_and_log_residue(luni_mapper_q.get_luni(), luni_mapper_q, 797);
  find_and_log_residue(luni_mapper_q.get_luni(), luni_mapper_q, 798);
  find_and_log_residue(luni_mapper_q.get_luni(), luni_mapper_q, 786);
  find_and_log_residue(luni_mapper_q.get_luni(), luni_mapper_q, 787);

  // 530 532 533  --- 521 522
  find_and_log_residue(luni_mapper_t.get_luni(), luni_mapper_t, 530);
  find_and_log_residue(luni_mapper_t.get_luni(), luni_mapper_t, 532);
  find_and_log_residue(luni_mapper_t.get_luni(), luni_mapper_t, 533);
  find_and_log_residue(luni_mapper_t.get_luni(), luni_mapper_t, 521);
  find_and_log_residue(luni_mapper_t.get_luni(), luni_mapper_t, 522);
}

// TODO: 1. Implement file-streaming for very large number of files.
//       2. Write a map processor that actually maps interactions.
//       3. Avoid re-computing Luni

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
