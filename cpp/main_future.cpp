#include "mapper.hpp"
#include "processor.hpp"
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
            << target.file_name << "_" << target.chain_name << " : " 
            << Matcher::results_to_string(ar.ar) << std::endl;
  /*mapper.print();*/
  /*alignment_computers::print_result(query, target, ar);*/
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
  spdlog::level::level_enum log_level = spdlog::level::warn;
  spdlog::set_level(log_level);

  DirectoryHandler dir_handler("/Users/bsejdiu/data/PDB_ARCHIVE_UNCOMP", ".cif", true);
  /*DirectoryHandler dir_handler("/Users/bsejdiu/data/mini_pdb_uncomp", ".cif", true);*/

  auto pdb_targets = dir_handler.get_all_files();

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
        "/Users/bsejdiu/data/gpcrs/small/b/2r4s.cif",
        "/Users/bsejdiu/data/gpcrs/small/b/3nya.cif",
        "/Users/bsejdiu/data/gpcrs/small/b/8jj8.cif"};

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

    FoldSeekOps ops;
    PrefilterOptions pf_ops;
    pf_ops.use_prefilter = true;

    std::unique_ptr<LahutaAligner> aligner = std::make_unique<LahutaAligner>(ops, pf_ops);
    aligner->set_computer(mapping_processor);

    LahutaProcessor processor(config); // , ops, pf_ops);
    processor.set_runner(std::move(aligner));
    /*processor.process_files(query_files);*/
    /*processor.process_files(query_files, target_files);*/
    /*processor.process_files({"/Users/bsejdiu/projects/lahuta/cpp/build/4ami.cif"}, pdb_targets);*/


    processor.process_files(query_files, pdb_targets);
    /*processor.process_files(query_files_t, target_files_t);*/

  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
  }
}
