#ifndef LAHUTA_ALIGNER_HPP
#define LAHUTA_ALIGNER_HPP

/*#include "processor.hpp"*/
#include "fseek/ops.hpp"
#include "fseek/utils.hpp"
#include "prefilter.hpp"
#include "seq.hpp"
#include "seq_aligner.hpp"
#include "seq_aligner_builder.hpp"
#include <functional>
#include <memory>
#include <spdlog/spdlog.h>

namespace lahuta {

namespace alignment_computers {
static void print_alignment(const SeqData &query, const SeqData &target, const AlignmentResult &AR) {
  std::cout << "Result: \n" << Matcher::results_to_string(AR.ar) << "\n"
  << "Q Alignment: " << AR.query_alignment(query) << "\n"
  << "T Alignment: " << AR.target_alignment(target) << "\n";
}

static void print_result(SeqData &query, SeqData &target, AlignmentResult &ar) {
  std::cout << "Result: " << Matcher::results_to_string(ar.ar) << "\n"
            << "Q data: " << ar.ar.front().qStartPos << " " << ar.ar.front().qEndPos << "\n"
            << "T data: " << ar.ar.front().dbStartPos << " " << ar.ar.front().dbEndPos << "\n"
            << "Q Alignment: " << ar.query_alignment(query) << "\n"
            << "T Alignment: " << ar.target_alignment(target) << "\n";

  auto scores = std::make_shared<AlignmentScores>(query, target, ar.ar[0]);
  std::cout << "RMSD: " << ar.scores.rmsd << "\n"
            << "LDDT: " << ar.scores.avgLddtScore << "\n"
            << "TMscore: " << ar.scores.tmscore << "\n"
            << "QTMscore: " << scores->get_query_score().tmscore << "\n"
            << "TTMscore: " << scores->get_target_score().tmscore << "\n"
            << "Prob: " << ar.scores.prob << "\n";
}
} // namespace alignment_computers

inline void log_sequence_(const SeqData &sd, const std::string &prefix = "") {
  spdlog::info(
      "{} Sequence: {} {} {} residues: {}",
      prefix,
      sd.file_name,
      sd.chain_name,
      sd.size(),
      format_seq(sd.SeqAA));
}

using FileList = std::vector<std::string>;
using AlignProcessFunc = std::function<void(SeqData &query, SeqData &target, AlignmentResult &ar)>;

class Runner {
public:
  virtual ~Runner() = default;
  virtual void run(FileList &query_files, FileList &target_files) = 0;
};

class LahutaAligner : public Runner {

public:
  LahutaAligner(FoldSeekOps &ops, PrefilterOptions &pf_ops) : ops_(ops), pf_ops_(pf_ops) {}

  LahutaAligner(FoldSeekOps *ops = nullptr, PrefilterOptions *pf_ops = nullptr)
      : ops_((ops) ? *ops : FoldSeekOps()), pf_ops_((pf_ops) ? *pf_ops : PrefilterOptions()) {}

  void set_computer(AlignProcessFunc func) { computer = func; }
  AlignProcessFunc get_default_computer() { return alignment_computers::print_alignment; }
  AlignProcessFunc get_computer() { return (computer) ? computer : get_default_computer(); }

  void run(FileList &query_files, FileList &target_files) override {
    read_files(query_files, target_files);
    build_resources();
    process_sequences();
  }

private:
  void read_files(FileList &query_files, FileList &target_files) {

    SeqCollection queries_ = extract_all(ops_, query_files);
    /*SeqCollection targets_ = extract_all(ops_, target_files);*/
    SeqCollection targets_ = extract_all_parallel(ops_, target_files);

    spdlog::warn(
        "Read {} queries and {} targets: {} ",
        queries_.size(),
        targets_.size(),
        targets_.total_length);

    queries = std::make_unique<SeqCollection>(std::move(queries_));
    targets = std::make_unique<SeqCollection>(std::move(targets_));
  }

  // NOTE: the motivation to separate this step is to allow:
  // - avoid buiding the index if not needed (implemented)
  // - to re-use the index for different chunks (not implemented)
  void build_resources() {
    if (pf_ops_.use_prefilter) {
      SeqFilter seq_filter_(*queries, *targets, pf_ops_);
      seq_filter_.build_index();
      seq_filter = std::make_unique<SeqFilter>(std::move(seq_filter_));
    }

    spdlog::warn("Building alignment size info: {} queries and {} targets", queries->size(), targets->size());
    aligner = SeqAlignerBuilder(ops_).build(*queries, *targets);
    aligner->set_needs_lddt(false);

    int num_threads = std::thread::hardware_concurrency();
    thread_aligners.resize(num_threads);
    for (int i = 0; i < num_threads; ++i) {
      thread_aligners[i] = SeqAlignerBuilder(ops_).build(*queries, *targets);
      /*thread_aligners[i]->set_needs_lddt(false);*/
    }

    pool = std::make_unique<ctpl::thread_pool>(num_threads);
  }

  void process_sequences() {
    for (auto &query : *queries) {
      if (spdlog::default_logger()->level() == spdlog::level::info) log_sequence_(query, "Q:");

      if (pf_ops_.use_prefilter) {
        Hits hits = seq_filter->filter(query);
        /*process_input_parallel(query, hits);*/
        process_input(query, hits);
      } else {
        process_input(query, *targets);
      }
    }
  }

  void process_input(SeqData &query, const Hits &hits) {
    for (const int &hit_id : hits) {
      SeqData &target = (*targets)[hit_id];

      if (spdlog::default_logger()->level() == spdlog::level::info) log_sequence_(target, "T:");

      auto alignment_result = aligner->align(query, target);
      if (!alignment_result.success) {
        /*spdlog::info("Alignment unsuccessful with {} - {}", target.file_name, target.chain_name);*/
        continue;
      };

      get_computer()(query, target, alignment_result);
    }
  }

  void process_input_parallel(SeqData &query, const Hits &hits) {
    std::vector<std::future<AlignmentResult>> futures;
    futures.reserve(hits.size());

    auto thread_aligners_ = &thread_aligners;
    for (const auto &hit_id : hits) {
      SeqData &target = (*targets)[hit_id];
      futures.emplace_back(pool->push([this, &query, &target](int thread_id) {
        return thread_aligners[thread_id]->align(query, target);
      }));
    }

    for (size_t i = 0; i < futures.size(); ++i) {
      auto alignment_result = futures[i].get();
      if (!alignment_result.success) {
        continue;
      }

      SeqData &target = (*targets)[hits[i]];
      get_computer()(query, target, alignment_result);
    }
  }

  void _process_input_parallel(SeqData &query, const Hits &hits) {

    struct alig_res {
      SeqData &query;
      SeqData &target;
      AlignmentResult alignment_result;
    };

    int num_threads = std::thread::hardware_concurrency();
    ctpl::thread_pool pool(num_threads);
    std::vector<std::future<alig_res>> futures;

    for (const auto &hit_id : hits) {
      SeqData &target = (*targets)[hit_id];
      futures.emplace_back(pool.push([this, &query, &target](int id) { //
        AlignmentResult alignment_result = aligner->align(query, target);
        return alig_res{query, target, std::move(alignment_result)};
      }));
    }

    for (auto &f : futures) {
      auto result = f.get();
      if (!result.alignment_result.success) {
        continue;
      };

      get_computer()(result.query, result.target, result.alignment_result);
    }
  }

  void process_input(SeqData &query, SeqCollection &targets) {
    for (auto &target : targets) {
      if (spdlog::default_logger()->level() == spdlog::level::info) log_sequence_(target, "T:");

      auto alignment_result = aligner->align(query, target);
      if (!alignment_result.success) {
        spdlog::info("Alignment unsuccessful with {} - {}", target.file_name, target.chain_name);
        continue;
      };

      get_computer()(query, target, alignment_result);
    }
  }

  FoldSeekOps ops_;
  PrefilterOptions pf_ops_;
  std::unique_ptr<SeqFilter> seq_filter;
  std::unique_ptr<SeqAligner> aligner;
  std::vector<std::unique_ptr<SeqAligner>> thread_aligners{1};
  std::unique_ptr<SeqCollection> queries;
  std::unique_ptr<SeqCollection> targets;
  AlignProcessFunc computer;

  std::unique_ptr<ctpl::thread_pool> pool;
};

} // namespace lahuta

#endif // LAHUTA_ALIGNER_HPP
