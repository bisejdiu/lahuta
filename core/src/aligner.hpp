#ifndef LAHUTA_ALIGNER_HPP
#define LAHUTA_ALIGNER_HPP

#include <memory>

#include <spdlog/spdlog.h>

#include "fseek/ops.hpp"
#include "fseek/utils.hpp"
#include "logging/logging.hpp"
#include "prefilter.hpp"
#include "procs/chunking/chunk_defs.hpp"
#include "seq.hpp"
#include "seq_aligner.hpp"
#include "seq_aligner_builder.hpp"

// clang-format off

namespace lahuta {

namespace alignment_computers {

static Matcher::result_t
print_alignment(const SeqData &query, const SeqData &target, const AlignmentResult &AR) {
  std::cout << "Result: \n"
            << Matcher::results_to_string(AR.ar) << "\n"
            << "Q Alignment: " << AR.query_alignment(query) << "\n"
            << "T Alignment: " << AR.target_alignment(target) << "\n";

  return AR.ar.front();
}

static void _print_alignment(const SeqData &query, const SeqData &target, const AlignmentResult &AR) {
  std::cout << "Result: \n"
            << Matcher::results_to_string(AR.ar) << "\n"
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
  Logger::get_logger()->info(
      "{} Sequence: {} {} {} residues: {}",
      prefix,
      sd.file_name,
      sd.chain_name,
      sd.size(),
      format_seq(sd.SeqAA));
}

// FIX: Should we have one Matcher::result_t per query-target pair?
struct AlignerResults {
  std::shared_ptr<SeqData> query;
  std::shared_ptr<SeqData> target;
  std::vector<Matcher::result_t> results;
};

struct AlignerResultsX {
    const lahuta::SeqData* query;          // no ownership
    const lahuta::SeqData* target;         // "
    std::vector<Matcher::result_t> results;
};


class LahutaAlignerBase {
public:
  virtual void run(FileList &query_files, FileList &target_files) = 0;
  virtual ~LahutaAlignerBase() = default;
};


class LahutaAligner: public LahutaAlignerBase {

public:
    LahutaAligner(FoldSeekOps ops = {}, PrefilterOptions pf_ops = {}, unsigned int n_threads = 0) 
      : ops_(std::move(ops)) , pf_ops_(std::move(pf_ops)), n_threads(n_threads) {}

  void run(std::vector<std::string> &query_files, std::vector<std::string> &target_files) override {
    load_sequences(query_files, target_files);
    build_resources();
    /*process_sequences();*/
    process_sequences();
  }

  std::vector<AlignerResults> &get_results() { return results; }

private:
  void load_sequences(std::vector<std::string> &query_files, std::vector<std::string> &target_files) {

    queries = extract_all(ops_, query_files);
    targets = extract_all_parallel(ops_, target_files);

    Logger::get_logger()->warn("Read {} queries and {} targets: {} ", queries.size(), targets.size(), targets.total_length);
  }

  // NOTE: the motivation to separate this step is to:
  // - avoid buiding the index if not needed (implemented)
  // - re-use the index for different chunks (not implemented)
  void build_resources() {

    if (pf_ops_.use_prefilter) {
      SeqFilter seq_filter_(queries, targets, pf_ops_);
      seq_filter_.build_index();
      seq_filter = std::make_unique<SeqFilter>(std::move(seq_filter_));
    }

    uint32_t max_size = static_cast<unsigned int>(std::max(queries.size(), targets.size()));
    uint32_t default_threads = std::min(std::thread::hardware_concurrency(), max_size);
    uint32_t num_threads =  n_threads > 0 ? std::min(n_threads, default_threads) : default_threads;

    thread_aligners.resize(num_threads);
    for (auto &aligner : thread_aligners) {
      aligner = SeqAlignerBuilder(ops_).build(queries, targets);
      aligner->set_needs_lddt(false);
    }

    thread_pool_ = std::make_unique<ctpl::thread_pool>(num_threads);
  }

  void process_sequences() {
    auto start = std::chrono::high_resolution_clock::now();

    for (std::size_t qi = 0; qi < queries.size(); ++qi) {
      if (pf_ops_.use_prefilter) {
        Hits hits = seq_filter->filter(queries[qi]);
        process_alignments(qi, hits);
      }
      else {
        process_alignments(qi, targets);
      }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto dur = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    Logger::get_logger()->info("Ali Processing time: {} us", dur);
  }

  template <typename TargetType>
  void process_alignments(std::size_t query_idx, TargetType &targets) {
    SeqData &query = get_target(targets, query_idx);

    std::vector<std::future<AlignmentResult>> futures;
    futures.reserve(targets.size() - query_idx - 1);

    std::vector<std::size_t> target_idxs;
    target_idxs.reserve(futures.capacity());

    for (std::size_t j = query_idx + 1; j < targets.size(); ++j) {
      SeqData &target = get_target(targets, j);
      futures.emplace_back(
        thread_pool_->push([this, &query, &target](int tid) {
            return thread_aligners[tid]->align(query, target);
          }
        )
      );
      target_idxs.push_back(j);
    }

    for (std::size_t k = 0; k < futures.size(); ++k) {
      auto ar = futures[k].get();
      if (!ar.success) continue;

      SeqData &tgt = get_target(targets, target_idxs[k]);
      results.push_back({std::make_shared<SeqData>(query), std::make_shared<SeqData>(tgt), ar.ar});
    }
  }

  SeqData &get_target(SeqCollection &targets, size_t index) const { return targets[index]; }
  SeqData &get_target(Hits &hits, size_t index) { return targets[hits[index]]; }

private:
  unsigned int n_threads;

  FoldSeekOps ops_;
  PrefilterOptions pf_ops_;

  SeqCollection queries;
  SeqCollection targets;

  std::unique_ptr<SeqFilter> seq_filter;
  std::vector<std::unique_ptr<SeqAligner>> thread_aligners{1};
  std::unique_ptr<ctpl::thread_pool> thread_pool_;

  std::vector<AlignerResults> results;
};

} // namespace lahuta

#endif // LAHUTA_ALIGNER_HPP
