#ifndef LAHUTA_ALIGNER_P_HPP
#define LAHUTA_ALIGNER_P_HPP

#include "logging.hpp"
#include "prefilter.hpp"
#include "seq.hpp"
#include "seq_aligner.hpp"
#include "seq_aligner_builder.hpp"
#include "aligner.hpp" // for AlignerResults
#include <memory>
#include <spdlog/spdlog.h>
#include <unordered_set>
#include <string>
#include <utility>

namespace lahuta {

/**
 * Pipeline-friendly version of LahutaAligner that works with pre-loaded sequences
 */
class LahutaAlignerP {
public:
    /**
     * Construct a LahutaAlignerP
     * 
     * @param ops FoldSeek operation parameters
     * @param pf_ops Prefilter options
     * @param n_threads Number of threads to use
     */
    LahutaAlignerP(FoldSeekOps ops = {}, PrefilterOptions pf_ops = {}, unsigned int n_threads = 0) 
        : ops_(std::move(ops)), pf_ops_(std::move(pf_ops)), n_threads(n_threads) {}

    /**
     * Run alignments on the provided sequences
     * 
     * @param queries Query sequences
     * @param targets Target sequences
     */
    void run(const SeqCollection& queries, const SeqCollection& targets) {
        // Clear previous results
        results_.clear();
        
        // Check if queries and targets are the same collection
        bool same_collection = &queries == &targets;
        
        Logger::get_instance().log(Logger::LogLevel::Info, 
            "Aligning {} queries against {} targets ({})", 
            queries.size(), targets.size(),
            same_collection ? "same collection" : "different collections");
        
        // Store references to the input sequences
        queries_ = queries;
        targets_ = targets;
        
        // Build resources
        build_resources();
        
        // Process sequences with pairwise alignments
        process_sequence_pairs(same_collection);
        
        Logger::get_instance().log(Logger::LogLevel::Info, 
            "Completed {} alignments", results_.size());
    }

    /**
     * Get the alignment results
     * 
     * @return Vector of alignment results
     */
    std::vector<AlignerResults>& get_results() { 
        return results_;
    }

private:
    /**
     * Build resources for sequence alignment
     */
    void build_resources() {
        if (pf_ops_.use_prefilter) {
            SeqFilter seq_filter_(queries_, targets_, pf_ops_);
            seq_filter_.build_index();
            seq_filter = std::make_unique<SeqFilter>(std::move(seq_filter_));
        }

        uint32_t max_size = static_cast<unsigned int>(std::max(queries_.size(), targets_.size()));
        uint32_t default_threads = std::min(std::thread::hardware_concurrency(), max_size);
        uint32_t num_threads = n_threads > 0 ? std::min(n_threads, default_threads) : default_threads;

        thread_aligners.resize(num_threads);
        for (auto &aligner : thread_aligners) {
            aligner = SeqAlignerBuilder(ops_).build(queries_, targets_);
            aligner->set_needs_lddt(false);
        }

        thread_pool_ = std::make_unique<ctpl::thread_pool>(num_threads);
    }

    /**
     * Process sequences to generate pairwise alignments
     * 
     * @param same_collection Whether queries and targets are the same collection
     */
    void process_sequence_pairs(bool same_collection) {
        auto start = std::chrono::high_resolution_clock::now();
        
        // Create all unique query-target pairs to align
        std::vector<std::pair<size_t, size_t>> pairs_to_align;
        std::unordered_set<std::string> pair_keys;
        
        // Generate all required alignments
        for (size_t i = 0; i < queries_.size(); ++i) {
            const auto& query = queries_[i];
            
            // Determine which targets to align with
            size_t start_j = 0;
            if (same_collection) {
                // If same collection, only align with sequences after this one
                // to avoid duplicate alignments (A-B and B-A)
                start_j = i + 1;
            }
            
            for (size_t j = start_j; j < targets_.size(); ++j) {
                // Skip self-comparison
                if (same_collection && i == j) {
                    continue;
                }
                
                const auto& target = targets_[j];
                
                // Create a unique key for this pair
                std::string pair_key = create_pair_key(
                    query.file_name, query.chain_name,
                    target.file_name, target.chain_name);
                
                // Skip if we've already seen this pair
                if (pair_keys.find(pair_key) != pair_keys.end()) {
                    continue;
                }
                
                // Add to our tracking structures
                pair_keys.insert(pair_key);
                pairs_to_align.emplace_back(i, j);
            }
        }
        
        Logger::get_instance().log(Logger::LogLevel::Info, 
            "Identified {} unique alignment pairs", pairs_to_align.size());
        
        // Execute all alignments in parallel
        std::vector<std::future<AlignmentResult>> futures;
        futures.reserve(pairs_to_align.size());
        
        for (const auto& [query_idx, target_idx] : pairs_to_align) {
            auto& query = queries_[query_idx];
            auto& target = targets_[target_idx];
            
            futures.emplace_back(
                thread_pool_->push([this, &query, &target](int tid) {
                    return thread_aligners[tid]->align(query, target);
                })
            );
        }
        
        // Process results
        for (size_t k = 0; k < futures.size(); ++k) {
            auto ar = futures[k].get();
            if (!ar.success) continue;
            
            const auto& [query_idx, target_idx] = pairs_to_align[k];
            const auto& query = queries_[query_idx];
            const auto& target = targets_[target_idx];
            
            Logger::get_instance().log(Logger::LogLevel::Info, 
                "Aligned {}::{} - {}::{}", 
                query.file_name, query.chain_name,
                target.file_name, target.chain_name);
            
            results_.push_back({
                std::make_shared<SeqData>(query), 
                std::make_shared<SeqData>(target), 
                ar.ar
            });
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        auto dur = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        Logger::get_instance().log(Logger::LogLevel::Info, 
            "Alignment processing time: {} us", dur);
    }

    /**
     * Create a unique key for a pair of chains to track which ones have been aligned
     * This properly handles the A-B vs B-A equivalence
     * 
     * @param file1 First file name
     * @param chain1 First chain name
     * @param file2 Second file name
     * @param chain2 Second chain name
     * @return A unique string key for this pair
     */
    std::string create_pair_key(const std::string& file1, const std::string& chain1,
                             const std::string& file2, const std::string& chain2) {
        // Always sort the pair to ensure consistent ordering
        if (file1 < file2 || (file1 == file2 && chain1 < chain2)) {
            return file1 + "::" + chain1 + "-" + file2 + "::" + chain2;
        } else {
            return file2 + "::" + chain2 + "-" + file1 + "::" + chain1;
        }
    }

private:
    unsigned int n_threads;

    FoldSeekOps ops_;
    PrefilterOptions pf_ops_;

    SeqCollection queries_;
    SeqCollection targets_;

    std::unique_ptr<SeqFilter> seq_filter;
    std::vector<std::unique_ptr<SeqAligner>> thread_aligners{1};
    std::unique_ptr<ctpl::thread_pool> thread_pool_;

    std::vector<AlignerResults> results_;
};

} // namespace lahuta

#endif // LAHUTA_ALIGNER_P_HPP 
