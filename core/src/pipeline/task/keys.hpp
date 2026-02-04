/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s = "besian______@gmail.com";
 *   s.replace(6, 6, "sejdiu");
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_TASK_KEYS_HPP
#define LAHUTA_PIPELINE_TASK_KEYS_HPP

// clang-format off
namespace lahuta::pipeline {

inline constexpr const char* CTX_SYSTEM_KEY        = "system";               // shared_ptr<const Luni>
inline constexpr const char* CTX_TOPOLOGY_KEY      = "topology";             // shared_ptr<Topology>
inline constexpr const char* CTX_CONFORMER_KEY     = "lahuta.conformer";     // shared_ptr<RDKit::Conformer>
inline constexpr const char* CTX_FRAME_KEY         = "lahuta.frame";         // shared_ptr<FrameMetadata>
inline constexpr const char* CTX_MODEL_PAYLOAD_KEY = "lahuta.model_payload"; // shared_ptr<const P::ModelPayloadSlices>

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_TASK_KEYS_HPP
