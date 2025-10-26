#ifndef LAHUTA_PIPELINE_DYNAMIC_KEYS_HPP
#define LAHUTA_PIPELINE_DYNAMIC_KEYS_HPP

// clang-format off
namespace lahuta::pipeline {

inline constexpr const char* CTX_SYSTEM_KEY      = "system";             // shared_ptr<const Luni>
inline constexpr const char* CTX_TOPOLOGY_KEY    = "topology";           // shared_ptr<Topology>
inline constexpr const char* CTX_CONFORMER_KEY   = "lahuta.conformer";   // shared_ptr<RDKit::Conformer>
inline constexpr const char* CTX_FRAME_KEY       = "lahuta.frame";       // shared_ptr<FrameMetadata>
inline constexpr const char* CTX_MODEL_METADATA_KEY = "lahuta.model_metadata"; // shared_ptr<ModelMetadata>

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_DYNAMIC_KEYS_HPP
