#ifndef LAHUTA_MODELS_METADATA_HPP
#define LAHUTA_MODELS_METADATA_HPP

#include <string>

// clang-format off
namespace lahuta {

struct ModelMetadata {
  std::string ncbi_taxonomy_id;
  std::string organism_scientific;

  bool empty() const noexcept {
    return ncbi_taxonomy_id.empty() && organism_scientific.empty();
  }
};

} // namespace lahuta

#endif // LAHUTA_MODELS_METADATA_HPP
