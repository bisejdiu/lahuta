#ifndef LAHUTA_CLI_TASKS_CONTACTS_MD_SOURCE_HPP
#define LAHUTA_CLI_TASKS_CONTACTS_MD_SOURCE_HPP

#include <memory>
#include <string>
#include <vector>

#include "sources/descriptor.hpp"

namespace lahuta::cli::contacts {

struct MdInputs {
  std::string structure_path;
  std::vector<std::string> trajectory_paths;
};

[[nodiscard]] MdInputs parse_md_inputs(const std::vector<std::string> &md_files);

[[nodiscard]] std::shared_ptr<lahuta::sources::IDescriptor>
make_md_source_descriptor(std::string structure_path, std::vector<std::string> trajectory_paths);

} // namespace lahuta::cli::contacts

#endif // LAHUTA_CLI_TASKS_CONTACTS_MD_SOURCE_HPP
