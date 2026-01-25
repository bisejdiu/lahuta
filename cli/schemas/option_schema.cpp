#include <utility>

#include "schemas/option_schema.hpp"

namespace lahuta::cli {

void OptionSchema::add(OptionDef def) { defs_.push_back(std::move(def)); }

void OptionSchema::append(const OptionSchema &other) {
  defs_.insert(defs_.end(), other.defs_.begin(), other.defs_.end());
}

const std::vector<OptionDef> &OptionSchema::defs() const noexcept { return defs_; }

std::vector<option::Descriptor> OptionSchema::build_descriptors() const {
  std::vector<option::Descriptor> descriptors;
  descriptors.reserve(defs_.size() + 1);
  for (const auto &def : defs_) {
    const char *short_opt = def.short_name.empty() ? "" : def.short_name.c_str();
    const char *long_opt  = def.long_name.empty() ? "" : def.long_name.c_str();
    descriptors.push_back({def.index, 0, short_opt, long_opt, def.validator, def.help.c_str()});
  }
  descriptors.push_back({0, 0, 0, 0, 0, 0});
  return descriptors;
}

} // namespace lahuta::cli
