#ifndef LAHUTA_ANALYSIS_CONTACTS_HOOKS_HPP
#define LAHUTA_ANALYSIS_CONTACTS_HOOKS_HPP

#include <string>

#include "pipeline/task/context.hpp"
#include "serialization/json.hpp"

namespace lahuta::analysis {
namespace P = lahuta::pipeline;

inline std::string build_contacts_summary_json(const P::TaskContext &ctx) {
  const std::string *ok   = ctx.get_text("contacts_success");
  const std::string *cnt  = ctx.get_text("contacts_count");
  const std::string *prov = ctx.get_text("contacts_provider");

  JsonBuilder builder;
  std::size_t n_contacts = 0;

  if (cnt) try {
      n_contacts = static_cast<std::size_t>(std::stoull(*cnt));
    } catch (...) {
    }

  builder.key("ok")
      .value(ok && *ok == "1")
      .key("num_contacts")
      .value(n_contacts)
      .key("provider")
      .value(prov ? *prov : std::string("unknown"));
  return builder.str();
}

} // namespace lahuta::analysis

#endif // LAHUTA_ANALYSIS_CONTACTS_HOOKS_HPP
