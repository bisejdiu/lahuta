#ifndef LAHUTA_SOURCES_REALIZER_HPP
#define LAHUTA_SOURCES_REALIZER_HPP

#include <optional>
#include <string>
#include <unordered_set>
#include <variant>

#include "sources/lmdb.hpp"
#include "sources/nmr.hpp"
#include "sources/trajectory.hpp"

// clang-format off
namespace lahuta::sources {

template <class... Ts>
struct Overloaded : Ts... { using Ts::operator()...; };
template <class... Ts>
Overloaded(Ts...)->Overloaded<Ts...>;

//
// DefaultRealizer expands a single in-flight IngestDescriptor into one or more
// PipelineItems. Currently ItemStream feeds descriptors serially, so we just
// need per-descriptor "emit once" flags. If we later want to to handle
// multiple descriptors concurrently, this must be replaced with a per-id
// state map, and correct synchronization. - Besian, October 2025
//
class Realizer {
public:
  void reset() {
    file_emitted_.clear();
    lmdb_.reset();
    nmr_ .reset();
    traj_.reset();
  }

  void reset(const std::string& id) {
    file_emitted_.erase(id);
    lmdb_.reset(id);
    nmr_ .reset(id);
    traj_.reset(id);
  }

  std::optional<PipelineItem> next(const IngestDescriptor& desc) {
    return std::visit(Overloaded{
      [&](const FileRef& ref) { return emit_file(desc, ref); },
      [&](const LMDBRef&)     { return lmdb_ .next(desc); },
      [&](const NMRRef&)      { return nmr_  .next(desc); },
      [&](const MDRef&)       { return traj_ .next(desc); }
    }, desc.origin);
  }

private:
  std::optional<PipelineItem> emit_file(const IngestDescriptor& desc, const FileRef& ref) {
    auto [it, inserted] = file_emitted_.insert(desc.id);
    if (!inserted) return std::nullopt;
    PipelineItem item;
    item.session_id = desc.id;
    item.item_path  = ref.path;
    return item;
  }

  std::unordered_set<std::string> file_emitted_;
  LMDBRealizer lmdb_;
  NMRRealizer  nmr_;
  TrajectoryRealizer traj_;
};

} // namespace lahuta::sources

#endif // LAHUTA_SOURCES_REALIZER_HPP
