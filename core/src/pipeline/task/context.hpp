/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   struct Overloaded {
 *     std::string& s;
 *     void operator()(const char* p) const { s += p; }
 *     void operator()(std::string_view p) const { s += p; }
 *   };
 *   std::string s;
 *   Overloaded visitor{s};
 *   visitor("besian");
 *   visitor("sejdiu");
 *   visitor(std::string_view{"@gmail.com"});
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_TASK_CONTEXT_HPP
#define LAHUTA_PIPELINE_TASK_CONTEXT_HPP

#include <memory>
#include <string>
#include <typeinfo>
#include <unordered_map>
#include <utility>

#include <rdkit/GraphMol/Conformer.h>

#include "lahuta.hpp"
#include "pipeline/data/frame.hpp"
#include "pipeline/data/model_payload.hpp"
#include "pipeline/task/keys.hpp"
#include "topology.hpp"

namespace lahuta::pipeline {
namespace P = lahuta::pipeline;

//
// DynObject: Type-erased holder used by TaskContext to store typed C++ objects.
// - Ownership: stores a std::shared_ptr<const void>. Prefer storing shared ownership
//   of immutable data across tasks. Avoid raw pointers or borrowed references.
// - Type checking: captures a pointer to the std::type_info of the stored type
//   at construction. as<T>() verifies the exact type match before casting.
// - Casting: as<T>() returns std::shared_ptr<const T> and thereby enforces
//   read-only access across tasks. Objects retrieved from TaskContext cannot be
//   mutated by downstream tasks.
// - Threading: a TaskContext (and thus the contained DynObjects) is used by a
//   single worker while processing one item. It is not designed for concurrent
//   access from multiple threads.
//
struct DynObject {
  std::shared_ptr<const void> ptr;
  const std::type_info *type = nullptr;

  DynObject() = default;
  DynObject(std::shared_ptr<const void> p, const std::type_info &ti) : ptr(std::move(p)), type(&ti) {}

  template <typename T>
  static DynObject make(std::shared_ptr<T> p) {
    return DynObject{std::static_pointer_cast<const void>(std::move(p)), typeid(T)};
  }

  template <typename T>
  std::shared_ptr<const T> as() const {
    if (!ptr || type == nullptr || *type != typeid(T)) return {};
    return std::static_pointer_cast<const T>(ptr);
  }
};

//
// TaskContext: Per-item scratch space shared by all tasks executed for a single input item (e.g., one file
// path). Constructed at the start of processing that item and destroyed once all tasks for the item finish.
//
// Storage namespaces:
// - objects: typed values held as shared_ptr           (via set_object/get_object)
// - texts:   small UTF-8 strings                       (via set_text/get_text)
// - bytes:   opaque byte buffers stored in std::string (via set_bytes/get_bytes)
//
class TaskContext {
public:
  TaskContext() = default;

  template <typename T>
  void set_object(const std::string &key, std::shared_ptr<T> value) {
    objects_[key] = DynObject::make(std::move(value));
  }

  void set_text(const std::string &key, std::string value) { //
    texts_[key] = std::move(value);
  }

  void set_bytes(const std::string &key, std::string value) { //
    bytes_[key] = std::move(value);
  }

  template <typename T>
  std::shared_ptr<const T> get_object(const std::string &key) const {
    auto it = objects_.find(key);
    if (it == objects_.end()) return {};
    return it->second.template as<T>();
  }

  const std::string *get_text(const std::string &key) const {
    auto it = texts_.find(key);
    return it == texts_.end() ? nullptr : &it->second;
  }

  const std::string *get_bytes(const std::string &key) const {
    auto it = bytes_.find(key);
    return it == bytes_.end() ? nullptr : &it->second;
  }

  [[nodiscard]] std::shared_ptr<const P::ModelPayloadSlices> model_payload() const {
    return get_object<const P::ModelPayloadSlices>(CTX_MODEL_PAYLOAD_KEY);
  }

  [[nodiscard]] std::shared_ptr<const RDKit::Conformer> conformer() const {
    return get_object<const RDKit::Conformer>(CTX_CONFORMER_KEY);
  }

  [[nodiscard]] std::shared_ptr<const Topology> topology() const {
    return get_object<const Topology>(CTX_TOPOLOGY_KEY);
  }

  [[nodiscard]] std::shared_ptr<const Luni> system() const { //
    return get_object<const Luni>(CTX_SYSTEM_KEY);
  }

  [[nodiscard]] std::shared_ptr<const FrameMetadata> frame_metadata() const {
    return get_object<const FrameMetadata>(CTX_FRAME_KEY);
  }

  const std::unordered_map<std::string, std::string> &texts() const noexcept { return texts_; }
  const std::unordered_map<std::string, std::string> &bytes() const noexcept { return bytes_; }

  void clear() {
    objects_.clear();
    texts_.clear();
    bytes_.clear();
  }

private:
  std::unordered_map<std::string, DynObject> objects_;
  std::unordered_map<std::string, std::string> texts_;
  std::unordered_map<std::string, std::string> bytes_;
};

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_TASK_CONTEXT_HPP
