#ifndef LAHUTA_PIPELINE_DYNAMIC_TYPES_HPP
#define LAHUTA_PIPELINE_DYNAMIC_TYPES_HPP

#include <memory>
#include <string>
#include <typeinfo>
#include <unordered_map>
#include <utility>
#include <vector>

// clang-format off
namespace lahuta::pipeline::dynamic {

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
  const std::type_info* type = nullptr;

  DynObject() = default;
  DynObject(std::shared_ptr<const void> p, const std::type_info& ti) : ptr(std::move(p)), type(&ti) {}

  template<typename T>
  static DynObject make(std::shared_ptr<T> p) {
    return DynObject{std::static_pointer_cast<const void>(std::move(p)), typeid(T)};
  }

  template<typename T>
  std::shared_ptr<const T> as() const {
    if (!ptr || type == nullptr || *type != typeid(T)) return {};
    return std::static_pointer_cast<const T>(ptr);
  }
};

//
// TaskContext: Per-item scratch space shared by all tasks executed for a single input item (e.g., one file path).
// Constructed at the start of processing that item and destroyed once all tasks for the item finish.
//
// Storage namespaces:
// - objects: typed values held as shared_ptr           (via set_object/get_object)
// - texts:   small UTF-8 strings                       (via set_text/get_text)
// - bytes:   opaque byte buffers stored in std::string (via set_bytes/get_bytes)
//
class TaskContext {
public:
  TaskContext() = default;

  template<typename T>
  void set_object(const std::string& key, std::shared_ptr<T> value) {
    objects_[key] = DynObject::make(std::move(value));
  }

  void set_text(const std::string& key, std::string value) {
    texts_[key] = std::move(value);
  }

  void set_bytes(const std::string& key, std::string value) {
    bytes_[key] = std::move(value);
  }

  template<typename T>
  std::shared_ptr<const T> get_object(const std::string& key) const {
    auto it = objects_.find(key);
    if (it == objects_.end()) return {};
    return it->second.template as<T>();
  }

  const std::string* get_text(const std::string& key) const {
    auto it = texts_.find(key);
    return it == texts_.end() ? nullptr : &it->second;
  }

  const std::string* get_bytes(const std::string& key) const {
    auto it = bytes_.find(key);
    return it == bytes_.end() ? nullptr : &it->second;
  }

  const std::unordered_map<std::string, std::string>& texts() const noexcept { return texts_; }
  const std::unordered_map<std::string, std::string>& bytes() const noexcept { return bytes_; }

  void clear() {
    objects_.clear();
    texts_  .clear();
    bytes_  .clear();
  }

private:
  std::unordered_map<std::string, DynObject> objects_;
  std::unordered_map<std::string, std::string> texts_;
  std::unordered_map<std::string, std::string> bytes_;
};

//
// Emission: One unit of output produced by a task.
// - channel: logical topic used by the ChannelMultiplexer to route to sinks.
//            Multiple tasks may emit to the same channel. Multiple sinks may subscribe to a channel (many-to-many).
// - payload: textual data (plain text, json). Tasks decide the serialization.
// Ordering & delivery:
// - Within an item, emissions are produced in the order tasks run.
//   Across items, emissions may interleave when running with multiple threads.
//
struct Emission {
  std::string channel;  // logical name, e.g. "contacts", "system"
  std::string payload;  // text/json
};

using EmissionList = std::vector<dynamic::Emission>;

//
// Lightweight view used internally by the backpressure/writer threads.
// Producers do not construct this. ChannelMultiplexer creates EmissionView from
// shared backing buffers and interned channel ids.
//
struct EmissionView {
  uint32_t         channel_id;   // interned id for routing/metrics
  std::string_view payload;      // view into a shared backing string
};

//
// TaskResult: Outcome of running a task on a single item.
// - ok=false: stop executing downstream tasks for that item. Emissions already
//             produced for the item are still forwarded to sinks. We don't do rollbacks.
// - emits:    zero or more Emission records to dispatch via the multiplexer.
//
// Implementations should avoid throwing. Convert recoverable failures to ok=false
// so the engine can continue with other items.
//
struct TaskResult {
  bool ok = true;
  std::vector<Emission> emits;
};

//
// ITask: Type-erased task interface implemented by built-in adapters and bindings.
// - Contract: run(item_path, ctx) returns a TaskResult. May read/write the TaskContext and emit payloads.
// - Thread-safety: when registering tasks, the StageManager records a
//   thread_safe hint. If any task is unsafe, the whole run is serialized.
// - Error handling: prefer returning ok=false over throwing exceptions across worker boundaries.
//
class ITask {
public:
  virtual ~ITask() = default;
  virtual TaskResult run(const std::string& item_path, TaskContext& ctx) = 0;
};

} // namespace lahuta::pipeline::dynamic

#endif // LAHUTA_PIPELINE_DYNAMIC_TYPES_HPP
