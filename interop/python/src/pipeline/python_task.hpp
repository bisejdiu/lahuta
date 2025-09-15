#ifndef LAHUTA_BINDINGS_PYTHON_TASK_HPP
#define LAHUTA_BINDINGS_PYTHON_TASK_HPP

#include <mutex>
#include <optional>
#include <string>
#include <utility>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pythoncapi_compat.h" // Backports public C-API to old Pythons

#include "pipeline/context.hpp"
#include "pipeline/dynamic/types.hpp"

// clang-format off
namespace py = pybind11;
namespace lahuta::bindings {
using namespace lahuta::pipeline::dynamic;

class PyCallableTask : public ITask {
public:
  // Callable contract: fn(ctx: PipelineContext) -> Any // not quite "Any"
  // Return value handling:
  // - None   : no emission
  // - str    : emit as text
  // - other  : serialize with json.dumps(..., ensure_ascii=False) and emit as text
  PyCallableTask(py::object fn, std::string store_key, std::optional<std::string> emit_channel, bool store, bool serialize)
    : func_(std::move(fn)), store_key_(std::move(store_key)), emit_channel_(std::move(emit_channel)),
      store_(store), serialize_(serialize) {

    if (!func_ || func_.is_none() || !PyCallable_Check(func_.ptr()))
      throw std::invalid_argument("PyCallableTask: fn must be callable");

    // Preload a fast JSON dumper once
    // try orjson before falling back to stdlib json
    try {
      py::gil_scoped_acquire gil;
      try {
        auto mod = py::module_::import("orjson");
        dumps_ = mod.attr("dumps");
        dumps_returns_bytes_ = true;
      } catch (...) {
        auto mod = py::module_::import("json");
        dumps_ = mod.attr("dumps");
        dumps_returns_bytes_ = false;
      }
    } catch (...) {
      // leave dumps_ empty. We'll lazy-import in run() if needed
    }
  }

  ~PyCallableTask() override {
    if (!Py_IsInitialized() || Py_IsFinalizing()) {
      try { (void) func_.release().ptr(); } catch (...) {}
      try { (void)dumps_.release().ptr(); } catch (...) {}
      return;
    }
    try {
      py::gil_scoped_acquire gil;
      func_  = py::object();
      dumps_ = py::object();
    } catch (...) {
      try { (void) func_.release().ptr(); } catch (...) {}
      try { (void)dumps_.release().ptr(); } catch (...) {}
    }
  }

  TaskResult run(const std::string& item_path, TaskContext& ctx) override {
    //
    // Optional per-task serialization: when enabled, only one thread executes
    // this Python callable at a time, preserving item-level parallelism without
    // collapsing the entire pipeline. When disabled, concurrency is limited only
    // by the GIL and user code behavior.
    //
    std::unique_lock<std::mutex> lk(mu_, std::defer_lock);
    if (serialize_) lk.lock();
    py::gil_scoped_acquire gil;
    try {
      //
      // We always construct and pass PipelineContext to the Python callable. This simplifies
      // things for us and for the user!  - Besian, September 2025
      //
      py::object arg = py::cast(PyPipelineContext(&ctx, item_path));

      py::object resp = func_(arg);
      if (resp.is_none()) return TaskResult{true, {}}; // skip emission

      //
      // Return handling: if str, emit as text, otherwise json serialize once via cached dumps
      // We should consider wrapping the str in a JSON object, which would make things more uniform, 
      // and simplify type annotations on the Python side. - Besian, September 2025
      //
      std::string payload;
      if (py::isinstance<py::str>(resp)) {
        payload = resp.cast<std::string>();
      } else {
        if (!dumps_) {
          // extremely rare - constructor failed - lazy import
          auto mod = py::module_::import("json");
          dumps_ = mod.attr("dumps");
          dumps_returns_bytes_ = false;
        }
        if (dumps_returns_bytes_) {
          py::bytes b = dumps_(resp);
          payload = b.cast<std::string>();
        } else {
          payload = dumps_(resp, py::arg("ensure_ascii") = false).cast<std::string>();
        }
      }
      if (store_) ctx.set_text(store_key_, payload);

      TaskResult tr; tr.ok = true;
      if (emit_channel_.has_value()) tr.emits.push_back(Emission{*emit_channel_, std::move(payload)});
      return tr;
    } catch (const py::error_already_set& ex) {
      // Extract a concise "Type: message" without traceback
      std::string type_name = "PythonError";
      std::string msg;
      try {
        py::gil_scoped_acquire gil2;
        if (ex.type()) {
          try { type_name = py::str(py::handle(ex.type()).attr("__name__")).cast<std::string>(); } catch (...) {}
        }
        if (ex.value()) {
          try { msg = py::str(py::handle(ex.value())).cast<std::string>(); } catch (...) {}
        }
      } catch (...) {
        // falling back to what() if we cannot introspect
        msg = ex.what() ? std::string(ex.what()) : std::string("error");
      }
      if (msg.empty()) {
        // Fallback: try what() and trim
        msg = ex.what() ? std::string(ex.what()) : std::string("error");
      }
      // Trim at first newline to remove traceback or extra context
      if (auto p = msg.find('\n'); p != std::string::npos) msg.erase(p);
      std::string short_msg = type_name + ": " + msg;
      return emit_error(item_path, ctx, "python", short_msg.c_str());
    } catch (const std::exception& ex) {
      // Keep only first line of std::exception
      std::string m = ex.what() ? std::string(ex.what()) : std::string("error");
      if (auto p = m.find('\n'); p != std::string::npos) m.erase(p);
      return emit_error(item_path, ctx, "cpp", m.c_str());
    } catch (...) {
      return emit_error(item_path, ctx, "unknown", "unknown error");
    }
  }

private:
  static std::string json_escape(std::string s) {
    static const char* hex = "0123456789ABCDEF";

    std::string out; out.reserve(s.size() + 8);
    for (char c : s) {
      unsigned char uc = static_cast<unsigned char>(c);
      switch (c) {
        case '\\': out += "\\\\"; break;
        case '"':  out += "\\\""; break;
        case '\n': out += "\\n";  break;
        case '\r': out += "\\r";  break;
        case '\t': out += "\\t";  break;
        default:
          if (uc < 0x20) {
            char buf[] = "\\u00XX";
            buf[4] = hex[(uc >> 4) & 0xF];
            buf[5] = hex[uc & 0xF];
            out.append(buf, 6);
          } else {
            out += c;
          }
          break;
      }
    }
    return out;
  }

  TaskResult emit_error(const std::string& item_path, TaskContext& ctx, const char* source, const char* message) {
    // Build a small JSON object with error info. Keep dependencies minimal.
    std::string src = source ? source : "unknown";
    std::string msg = message ? message : "";
    std::string payload = std::string("{\"error\":{\"source\":\"") + json_escape(src) +
                          "\",\"message\":\"" + json_escape(msg) +
                          "\"},\"path\":\"" + json_escape(item_path) + "\"}";
    if (store_) ctx.set_text(store_key_, payload);
    TaskResult tr; tr.ok = true;
    if (emit_channel_.has_value()) tr.emits.push_back(Emission{*emit_channel_, std::move(payload)});
    return tr;
  }

  py::object func_;
  py::object dumps_;
  bool dumps_returns_bytes_ = false;
  std::string store_key_;
  std::optional<std::string> emit_channel_;
  bool store_;
  bool serialize_ = false;
  mutable std::mutex mu_;
};

} // namespace lahuta::bindings

#endif // LAHUTA_BINDINGS_PYTHON_TASK_HPP
