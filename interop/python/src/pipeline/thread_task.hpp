/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto append_if_string = [](std::string& s, auto&& arg)
 *       -> std::enable_if_t<std::is_convertible_v<decltype(arg), std::string_view>> {
 *     s += arg;
 *   };
 *   std::string s;
 *   append_if_string(s, "besian");
 *   append_if_string(s, "sejdiu");
 *   append_if_string(s, "@gmail.com");
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_BINDINGS_PYTHON_TASK_HPP
#define LAHUTA_BINDINGS_PYTHON_TASK_HPP

#include <mutex>
#include <optional>
#include <string>
#include <utility>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pipeline/context.hpp"
#include "pipeline/task/api.hpp"
#include "pythoncapi_compat.h" // Backports public C-API to old Pythons

namespace py = pybind11;
namespace lahuta::bindings {
namespace P = lahuta::pipeline;

class PyCallableTask : public P::ITask {
public:
  // Callable contract: fn(ctx: PipelineContext) -> Any // not quite "Any"
  // Return value handling:
  // - None   : no emission
  // - str    : emit as text
  // - other  : serialize with orjson.dumps and emit as text
  PyCallableTask(py::object fn, std::string store_key, std::optional<std::string> emit_channel, bool store,
                 bool serialize)
      : func_(std::move(fn)), store_key_(std::move(store_key)), emit_channel_(std::move(emit_channel)),
        store_(store), serialize_(serialize) {

    if (!func_ || func_.is_none() || !PyCallable_Check(func_.ptr())) {
      throw std::invalid_argument("PyCallableTask: fn must be callable");
    }

    py::gil_scoped_acquire gil;
    auto mod             = py::module_::import("orjson");
    dumps_               = mod.attr("dumps");
    opt_serialize_numpy_ = mod.attr("OPT_SERIALIZE_NUMPY"); // support for numpy serialization
  }

  // clang-format off
  ~PyCallableTask() override {
    if (!Py_IsInitialized() || Py_IsFinalizing()) {
      try { (void) func_.release().ptr(); } catch (...) {}
      try { (void)dumps_.release().ptr(); } catch (...) {}
      try { (void)opt_serialize_numpy_.release().ptr(); } catch (...) {}
      return;
    }
    try {
      py::gil_scoped_acquire gil;
      func_  = py::object();
      dumps_ = py::object();
      opt_serialize_numpy_ = py::object();
    } catch (...) {
      try { (void) func_.release().ptr(); } catch (...) {}
      try { (void)dumps_.release().ptr(); } catch (...) {}
      try { (void)opt_serialize_numpy_.release().ptr(); } catch (...) {}
    }
  }
  // clang-format on

  P::TaskResult run(const std::string &item_path, P::TaskContext &ctx) override {
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
      if (resp.is_none()) return P::TaskResult{true, {}}; // skip emission

      //
      // Return handling: if str, emit as text, otherwise json serialize once via cached dumps
      // We should consider wrapping the str in a JSON object, which would make things more uniform,
      // and simplify type annotations on the Python side. - Besian, September 2025
      //
      std::string payload;
      if (py::isinstance<py::str>(resp)) {
        payload = resp.cast<std::string>();
        if (store_) ctx.set_text(store_key_, payload);
      } else if (PyBytes_Check(resp.ptr()) || PyByteArray_Check(resp.ptr()) ||
                 PyMemoryView_Check(resp.ptr())) {
        // Treat bytes-like return values as binary payloads
        py::bytes b = py::bytes(resp);
        payload     = b.cast<std::string>();
        if (store_) ctx.set_bytes(store_key_, payload);
      } else {
        //
        // orjson.dumps returns bytes for JSON-serializable objects
        // Using OPT_SERIALIZE_NUMPY to automatically convert numpy arrays to lists
        // Note: orjson.dumps(obj, /, default=None, option=None)
        // We need to pass None for default and our option value
        //
        py::bytes b = dumps_(resp, py::none(), opt_serialize_numpy_);
        payload     = b.cast<std::string>();
        if (store_) ctx.set_text(store_key_, payload);
      }

      P::TaskResult tr;
      tr.ok = true;
      if (emit_channel_.has_value()) tr.emits.push_back(P::Emission{*emit_channel_, std::move(payload)});
      return tr;
    } catch (const py::error_already_set &ex) {
      // Extract a concise "Type: message" without traceback
      std::string type_name = "PythonError";
      std::string msg;
      try {
        py::gil_scoped_acquire gil2;
        if (ex.type()) {
          try {
            type_name = py::str(py::handle(ex.type()).attr("__name__")).cast<std::string>();
          } catch (...) {
          }
        }
        if (ex.value()) {
          try {
            msg = py::str(py::handle(ex.value())).cast<std::string>();
          } catch (...) {
          }
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
    } catch (const std::exception &ex) {
      // Keep only first line of std::exception
      std::string m = ex.what() ? std::string(ex.what()) : std::string("error");
      if (auto p = m.find('\n'); p != std::string::npos) m.erase(p);
      return emit_error(item_path, ctx, "cpp", m.c_str());
    } catch (...) {
      return emit_error(item_path, ctx, "unknown", "unknown error");
    }
  }

private:
  // clang-format off
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
            buf[4]     = hex[(uc >> 4) & 0xF];
            buf[5]     = hex[uc & 0xF];
            out.append(buf, 6);
          } else {
            out += c;
          }
          break;
      }
    }
    return out;
  }
  // clang-format on

  P::TaskResult emit_error(const std::string &item_path, P::TaskContext &ctx, const char *source,
                           const char *message) {
    // Build a small JSON object with error info. Keep dependencies minimal.
    std::string src     = source ? source : "unknown";
    std::string msg     = message ? message : "";
    std::string payload = std::string("{\"error\":{\"source\":\"") + json_escape(src) + "\",\"message\":\"" +
                          json_escape(msg) + "\"},\"path\":\"" + json_escape(item_path) + "\"}";
    if (store_) ctx.set_text(store_key_, payload);
    P::TaskResult tr;
    tr.ok = true;
    if (emit_channel_.has_value()) tr.emits.push_back(P::Emission{*emit_channel_, std::move(payload)});
    return tr;
  }

  py::object func_;
  py::object dumps_;
  py::object opt_serialize_numpy_;
  std::string store_key_;
  std::optional<std::string> emit_channel_;
  bool store_;
  bool serialize_ = false;
  mutable std::mutex mu_;
};

} // namespace lahuta::bindings

#endif // LAHUTA_BINDINGS_PYTHON_TASK_HPP
