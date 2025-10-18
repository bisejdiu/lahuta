#ifndef LAHUTA_BINDINGS_PYTHON_PROCESS_TASK_HPP
#define LAHUTA_BINDINGS_PYTHON_PROCESS_TASK_HPP

#include <condition_variable>
#include <cstdint>
#include <memory>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "pipeline/compute/parameters.hpp"
#include "pipeline/dynamic/types.hpp"
#include "pipeline/process_pool.hpp"
#include "pipeline/worker_protocol.hpp"

// clang-format off
namespace py = pybind11;

namespace lahuta::bindings {
using namespace lahuta::pipeline::dynamic;
using lahuta::pipeline::compute::BuildTopologyParams;
using lahuta::pipeline::compute::SystemReadParams;

namespace {

//
// ConcurrencyLimiter bounds the number of concurrent Python tasks submitted to
// the process pool. Preserves StageManager's threaded execution for C++ tasks
// while preventing the Python pool from being oversubscribed.
//
class ConcurrencyLimiter {
public:
  explicit ConcurrencyLimiter(std::size_t capacity) : capacity_(capacity), available_(capacity) {}

  class Permit {
  public:
    explicit Permit(ConcurrencyLimiter& parent) : parent_(parent) {}
    Permit(const Permit&) = delete;
    Permit& operator=(const Permit&) = delete;
    Permit(Permit&&) = delete;
    Permit& operator=(Permit&&) = delete;
    ~Permit() { parent_.release(); }

  private:
    ConcurrencyLimiter& parent_;
  };

  std::unique_ptr<Permit> acquire() {
    if (capacity_ == 0) return nullptr;
    std::unique_lock<std::mutex> lock(mutex_);
    cv_.wait(lock, [&] { return available_ > 0; });
    --available_;
    return std::make_unique<Permit>(*this);
  }

  void resize(std::size_t capacity) {
    std::unique_lock<std::mutex> lock(mutex_);
    std::size_t running = capacity_ >= available_ ? capacity_ - available_ : 0;
    capacity_ = capacity;
    if (capacity_ <= running) {
      available_ = 0;
    } else {
      available_ = capacity_ - running;
    }
    lock.unlock();
    cv_.notify_all();
  }

private:
  void release() {
    std::unique_lock<std::mutex> lock(mutex_);
    if (available_ < capacity_) ++available_;
    lock.unlock();
    cv_.notify_one();
  }

  std::mutex mutex_;
  std::condition_variable cv_;
  std::size_t capacity_;
  std::size_t available_;
};

} // namespace

class PyProcessTask : public ITask {
public:
  PyProcessTask(std::string name,
                std::string module,
                std::string qualname,
                py::object serialized_callable,
                std::optional<std::string> emit_channel,
                bool store,
                const SystemReadParams*    sys_params,
                const BuildTopologyParams* topo_params)
      : store_key_(std::move(name)),
        module_(std::move(module)),
        qualname_(std::move(qualname)),
        serialized_callable_(std::move(serialized_callable)),
        emit_channel_(std::move(emit_channel)),
        store_(store),
        sys_params_(sys_params),
        topo_params_(topo_params) {
    if ((module_.empty() || qualname_.empty()) && serialized_callable_.is_none()) {
      throw std::invalid_argument("PyProcessTask: module/qualname missing and no serialized callable provided");
    }
  }

  TaskResult run(const std::string& item_path, TaskContext& ctx) override {
    try {
      py::gil_scoped_acquire gil;

      py::dict initial_store;
      for (const auto& kv : ctx.texts()) {
        initial_store[py::str(kv.first)] = py::str(kv.second);
      }

      py::dict initial_bytes;
      for (const auto& kv : ctx.bytes()) {
        // Cast to Python bytes to preserve binary payloads
        initial_bytes[py::str(kv.first)] = py::bytes(kv.second);
      }

      //
      // TODO: Parameter updates require more work. This is minimal for now, but if we make changes
      // to the params structs, we need to update this code as well - which we'll definitely forget.
      //
      py::dict system_params;
      if (sys_params_) {
        system_params["is_model"] = sys_params_->is_model;
      }

      py::dict topology_params;
      if (topo_params_) {
        topology_params["flags"] = static_cast<std::uint32_t>(topo_params_->flags);
        topology_params["atom_typing_method"] = static_cast<int>(topo_params_->atom_typing_method);
      }

      std::shared_ptr<ConcurrencyLimiter> limiter_snapshot;
      {
        std::lock_guard<std::mutex> lock(limiter_mutex_);
        limiter_snapshot = limiter_;
      }
      std::unique_ptr<ConcurrencyLimiter::Permit> permit;
      if (limiter_snapshot) {
        permit = limiter_snapshot->acquire();
      }

      auto response = PyProcessPool::instance().invoke(
        module_, qualname_, serialized_callable_, item_path, initial_store, initial_bytes, system_params, topology_params, timeout_seconds_);

      if (response.timed_out)      { return emit_error(item_path, ctx, "python", "Timeout");        }
      if (response.worker_crashed) { return emit_error(item_path, ctx, "python", "Worker crashed"); }
      if (!response.success)       { return emit_error(item_path, ctx, "python", "unknown error");  }

      py::object response_obj = response.value;
      if (!py::isinstance<py::dict>(response_obj)) {
        return emit_error(item_path, ctx, "python", "worker returned non-dict response");
      }
      py::dict response_dict = response_obj.cast<py::dict>();

      worker_protocol::apply_store_text (response_dict, ctx);
      worker_protocol::apply_store_bytes(response_dict, ctx);

      auto payload_result = worker_protocol::extract_payload(response_dict);
      bool ok_flag        = worker_protocol::extract_ok_flag(response_dict);

      // Maybe store payload
      if (store_ && payload_result.found) {
        if (payload_result.is_binary) {
          ctx.set_bytes(store_key_, payload_result.data);
        } else {
          ctx.set_text (store_key_, payload_result.data);
        }
      }


      TaskResult tr; // Build result
      tr.ok = ok_flag;
      if (payload_result.found && emit_channel_.has_value()) {
        tr.emits.push_back(Emission{*emit_channel_, std::move(payload_result.data)});
      }

      return tr;
    } catch (const py::error_already_set& ex) {
      std::string type_name = "PythonError";
      std::string msg;
      try {
        if (ex.type()) {
          type_name = py::str(py::handle(ex.type()).attr("__name__")).cast<std::string>();
        }
        if (ex.value()) {
          msg = py::str(py::handle(ex.value())).cast<std::string>();
        }
      } catch (...) {
        msg = ex.what() ? std::string(ex.what()) : std::string("error");
      }
      if (msg.empty()) {
        msg = ex.what() ? std::string(ex.what()) : std::string("error");
      }
      if (auto p = msg.find('\n'); p != std::string::npos) msg.erase(p);
      return emit_error(item_path, ctx, "python", (type_name + ": " + msg).c_str());
    } catch (const std::exception& ex) {
      std::string msg = ex.what() ? std::string(ex.what()) : std::string("error");
      if (auto p = msg.find('\n'); p != std::string::npos) msg.erase(p);
      return emit_error(item_path, ctx, "cpp", msg.c_str());
    } catch (...) {
      return emit_error(item_path, ctx, "unknown", "unknown error");
    }
  }

private:
  static std::string json_escape(std::string s) {
    static const char* hex = "0123456789ABCDEF";

    std::string out;
    out.reserve(s.size() + 8);
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

  TaskResult emit_error(const std::string& item_path,
                        TaskContext& ctx,
                        const char* source,
                        const char* message) {
    std::string src = source ? source : "unknown";
    std::string msg = message ? message : "";
    std::string payload = std::string("{\"error\":{\"source\":\"") + json_escape(src) +
                          "\",\"message\":\"" + json_escape(msg) +
                          "\"},\"path\":\"" + json_escape(item_path) + "\"}";
    if (store_) {
      ctx.set_text(store_key_, payload);
    }
    TaskResult tr; tr.ok = true;
    if (emit_channel_.has_value()) {
      tr.emits.push_back(Emission{*emit_channel_, std::move(payload)});
    }
    return tr;
  }

  std::string store_key_;
  std::string module_;
  std::string qualname_;
  py::object serialized_callable_;
  std::optional<std::string> emit_channel_;
  bool store_;
  const SystemReadParams*    sys_params_  = nullptr;
  const BuildTopologyParams* topo_params_ = nullptr;

public:
  static void set_concurrency_limit(std::size_t limit) {
    std::lock_guard<std::mutex> lock(limiter_mutex_);
    if (limit == 0) {
      limiter_.reset();
      return;
    }
    if (!limiter_) {
      limiter_ = std::make_shared<ConcurrencyLimiter>(limit);
    } else {
      limiter_->resize(limit);
    }
  }

  static void set_timeout(double seconds) {
    timeout_seconds_ = seconds > 0.0 ? seconds : 0.0;
  }

private:
  inline static std::shared_ptr<ConcurrencyLimiter> limiter_;
  inline static std::mutex limiter_mutex_;
  inline static double timeout_seconds_ = 300.0;
};

} // namespace lahuta::bindings

#endif // LAHUTA_BINDINGS_PYTHON_PROCESS_TASK_HPP
