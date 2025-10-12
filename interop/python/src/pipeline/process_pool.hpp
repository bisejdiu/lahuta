#ifndef LAHUTA_BINDINGS_PROCESS_POOL_HPP
#define LAHUTA_BINDINGS_PROCESS_POOL_HPP

#include <mutex>
#include <stdexcept>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

// clang-format off
namespace py = pybind11;

namespace lahuta::bindings {

//
// Lightweight wrapper around Python's multiprocessing.Pool that's controlled
// from our C++ pipeline runtime. The pool lives for the duration of a pipeline
// run and is shared by all PyProcessTask instances.
//
class PyProcessPool {
public:
  //
  // InvokeResult captures the outcome of a single apply_async/get call.
  // We distinguish between timeout and worker crash, but Python's
  // multiprocessing.Pool does not seem to reliably detect the difference in practice:
  // when a worker crashes (e.g., os._exit(1)), Pool.apply_async().get(timeout)
  // waits the full timeout period and then raises TimeoutError, exactly like a
  // genuine timeout. We check for both builtin TimeoutError and
  // multiprocessing.context.TimeoutError, and only detect crashes via
  // EOFError/MaybeEncodingError. Callers should treat both timed_out and worker_crashed
  // as fatal errors for the current task.        - Besian, October 2025
  //
  struct InvokeResult {
    bool success        = false;
    bool timed_out      = false;
    bool worker_crashed = false;
    py::object value;
  };

  static PyProcessPool& instance() {
    static PyProcessPool pool;
    return pool;
  }

  void configure(std::size_t processes) {
    if (processes == 0) {
      throw std::invalid_argument("PyProcessPool.configure: processes must be >= 1");
    }

    std::lock_guard<std::mutex> lock(mutex_);

    py::gil_scoped_acquire gil;
    if (worker_may_have_crashed_) {
      shutdown_unlocked(true, true);
    }
    if (!pool_.is_none() && configured_processes_ == processes) return;

    shutdown_unlocked(true, worker_may_have_crashed_);

    py::module_ mp = py::module_::import("multiprocessing");
    py::object ctx = mp.attr("get_context")("spawn");
    pool_ = ctx.attr("Pool")(processes);

    py::module_ worker_mod = py::module_::import("lahuta.pipeline._process_worker");
    worker_func_ = worker_mod.attr("execute_callable");
    ctx_ = ctx;

    configured_processes_ = processes;
    worker_may_have_crashed_ = false;
  }

  void shutdown() {
    std::lock_guard<std::mutex> lock(mutex_);
    py::gil_scoped_acquire gil;
    shutdown_unlocked(true, worker_may_have_crashed_);
  }

  InvokeResult invoke(const std::string& module,
                      const std::string& qualname,
                      const py::object& serialized_callable,
                      const std::string& item_path,
                      const py::dict& initial_store,
                      const py::dict& initial_bytes,
                      const py::dict& system_params,
                      const py::dict& topology_params,
                      double timeout_seconds) {
    py::gil_scoped_acquire gil;

    py::object pool;
    py::object worker;
    {
      std::lock_guard<std::mutex> lock(mutex_);
      if (pool_.is_none() || worker_func_.is_none()) {
        throw std::runtime_error("PyProcessPool.invoke: pool not configured");
      }
      pool   = pool_;
      worker = worker_func_;
    }

    py::tuple args = py::make_tuple(module,
                                    qualname,
                                    serialized_callable,
                                    item_path,
                                    initial_store,
                                    initial_bytes,
                                    system_params,
                                    topology_params);
    py::object async_result = pool.attr("apply_async")(worker, args); // pool.apply_async(worker, args)
    try {
      py::object result;
      if (timeout_seconds > 0.0) {
        result = async_result.attr("get")(timeout_seconds);
      } else {
        result = async_result.attr("get")();
      }
      return InvokeResult{true, false, false, result};
    } catch (const py::error_already_set& ex) {

      //
      // Pythons multiprocessing raises multiprocessing.context.TimeoutError
      // for both genuine timeouts and for worker crashes. We check
      // both exception types and return timed_out=true for either.  - Besian, October 2025
      //
      py::object mp_timeout_error;
      try {
        mp_timeout_error = py::module_::import("multiprocessing.context").attr("TimeoutError");
      } catch (...) {
        mp_timeout_error = py::none();
      }

      if (ex.matches(PyExc_TimeoutError) || (!mp_timeout_error.is_none() && ex.matches(mp_timeout_error.ptr()))) {
        // Worker may have crashed (os._exit doesn't send proper signals)
        std::lock_guard<std::mutex> lock(mutex_);
        worker_may_have_crashed_ = true;
        return InvokeResult{false, true, false, py::none()};
      }
      py::object maybe_encoding_error;
      try {
        maybe_encoding_error = py::module_::import("multiprocessing.pool").attr("MaybeEncodingError");
      } catch (...) {
        maybe_encoding_error = py::none();
      }

      if (ex.matches(PyExc_EOFError) || (!maybe_encoding_error.is_none() && ex.matches(maybe_encoding_error.ptr()))) {
        // Only restart on actual worker crash (EOF indicates worker died)
        {
          std::lock_guard<std::mutex> lock(mutex_);
          worker_may_have_crashed_ = true;
        }
        restart_pool();
        return InvokeResult{false, false, true, py::none()};
      }

      throw; // re-raise other exceptions
    }
  }

private:
  PyProcessPool() : pool_(py::none()), worker_func_(py::none()), ctx_(py::none()) {}

  void shutdown_unlocked(bool reset_state = true, bool force_terminate = false) {
    if (pool_.is_none()) {
      if (reset_state) {
        configured_processes_ = 0;
        worker_func_ = py::none();
        ctx_         = py::none();
        worker_may_have_crashed_ = false;
      }
      return;
    }

    try {
      if (force_terminate) {
        // When workers may have crashed (os._exit), terminate immediately to avoid hanging on join()
        pool_.attr("terminate")();
        pool_.attr("join")();
      } else {
        pool_.attr("close")();
        pool_.attr("join")();
      }
    } catch (...) {
      try {
        pool_.attr("terminate")();
      } catch (...) {
        // ignore
      }
    }
    pool_ = py::none();
    if (reset_state) {
      worker_func_ = py::none();
      ctx_         = py::none();
      configured_processes_ = 0;
      worker_may_have_crashed_ = false;
    }
  }

  void restart_pool() {
    std::lock_guard<std::mutex> lock(mutex_);
    if (configured_processes_ == 0 || ctx_.is_none()) {
      configured_processes_ = 0;
      worker_func_ = py::none();
      pool_        = py::none();
      ctx_         = py::none();
      return;
    }

    auto processes = configured_processes_;
    shutdown_unlocked(false, worker_may_have_crashed_);
    pool_ = ctx_.attr("Pool")(processes);
    worker_may_have_crashed_ = false;
  }

  std::mutex mutex_;
  std::size_t configured_processes_ = 0;
  py::object pool_;
  py::object worker_func_;
  py::object ctx_;
  bool worker_may_have_crashed_ = false;
};

} // namespace lahuta::bindings

#endif // LAHUTA_BINDINGS_PROCESS_POOL_HPP
