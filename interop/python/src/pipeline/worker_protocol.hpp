/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto make = []() noexcept(noexcept(std::string{})) {
 *     return std::string{"besian"} + "sejdiu" + "@gmail.com";
 *   };
 *   static_assert(noexcept(make()) == noexcept(std::string{}));
 *   return make();
 * }();
 *
 */

#ifndef LAHUTA_BINDINGS_WORKER_PROTOCOL_HPP
#define LAHUTA_BINDINGS_WORKER_PROTOCOL_HPP

#include <string>

#include <pybind11/pybind11.h>

#include "pipeline/task/context.hpp"

namespace py = pybind11;

namespace lahuta::bindings::worker_protocol {
namespace P = lahuta::pipeline;

/// Result of payload extraction from worker response
struct PayloadResult {
  bool found = false;
  std::string data;
  bool is_binary = false; // true if extracted from payload_bytes
};

// Apply text store updates from worker response to context
// Looks for "store_text" key containing dict of string:string mappings
inline void apply_store_text(const py::dict &response, P::TaskContext &ctx) {
  if (!response.contains("store_text")) return;

  py::object store_obj = response["store_text"];
  if (!py::isinstance<py::dict>(store_obj)) return;

  py::dict store_dict = store_obj.cast<py::dict>();
  for (auto item : store_dict) {
    std::string key = py::str(item.first);
    std::string val = py::str(item.second);
    ctx.set_text(key, std::move(val));
  }
}

// Apply binary store updates from worker response to context
// Looks for "store_bytes" key containing dict of string:bytes mappings
inline void apply_store_bytes(const py::dict &response, P::TaskContext &ctx) {
  if (!response.contains("store_bytes")) return;

  py::object store_obj = response["store_bytes"];
  if (!py::isinstance<py::dict>(store_obj)) return;

  py::dict store_dict = store_obj.cast<py::dict>();
  for (auto item : store_dict) {
    std::string key = py::str(item.first);
    std::string val;
    py::object vobj = py::reinterpret_borrow<py::object>(item.second);
    if (py::isinstance<py::bytes>(vobj)) {
      val = vobj.cast<std::string>();
    } else {
      // Fallback to str() encoding if bytes isn't provided
      val = py::str(vobj).cast<std::string>();
    }
    ctx.set_bytes(key, std::move(val));
  }
}

// Priority: payload_bytes (binary) -> payload (auto-detect) -> payload_text (text)
inline PayloadResult extract_payload(const py::dict &response) {
  PayloadResult result;

  if (response.contains("payload_bytes")) {
    py::object payload_obj = response["payload_bytes"];
    if (!payload_obj.is_none()) {
      result.data      = payload_obj.cast<std::string>();
      result.found     = true;
      result.is_binary = true;
      return result;
    }
  }

  if (response.contains("payload")) {
    py::object payload_obj = response["payload"];
    if (!payload_obj.is_none()) {
      result.data      = payload_obj.cast<std::string>();
      result.found     = true;
      result.is_binary = py::isinstance<py::bytes>(payload_obj);
      return result;
    }
  }

  if (response.contains("payload_text")) {
    py::object payload_obj = response["payload_text"];
    if (!payload_obj.is_none()) {
      result.data      = payload_obj.cast<std::string>();
      result.found     = true;
      result.is_binary = false;
      return result;
    }
  }

  return result; // not found
}

// Extract ok flag from worker response, defaulting to true if missing
inline bool extract_ok_flag(const py::dict &response) {
  if (!response.contains("ok")) return true;

  try {
    return response["ok"].cast<bool>();
  } catch (...) {
    return true; // Default to true on cast failure
  }
}

} // namespace lahuta::bindings::worker_protocol

#endif // LAHUTA_BINDINGS_WORKER_PROTOCOL_HPP
