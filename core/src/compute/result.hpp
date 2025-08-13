#pragma once

#include "error.hpp"
#include <memory>
#include <optional>
#include <stdexcept>
#include <typeinfo>

// clang-format off
namespace lahuta::topology::compute {

/// encodes the result of a computation
class ComputationResult {
private:
  class ResultBase {
  public:
    virtual ~ResultBase() = default;
    virtual const std::type_info &get_type() const = 0;
  };

  template <typename T>
  class TypedResult : public ResultBase {
  public:
    explicit TypedResult(T value) : value_(std::move(value)) {}
    const std::type_info &get_type() const override { return typeid(T); }
    const T &get_value() const { return value_; }
    T move_value() { return std::move(value_); } // NOTE: not tested, likely very dangerous

  private:
    T value_;
  };

  std::shared_ptr<ResultBase>     result_;
  std::optional<ComputationError> error_;

public:
  // by default constructor creates an uninitialized result with an error
  ComputationResult()
      : result_(nullptr),
        error_(ComputationError("Uninitialized result", ComputationError::Severity::Error)) {}

  template <typename T>
  explicit ComputationResult(T value)
      : result_(std::make_shared<TypedResult<T>>(std::move(value))),
        error_(std::nullopt) {}

  explicit ComputationResult(ComputationError error)
      : result_(nullptr), error_(std::move(error)) {}

  bool is_success() const { return !has_error(); }
  bool has_error()  const { return error_.has_value(); }
  bool has_value()  const { return result_ != nullptr; }

  /// Get the error if one exists
  const ComputationError &error() const {
    if (!error_.has_value()) throw std::runtime_error("No error present in result");
    return error_.value();
  }

  /// Get the value. Throws if the result is an error or uninitialized.
  template <typename T>
  const T &get_value() const {
    if (has_error()) throw std::runtime_error("Cannot access value on error result: " + error().get_message());
    if (!result_)    throw std::runtime_error("Accessing uninitialized result");
    if (result_->get_type() != typeid(T)) throw std::bad_cast();

    auto *typed = static_cast<const TypedResult<T> *>(result_.get());
    return typed->get_value();
  }

  /// Get the type of the value. Throws if the result is an error or uninitialized.
  const std::type_info &get_type() const {
    if (has_error()) throw std::runtime_error("Cannot get type of error result");
    if (!result_)    throw std::runtime_error("Accessing uninitialized result");

    return result_->get_type();
  }

  /// Move the value out of the result. NOte that ComputationResult does not own the value.
  /// So this is a dangerous operation. Use with caution.
  template <typename T>
  T move_value() {
    if (has_error()) throw std::runtime_error("Cannot move value from error result");
    if (!result_)    throw std::runtime_error("Accessing uninitialized result");
    if (result_->get_type() != typeid(T)) throw std::bad_cast();

    auto *typed = static_cast<TypedResult<T> *>(result_.get());
    return typed->move_value();
  }
};

} // namespace lahuta::topology::compute
