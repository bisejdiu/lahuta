#ifndef LAHUTA_COMPUTE_PARAMETERS_HPP
#define LAHUTA_COMPUTE_PARAMETERS_HPP

#include <memory>

#include "compute/_defs.hpp"

// clang-format off
namespace lahuta::topology::compute {

// Interface class encoding the parameters of a computation
struct ParameterInterface {
  using TypeId = u8;

  virtual ~ParameterInterface() = default;
  virtual std::unique_ptr<ParameterInterface> clone() const = 0;

  TypeId type_id() const noexcept { return type_id_; }

protected:
  explicit ParameterInterface(TypeId id) noexcept : type_id_(id) {}

private:
  TypeId type_id_;
};

/// Base class for all parameter types. Forces Derived to provide a TYPE_ID (used to cheaply identify the parameter type)
template <typename Derived>
struct ParameterBase : public ParameterInterface {
  ParameterBase() noexcept : ParameterInterface(Derived::TYPE_ID) {}
  explicit ParameterBase(TypeId id) noexcept : ParameterInterface(id) {}

  std::unique_ptr<ParameterInterface> clone() const override {
    return std::make_unique<Derived>(static_cast<const Derived&>(*this));
  }
};

} // namespace lahuta::topology::compute

#endif // LAHUTA_COMPUTE_PARAMETERS_HPP
