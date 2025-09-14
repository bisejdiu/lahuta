#pragma once

#include <type_traits>

// clang-format off
namespace lahuta::topology::compute {

enum class Mut { ReadOnly, ReadWrite };
template <typename D, Mut M> class ComputeEngine;

/// provides access to the data being operated on and the engine that is executing the computation
template <typename DataT, Mut M>
class DataContext {
public:
  explicit DataContext(DataT& data) 
    : data_(data), engine_(nullptr) {}

  DataContext(DataT& data, ComputeEngine<DataT, M>* engine) 
    : data_(data), engine_(engine) {}

  ComputeEngine<DataT, M>* get_engine() const { return engine_; }

  /// access to data with appropriate constness
  typename std::conditional<M == Mut::ReadOnly, const DataT&, DataT&>::type data() const { return data_; }

  template <Mut MM = M>
  operator typename std::enable_if<MM == Mut::ReadWrite, DataContext<DataT, Mut::ReadOnly>>::type() const noexcept {
    // reinterpret_cast is safe because we're reading though the engine pointer
    return DataContext<DataT, Mut::ReadOnly>(data_, reinterpret_cast<ComputeEngine<DataT, Mut::ReadOnly>*>(engine_));
  }

private:
  DataT& data_;
  ComputeEngine<DataT, M>* engine_;
};

} // namespace lahuta::topology::compute
