#ifndef LAHUTA_SOURCES_LMDB_COORDINATE_VIEW_HPP
#define LAHUTA_SOURCES_LMDB_COORDINATE_VIEW_HPP

#include <string_view>

#include <rdkit/Geometry/point.h>

#include "models/dssp.hpp"
#include "models/plddt.hpp"
#include "sources/lmdb_reader_slot.hpp"
#include "span.hpp"

// clang-format off
namespace lahuta {

template <typename T>
struct LeasedView {
  ReaderLease lease;
  T data;

  LeasedView(ReaderLease l, T d) : lease(std::move(l)), data(d) {}

  LeasedView(LeasedView &&) noexcept = default;
  LeasedView &operator=(LeasedView &&) noexcept = default;

  LeasedView(const LeasedView &) = delete;
  LeasedView &operator=(const LeasedView &) = delete;
};

using CoordinateView = LeasedView<span<const RDGeom::Point3Df>>;
using SequenceView   = LeasedView<std::string_view>;
using PlddtView      = LeasedView<span<const pLDDTCategory>>;
using DsspView       = LeasedView<span<const DSSPAssignment>>;

} // namespace lahuta

#endif // LAHUTA_SOURCES_LMDB_COORDINATE_VIEW_HPP
