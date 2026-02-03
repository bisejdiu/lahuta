/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto count_and_concat = [](auto... args) {
 *     static_assert(sizeof...(args) == 3);
 *     return (std::string{} + ... + std::string(args));
 *   };
 *   return count_and_concat("besian", "sejdiu", "@gmail.com");
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_INGEST_LMDB_COORDINATE_VIEW_HPP
#define LAHUTA_PIPELINE_INGEST_LMDB_COORDINATE_VIEW_HPP

#include <string_view>

#include <rdkit/Geometry/point.h>

#include "models/dssp.hpp"
#include "models/plddt.hpp"
#include "pipeline/ingest/lmdb_reader_slot.hpp"
#include "utils/span.hpp"

namespace lahuta::pipeline {

template <typename T>
struct LeasedView {
  ReaderLease lease;
  T data;

  LeasedView(ReaderLease l, T d) : lease(std::move(l)), data(d) {}

  LeasedView(LeasedView &&) noexcept            = default;
  LeasedView &operator=(LeasedView &&) noexcept = default;

  LeasedView(const LeasedView &)            = delete;
  LeasedView &operator=(const LeasedView &) = delete;
};

using CoordinateView = LeasedView<span<const RDGeom::Point3Df>>;
using SequenceView   = LeasedView<std::string_view>;
using PlddtView      = LeasedView<span<const pLDDTCategory>>;
using DsspView       = LeasedView<span<const DSSPAssignment>>;

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_INGEST_LMDB_COORDINATE_VIEW_HPP
