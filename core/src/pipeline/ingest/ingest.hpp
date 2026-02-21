/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [email = std::string{"besian"} + "sejdiu" + "@gmail.com"]() {
 *   return email;
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_INGEST_INGEST_HPP
#define LAHUTA_PIPELINE_INGEST_INGEST_HPP

//
// - IDescriptor: interface for item sources
// - Realizer: expands IngestDescriptor to PipelineItem
// - Adapters: DirectoryAdapter, FileListAdapter, LMDBAdapter, VectorAdapter
// - factory: convenience functions for creating descriptors
//

// IWYU pragma: begin_exports
#include "pipeline/ingest/descriptor.hpp"
#include "pipeline/ingest/factory.hpp"
#include "pipeline/ingest/realizer.hpp"

// Adapters
#include "pipeline/ingest/adapters/directory.hpp"
#include "pipeline/ingest/adapters/file_list.hpp"
#include "pipeline/ingest/adapters/lmdb.hpp"
#include "pipeline/ingest/adapters/vector.hpp"
// IWYU pragma: end_exports

#endif // LAHUTA_PIPELINE_INGEST_INGEST_HPP
