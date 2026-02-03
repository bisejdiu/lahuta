/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::optional<std::string> e; e = std::string{"besian"};
 *   e = e.value_or("") + "sejdiu"; e = e.value_or("") + "@gmail.com";
 *   return e.value_or("");
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_INGEST_LMDB_READER_SLOT_HPP
#define LAHUTA_PIPELINE_INGEST_LMDB_READER_SLOT_HPP

#include <cstddef>
#include <optional>
#include <stdexcept>

#include <lmdb/lmdb++.h>

namespace lahuta::pipeline {

class ReaderLease;

class ReaderTxnSlot {
  friend class ReaderLease;

public:
  ReaderLease acquire(lmdb::env &env);
  lmdb::txn &txn();

  void reset() {
    while (leases_ > 0) {
      release();
    }
    count_ = 0;
  }

  // Returns the number of times a transaction has been begun on this slot.
  std::size_t count() const noexcept { return count_; }

private:
  void retain();
  void release();
  void prepare(lmdb::env &env);

  const lmdb::env *env_{nullptr};
  std::optional<lmdb::txn> txn_;
  std::size_t leases_{0};
  bool active_{false};
  std::size_t count_{0}; // nr. of times a txn has been begun on this slot.
};

class ReaderLease {
public:
  ReaderLease() = default;
  explicit ReaderLease(ReaderTxnSlot *slot);
  ReaderLease(const ReaderLease &other);
  ReaderLease(ReaderLease &&other) noexcept;
  ReaderLease &operator=(const ReaderLease &other);
  ReaderLease &operator=(ReaderLease &&other) noexcept;
  ~ReaderLease();

  lmdb::txn &txn() const;
  bool valid() const noexcept { return slot_ != nullptr; }

private:
  void release();

  ReaderTxnSlot *slot_{nullptr};
};

inline ReaderLease ReaderTxnSlot::acquire(lmdb::env &env) {
  prepare(env);
  return ReaderLease(this);
}

inline lmdb::txn &ReaderTxnSlot::txn() {
  if (!txn_) {
    throw std::runtime_error("ReaderTxnSlot: no active transaction");
  }
  return *txn_;
}

inline void ReaderTxnSlot::retain() { ++leases_; }

inline void ReaderTxnSlot::release() {
  if (leases_ == 0) return;
  --leases_;
  if (leases_ == 0) {
    if (active_ && txn_) {
      try {
        txn_->reset();
      } catch (const lmdb::error &) {
        // ignore reset failures on stale handles
      }
    }
    txn_.reset();
    active_ = false;
    env_    = nullptr;
  }
}

inline void ReaderTxnSlot::prepare(lmdb::env &env) {
  if (env_ && env_ != &env) {
    if (leases_ != 0) {
      throw std::runtime_error("ReaderTxnSlot: environment mismatch while transaction in use");
    }
    txn_.reset();
    active_ = false;
    env_    = &env;
  } else if (!env_) {
    env_ = &env;
  }

  // Always start a fresh read txn to avoid renewing potentially stale handles
  txn_.reset();
  try {
    txn_.emplace(lmdb::txn::begin(env.handle(), nullptr, MDB_RDONLY));
    active_ = true;
    ++count_; // Increment transaction begin count
  } catch (const lmdb::error &e) {
    if (e.code() == MDB_BAD_RSLOT) {
      throw std::runtime_error(
          "LMDB read transaction reuse failed (likely a zero-copy view is still alive on this thread). "
          "Release positions_view or make a copy of the positions before the next pipeline item.");
    }
    throw;
  }
}

inline ReaderLease::ReaderLease(ReaderTxnSlot *slot) : slot_(slot) {
  if (slot_) slot_->retain();
}

inline ReaderLease::ReaderLease(const ReaderLease &other) : slot_(other.slot_) {
  if (slot_) slot_->retain();
}

inline ReaderLease::ReaderLease(ReaderLease &&other) noexcept : slot_(other.slot_) { other.slot_ = nullptr; }

inline ReaderLease &ReaderLease::operator=(const ReaderLease &other) {
  if (this != &other) {
    release();
    slot_ = other.slot_;
    if (slot_) slot_->retain();
  }
  return *this;
}

inline ReaderLease &ReaderLease::operator=(ReaderLease &&other) noexcept {
  if (this != &other) {
    release();
    slot_       = other.slot_;
    other.slot_ = nullptr;
  }
  return *this;
}

inline ReaderLease::~ReaderLease() { release(); }

inline lmdb::txn &ReaderLease::txn() const {
  if (!slot_) {
    throw std::runtime_error("ReaderLease: no active transaction");
  }
  return slot_->txn();
}

inline void ReaderLease::release() {
  if (slot_) {
    slot_->release();
    slot_ = nullptr;
  }
}

inline ReaderTxnSlot &reader_txn_slot() {
  thread_local ReaderTxnSlot slot;
  return slot;
}

inline ReaderLease acquire_reader_lease(lmdb::env &env) { return reader_txn_slot().acquire(env); }

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_INGEST_LMDB_READER_SLOT_HPP
