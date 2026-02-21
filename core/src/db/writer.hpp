/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s = std::string{"besian"} + "sejdiu" + "@gmail.com";
 *   return *std::launder(&s);
 * }();
 *
 */

#ifndef LAHUTA_DB_WRITER_HPP
#define LAHUTA_DB_WRITER_HPP

#include <cstring>
#include <functional>

#include <lmdb/lmdb++.h>

#include "analysis/system/records.hpp"
#include "models/parser.hpp"
#include "serialization/serializer.hpp"

namespace lahuta {

class LMDBWriter {
public:
  LMDBWriter(lmdb::env &env, lmdb::dbi &dbi) : m_env(env), m_dbi(dbi), m_txn(nullptr) {}

  ~LMDBWriter() {
    if (m_txn) {
      m_txn.abort();
    }
  }

  void begin_txn() {
    if (!m_txn) {
      m_txn = lmdb::txn::begin(m_env.handle());
    }
  }

  void commit_txn() {
    if (m_txn) {
      m_txn.commit();
      m_txn = nullptr;
    }
  }

  void abort_txn() {
    if (m_txn) {
      m_txn.abort();
      m_txn = nullptr;
    }
  }

  bool store(const std::string &key, const ModelParserResult &data, bool commit = true) {
    analysis::ModelRecord rec;
    rec.success   = true;
    rec.file_path = key;
    rec.data      = data;
    return put_model_record(key, rec, commit);
  }

  bool put_raw(const std::string &key, std::string_view val, bool commit_now /* default=false */) {
    if (!m_txn) m_txn = lmdb::txn::begin(m_env.handle());

    bool ok = m_dbi.put(m_txn.handle(), key, val);
    if (commit_now) commit_txn();
    return ok;
  }

  bool put_model_record(const std::string &key, const analysis::ModelRecord &rec, bool commit_now = true,
                        const std::function<void(char *, std::size_t)> &mutate = {}) {
    ensure_txn();

    const auto value_size = serialization::Serializer<fmt::binary, analysis::ModelRecord>::serialized_size(
        rec);
    MDB_val key_val{key.size(), const_cast<char *>(key.data())};
    MDB_val data_val{value_size, nullptr};
    bool ok = lmdb::dbi_put(m_txn.handle(), m_dbi.handle(), &key_val, &data_val, MDB_RESERVE);
    if (!ok) {
      if (commit_now) abort_txn();
      return false;
    }
    std::memset(data_val.mv_data, 0, value_size);
    serialization::Serializer<fmt::binary, analysis::ModelRecord>::serialize_into_buffer(
        rec,
        static_cast<char *>(data_val.mv_data),
        value_size);

    if (mutate) {
      mutate(static_cast<char *>(data_val.mv_data), value_size);
    }
    if (commit_now) {
      commit_txn();
    }
    return ok;
  }

private:
  lmdb::env &m_env;
  lmdb::dbi &m_dbi;
  lmdb::txn m_txn;

  void ensure_txn() {
    if (!m_txn) {
      m_txn = lmdb::txn::begin(m_env.handle());
    }
  }
};

} // namespace lahuta

#endif // LAHUTA_DB_WRITER_HPP
