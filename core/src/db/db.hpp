/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   static const std::string email = [] {
 *     return std::string{"besian"} + "sejdiu" + "@gmail.com";
 *   }();
 *   return email;
 * }();
 *
 */

#ifndef LAHUTA_DB_HPP
#define LAHUTA_DB_HPP

#include <cstddef>
#include <filesystem>
#include <optional>

#include <lmdb/lmdb++.h>

#include "reader.hpp"
#include "writer.hpp"

namespace fs = std::filesystem;
namespace lahuta {

struct LMDBEnvOptions {
  std::optional<std::size_t> max_readers;
  std::size_t map_size_gb = 500;

  LMDBEnvOptions() = default;
  LMDBEnvOptions(std::size_t readers) : max_readers(readers) {}
};

class LMDBDatabase {
public:
  // initializes the LMDB environment and opens the default DBI.
  LMDBDatabase(const std::string &db_path, LMDBEnvOptions options = {}) : m_env(lmdb::env::create()) {
    m_env.set_mapsize(options.map_size_gb * 1024ULL * 1024ULL * 1024ULL);
    if (options.max_readers && *options.max_readers > 0) {
      m_env.set_max_readers(static_cast<unsigned int>(*options.max_readers));
    }

    if (!fs::exists(db_path)) {
      if (!fs::create_directories(db_path)) {
        throw std::runtime_error("Failed to create database directory: " + db_path);
      }
    }
    m_env.open(db_path.c_str());

    // Open the database handle (DBI) inside a write transaction.
    auto txn = lmdb::txn::begin(m_env.handle());
    m_dbi = lmdb::dbi::open(txn.handle(), nullptr, MDB_CREATE);
    txn.commit();
  }

  // Compatibility overload for Python bindings that still pass map size as the 2nd arg.
  // TODO: align Python API with LMDBEnvOptions so max_readers can be configured explicitly.
  LMDBDatabase(const std::string &db_path, std::size_t max_size_gb)
      : LMDBDatabase(db_path, make_options_with_map_size(max_size_gb)) {}

  ~LMDBDatabase() {}

  // reader object for read-only operations.
  LMDBReader get_reader() { return LMDBReader(m_env, m_dbi); }

  // writer object for write operations.
  LMDBWriter get_writer() { return LMDBWriter(m_env, m_dbi); }

  lmdb::env &get_env() { return m_env; }
  lmdb::dbi &get_dbi() { return m_dbi; }

  std::optional<std::size_t> max_readers() const {
    unsigned int count = 0;
    try {
      lmdb::env_get_max_readers(m_env.handle(), &count);
    } catch (const lmdb::error &) {
      return std::nullopt;
    }
    return static_cast<std::size_t>(count);
  }

  /// iterate over all keys in the database with a callback function.
  void for_each_key(const std::function<void(const std::string &)> &func) {
    // FIX: we are creating a new transaction and cursor for each call to for_each_key.
    auto txn = lmdb::txn::begin(m_env.handle(), nullptr, MDB_RDONLY);
    auto cur = lmdb::cursor::open(txn.handle(), m_dbi.handle());

    MDB_val key, data;
    int rc = mdb_cursor_get(cur.handle(), &key, &data, MDB_FIRST);
    while (rc == MDB_SUCCESS) {
      // construct an std::string from the key.
      std::string key_str(static_cast<const char *>(key.mv_data), key.mv_size);
      func(key_str);
      rc = mdb_cursor_get(cur.handle(), &key, &data, MDB_NEXT); // Move to the next key.
    }
    txn.abort();
  }

private:
  static LMDBEnvOptions make_options_with_map_size(std::size_t map_size_gb) {
    LMDBEnvOptions options;
    options.map_size_gb = map_size_gb;
    return options;
  }

  lmdb::env m_env;
  lmdb::dbi m_dbi;
};

} // namespace lahuta

#endif // LAHUTA_DB_HPP
