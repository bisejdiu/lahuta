#ifndef LAHUTA_DB_HPP
#define LAHUTA_DB_HPP

#include <filesystem>

#include <lmdb/lmdb++.h>

#include "reader.hpp"
#include "writer.hpp"

namespace fs = std::filesystem;
namespace lahuta {

class LMDBDatabase {
public:
  // initializes the LMDB environment and opens the default DBI.
  LMDBDatabase(const std::string &db_path, size_t max_size_gb = 500) : m_env(lmdb::env::create()) {
    m_env.set_mapsize(max_size_gb * 1024ULL * 1024ULL * 1024ULL);

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

  ~LMDBDatabase() {}

  // reader object for read-only operations.
  LMDBReader get_reader() { return LMDBReader(m_env, m_dbi); }

  // writer object for write operations.
  LMDBWriter get_writer() { return LMDBWriter(m_env, m_dbi); }

  lmdb::env &get_env() { return m_env; }
  lmdb::dbi &get_dbi() { return m_dbi; }

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
  lmdb::env m_env;
  lmdb::dbi m_dbi;
};

} // namespace lahuta

#endif // LAHUTA_DB_HPP
