#ifndef LAHUTA_IO_LMDB_BACKEND_HPP
#define LAHUTA_IO_LMDB_BACKEND_HPP

#include "db/writer.hpp"
#include <string>

namespace lahuta {

// clang-format off
template<class SerializerT>
class LMDBBackend {
  lahuta::LMDBWriter& w_;
public:
  explicit LMDBBackend(lahuta::LMDBWriter& w) : w_(w) {}

  void write(std::string_view key, std::string_view val) {
    w_.put_raw(std::string(key), val, /*commit=*/false);
  }

  void begin_txn()  noexcept { w_.begin_txn();  }
  void commit_txn() noexcept { w_.commit_txn(); }
};

} // namespace lahuta

#endif // LAHUTA_IO_LMDB_BACKEND_HPP
