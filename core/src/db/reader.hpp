#ifndef LAHUTA_DB_READER_HPP
#define LAHUTA_DB_READER_HPP

#include "lmdb/lmdb++.h"
#include "models/parser.hpp"
#include "serialization/serializer.hpp"

// clang-format off
namespace lahuta {

class LMDBReader {
public:
  LMDBReader(lmdb::env &env, lmdb::dbi &dbi) : m_env(env), m_dbi(dbi) {}

  /// Retrieve a ModelParserResult based on a key. Returns false if key not found.
  bool retrieve(const std::string &key, ModelParserResult &result) {
    std::string raw;
    if (!fetch(key, raw)) {
      return false;
    }
    result = serialization::Serializer<fmt::binary, ModelParserResult>::deserialize(raw.data(), raw.size());
    return true;
  }

  bool fetch(std::string_view key, std::string &value) {
    auto txn = lmdb::txn::begin(m_env.handle(), nullptr, MDB_RDONLY);
    std::string_view view;
    if (!m_dbi.get(txn.handle(), key, view)) {
      txn.abort();
      return false;
    }
    value.assign(view.data(), view.size());
    txn.abort();
    return true;
  }

  template <class Rec, class FormatTag> bool fetch_typed(std::string_view key, Rec &out) {
    std::string raw;
    if (!fetch(key, raw)) return false;
    out = serialization::Serializer<FormatTag, Rec>::deserialize(raw.data(), raw.size());
    return true;
  }

private:
  lmdb::env &m_env;
  lmdb::dbi &m_dbi;
};

} // namespace lahuta

#endif // LAHUTA_DB_READER_HPP
