#ifndef LAHUTA_PIPELINE_SOURCES_DB_KEY_SOURCE_HPP
#define LAHUTA_PIPELINE_SOURCES_DB_KEY_SOURCE_HPP

#include "db/db.hpp"
#include "logging.hpp"
#include <optional>
#include <string>
#include <vector>

// clang-format off
namespace lahuta::sources {

/// provides keys from an LMDB database
class DBKeySource {
public:
  using value_type = std::string;

  explicit DBKeySource(LMDBDatabase &database, std::size_t batch_size = 1024)
    : db_(database), batch_size_(batch_size) {
    load_next_batch();
  }

  [[nodiscard]] std::optional<std::string> next() {
    if (current_index_ == current_batch_.size() && !load_next_batch())
      return std::nullopt;

    return current_batch_[current_index_++];
  }

  [[nodiscard]] std::size_t size() const noexcept { return current_batch_.size(); }

private:
  bool load_next_batch() {
    current_batch_.clear();
    current_index_ = 0;

    try {
      current_batch_.reserve(batch_size_);

      // ro transaction
      auto txn    = lmdb::txn::begin(db_.get_env(), nullptr, MDB_RDONLY);
      auto dbi    = lmdb::dbi::open(txn, nullptr);
      auto cursor = lmdb::cursor::open(txn, dbi);

      std::string_view key, data;
      if (last_key_.empty()) {
        key = "";
        if (!cursor.get(key, data, MDB_FIRST)) {
          cursor.close();
          txn.abort();
          return false; // empty db
        }
      } else {
        // start after the last key we processed
        key = last_key_;
        if (!cursor.get(key, data, MDB_SET)) {
          cursor.close();
          txn.abort();
          return false; // should not happen?!
        }

        // go to the next key
        if (!cursor.get(key, data, MDB_NEXT)) {
          cursor.close();
          txn.abort();
          return false; // no more keys
        }
      }

      // collect keys up to batch_size
      std::size_t count = 0;
      do {
        current_batch_.push_back(std::string(key));
        count++;

        if (count >= batch_size_) {
          break;
        }
      } while (cursor.get(key, data, MDB_NEXT));

      if (!current_batch_.empty()) {
        last_key_ = current_batch_.back();
      }

      cursor.close();
      txn.abort();

      return !current_batch_.empty();

    } catch (const lmdb::error &e) {
      Logger::get_logger()->error("LMDB error loading batch of keys: {}", e.what());
    } catch (const std::exception &e) {
      Logger::get_logger()->error("Error loading batch of keys: {}", e.what());
    }

    return false;
  }

  LMDBDatabase &db_;
  std::size_t  batch_size_;
  std::size_t  current_index_ = 0;
  std::string  last_key_; // pagination
  std::vector<std::string>  current_batch_;
};

} // namespace lahuta::sources

#endif // LAHUTA_PIPELINE_SOURCES_DB_KEY_SOURCE_HPP
