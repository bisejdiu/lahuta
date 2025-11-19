#ifndef LAHUTA_PIPELINE_DYNAMIC_SINK_LMDB_HPP
#define LAHUTA_PIPELINE_DYNAMIC_SINK_LMDB_HPP

#include <cstdint>
#include <cstring>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <string_view>

#include "db/db.hpp"
#include "db/writer.hpp"
#include "pipeline/dynamic/sink_iface.hpp"
#include "pipeline/dynamic/types.hpp"

// clang-format off
namespace lahuta::pipeline::dynamic {

//
// LmdbSink: Writes emissions to an LMDB database.
//
// - Payload is a binary serialization of analysis::system::ModelRecord
//   produced by serialization::Serializer<fmt::binary, ModelRecord>.
// - The key is extracted by reading the file_path field from the serialized
//   payload header without fully deserializing the model blob.
//
class LmdbSink : public IDynamicSink {
public:
  explicit LmdbSink(std::shared_ptr<lahuta::LMDBDatabase> db, std::size_t batch_size)
    : writer_(db->get_env(), db->get_dbi()), batch_size_(batch_size) {
    if (batch_size_ == 0) batch_size_ = 1;
  }

  void write(EmissionView e) override {
    std::lock_guard<std::mutex> lk(mu_);
    // Payload layout (see serialization/specializations/model.hpp):
    // [1 byte success][4 bytes path_len][path][4 bytes blob_len][blob]
    const char* data = e.payload.data();
    const std::size_t n = e.payload.size();
    if (n < 1 + sizeof(uint32_t)) {
      throw std::runtime_error("LmdbSink: payload too small");
    }
    std::size_t off = 1; // skip success
    uint32_t path_len = 0;
    std::memcpy(&path_len, data + off, sizeof(path_len));
    off += sizeof(path_len);
    if (off + path_len > n) {
      throw std::runtime_error("LmdbSink: corrupted payload (path_len)");
    }
    const std::string_view key{data + off, path_len};

    // Begin txn on first write
    if (!started_) {
      writer_.begin_txn();
      started_ = true;
    }

    // Write full payload as value under extracted key
    writer_.put_raw(std::string(key), e.payload, /*commit_now=*/false);
    ++since_commit_;
    if (since_commit_ >= batch_size_) {
      writer_.commit_txn();
      since_commit_ = 0;
      // Start a new txn for subsequent writes
      writer_.begin_txn();
    }
  }

  void flush() override {
    std::lock_guard<std::mutex> lk(mu_);
    flush_locked();
  }

  void close() override {
    std::lock_guard<std::mutex> lk(mu_);
    flush_locked();
  }

private:
  void flush_locked() {
    if (started_) {
      writer_.commit_txn();
      since_commit_ = 0;
      started_ = false; // No open txn
    }
  }

  lahuta::LMDBWriter writer_;
  std::size_t batch_size_ = 1024;
  std::size_t since_commit_ = 0;
  bool started_ = false;
  std::mutex mu_;
};

} // namespace lahuta::pipeline::dynamic

#endif // LAHUTA_PIPELINE_DYNAMIC_SINK_LMDB_HPP
