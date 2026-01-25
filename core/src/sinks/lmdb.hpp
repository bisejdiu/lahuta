#ifndef LAHUTA_PIPELINE_DYNAMIC_SINK_LMDB_HPP
#define LAHUTA_PIPELINE_DYNAMIC_SINK_LMDB_HPP

#include <cstring>
#include <memory>
#include <mutex>
#include <string>
#include <string_view>
#include <utility>

#include "db/db.hpp"
#include "db/writer.hpp"
#include "pipeline/dynamic/sink_iface.hpp"
#include "pipeline/dynamic/types.hpp"

// clang-format off
namespace lahuta::pipeline::dynamic {

//
// LmdbSink: Writes emissions to an LMDB database.
//
// - Payload is a binary serialization of analysis::system::ModelRecord.
// - Sink deserializes and reserializes into the LMDB value using MDB_RESERVE
//   so the final value address is known at serialization time. This allows the
//   ModelPayloadHeader (magic, version, field slices) to be written correctly.
//
class LmdbSink : public IDynamicSink {
public:
  explicit LmdbSink(std::shared_ptr<lahuta::LMDBDatabase> db, std::size_t batch_size)
    : db_(std::move(db)), writer_(db_->get_env(), db_->get_dbi()), batch_size_(batch_size) {
    if (batch_size_ == 0) batch_size_ = 1;
  }

  void write(EmissionView e) override {
    using analysis::system::ModelRecord;
    std::lock_guard<std::mutex> lk(mu_);
    const auto rec = serialization::Serializer<fmt::binary, ModelRecord>::deserialize(e.payload.data(), e.payload.size());
    const std::string& key = rec.file_path;

    // Begin txn on first write
    if (!started_) {
      writer_.begin_txn();
      started_ = true;
    }

    writer_.put_model_record(key, rec, /*commit_now=*/false);
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

  std::shared_ptr<lahuta::LMDBDatabase> db_;
  lahuta::LMDBWriter writer_;
  std::size_t batch_size_ = 1024;
  std::size_t since_commit_ = 0;
  bool started_ = false;
  std::mutex mu_;
};

} // namespace lahuta::pipeline::dynamic

#endif // LAHUTA_PIPELINE_DYNAMIC_SINK_LMDB_HPP
