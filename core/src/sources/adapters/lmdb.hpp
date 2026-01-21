#ifndef LAHUTA_SOURCES_ADAPTERS_LMDB_HPP
#define LAHUTA_SOURCES_ADAPTERS_LMDB_HPP

#include <memory>
#include <optional>
#include <string>

#include "sources/descriptor.hpp"
#include "sources/lmdb_descriptor.hpp"

// clang-format off
namespace lahuta::sources {

class LMDBAdapter final : public IDescriptor {
public:
  LMDBAdapter(std::string env_path, std::string db_name, std::size_t batch)
      : inner_(std::move(env_path), std::move(db_name), batch) {}

  LMDBAdapter(std::string env_path, std::string db_name, std::size_t batch, LMDBEnvOptions options)
      : inner_(std::move(env_path), std::move(db_name), batch, std::move(options)) {}

  LMDBAdapter(std::shared_ptr<LMDBDatabase> db, std::string db_name, std::size_t batch)
      : inner_(std::move(db), std::move(db_name), batch) {}

  std::optional<IngestDescriptor> next() override { return inner_.next(); }
  void reset() override { inner_.reset(); }
  std::optional<std::size_t> max_readers() const { return inner_.max_readers(); }

private:
  sources::LMDBDescriptor inner_;
};

} // namespace lahuta::sources

#endif // LAHUTA_SOURCES_ADAPTERS_LMDB_HPP
