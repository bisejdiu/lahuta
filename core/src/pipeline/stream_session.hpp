#ifndef LAHUTA_PIPELINE_STREAM_SESSION_HPP
#define LAHUTA_PIPELINE_STREAM_SESSION_HPP

#include <cstddef>
#include <functional>
#include <memory>
#include <string_view>

#include "pipeline/data_requirements.hpp"
#include "pipeline/model_payload.hpp"
#include "topology.hpp"

// clang-format off
namespace lahuta {

class Luni;
struct ModelMetadata;

struct StreamSession {
  virtual ~StreamSession() = default;
  virtual std::string_view get_session_id() const = 0;
  virtual std::shared_ptr<const Luni>     get_or_load_system() const = 0;
  virtual std::shared_ptr<const Topology> get_or_load_topology(const TopologyBuildingOptions &) const = 0;
  virtual pipeline::ModelPayloadSlices model_payload(pipeline::DataFieldSet /*requested*/) const { return {}; }

  class Permit { // RAII counting-semaphore token
  public:
    using ReleaseCallback = std::function<void()>;
    explicit Permit(ReleaseCallback releaser) : releaser_(std::move(releaser)) {}
    ~Permit() { if (releaser_) releaser_(); }
    Permit(const Permit &) = delete;
    Permit &operator=(const Permit &) = delete;
    Permit(Permit &&) noexcept = default;
    Permit &operator=(Permit &&) noexcept = default;

  private:
    ReleaseCallback releaser_;
  };

  virtual Permit acquire_permit() const = 0;
  virtual std::size_t max_inflight_frames() const noexcept = 0;
};

} // namespace lahuta

#endif // LAHUTA_PIPELINE_STREAM_SESSION_HPP
