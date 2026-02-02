#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>
#include <string_view>
#include <vector>

#include <Geometry/point.h>

#include "analysis/system/model_parse_task.hpp"
#include "hash/fnv1a.hpp"
#include "logging/logging.hpp"
#include "pipeline/data/data_requirements.hpp"
#include "tasks/positions_task.hpp"

namespace lahuta::cli::positions {
namespace A = lahuta::analysis;
namespace P = lahuta::pipeline;
namespace {

std::string sanitize_model_name(std::string_view model_path) {
  std::filesystem::path path(model_path);
  auto name = path.filename().string();
  if (!name.empty()) {
    return name;
  }
  return std::string(model_path);
}

std::string hex_byte(uint8_t value) {
  constexpr char Hex[] = "0123456789abcdef";
  std::string out(2, '0');
  out[0] = Hex[(value >> 4) & 0x0F];
  out[1] = Hex[value & 0x0F];
  return out;
}

std::filesystem::path output_path_for_model(const std::filesystem::path &base_dir,
                                            const std::string &model_name, int tree_depth) {
  if (tree_depth <= 0) {
    return base_dir / (model_name + ".bin");
  }

  const uint64_t hash = fnv1a_64(model_name);
  const auto level1   = hex_byte(static_cast<uint8_t>((hash >> 8) & 0xFF));
  if (tree_depth == 1) {
    return base_dir / level1 / (model_name + ".bin");
  }

  const auto level2 = hex_byte(static_cast<uint8_t>(hash & 0xFF));
  return base_dir / level1 / level2 / (model_name + ".bin");
}

void write_positions_as_float(const RDGeom::POINT3D_VECT &pts, const std::filesystem::path &path) {
  std::ofstream ofs(path, std::ios::binary);
  if (!ofs) {
    throw std::runtime_error("Failed to open file for writing: " + path.string());
  }

  std::vector<float> data;
  data.reserve(pts.size() * 3);
  for (const auto &pt : pts) {
    data.push_back(static_cast<float>(pt.x));
    data.push_back(static_cast<float>(pt.y));
    data.push_back(static_cast<float>(pt.z));
  }

  ofs.write(reinterpret_cast<const char *>(data.data()),
            static_cast<std::streamsize>(data.size() * sizeof(float)));
}

class PositionsTask final : public P::ITask {
public:
  explicit PositionsTask(std::filesystem::path output_dir, int tree_depth)
      : output_dir_(std::move(output_dir)), tree_depth_(tree_depth) {}

  P::TaskResult run(const std::string &item_path, P::TaskContext &ctx) override {
    auto payload                          = ctx.model_payload();
    const RDGeom::POINT3D_VECT *positions = nullptr;
    std::shared_ptr<const ModelParserResult> parsed;
    if (payload && payload->positions && !payload->positions->empty()) {
      positions = payload->positions.get();
    } else if (!payload) {
      parsed = A::get_parsed_model_result(ctx);
      if (parsed && parsed->coords_size() > 0) {
        positions = &parsed->coords;
      }
    }

    if (!positions || positions->empty()) {
      Logger::get_logger()->warn("[positions:input] Missing position data for '{}'", item_path);
      return {};
    }

    return write_positions_from_vector(item_path, *positions);
  }

  P::DataFieldSet data_requirements() const override {
    return P::DataFieldSet::of({P::DataField::Positions});
  }

private:
  P::TaskResult write_positions_from_vector(const std::string &item_path,
                                            const RDGeom::POINT3D_VECT &positions) {
    const std::string model_name = sanitize_model_name(item_path);
    if (model_name.empty()) {
      Logger::get_logger()->error("[positions:input] Empty model name for '{}'", item_path);
      return P::TaskResult{false, {}};
    }

    const auto output_path = output_path_for_model(output_dir_, model_name, tree_depth_);

    P::TaskResult result;
    try {
      std::filesystem::create_directories(output_path.parent_path());
      write_positions_as_float(positions, output_path);
      Logger::get_logger()->debug("[positions:write] Wrote {} atoms to {}",
                                  positions.size(),
                                  output_path.string());
      result.ok = true;
    } catch (const std::exception &e) {
      Logger::get_logger()->error("[positions:write] Failed to write '{}': {}",
                                  output_path.string(),
                                  e.what());
      result.ok = false;
    }
    return result;
  }

  std::filesystem::path output_dir_;
  int tree_depth_ = 1;
};

} // namespace

std::shared_ptr<P::ITask> make_positions_task(std::filesystem::path output_dir, int tree_depth) {
  return std::make_shared<PositionsTask>(std::move(output_dir), tree_depth);
}

} // namespace lahuta::cli::positions
