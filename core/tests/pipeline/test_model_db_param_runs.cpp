#include <filesystem>
#include <memory>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

#include <gtest/gtest.h>

#include "analysis/contacts/computation.hpp"
#include "analysis/system/model_pack_task.hpp"
#include "db/db.hpp"
#include "io/sinks/lmdb.hpp"
#include "pipeline/compute/parameters.hpp"
#include "pipeline/dynamic/manager.hpp"
#include "pipeline/dynamic/sink_iface.hpp"
#include "pipeline/dynamic/sources.hpp"
#include "runtime.hpp"

using namespace lahuta;
using namespace lahuta::pipeline;

// clang-format off
namespace {
namespace fs = std::filesystem;

static fs::path make_unique_temp_directory(const std::string& prefix) {
  auto temp_dir = fs::temp_directory_path();
  std::random_device rd;
  std::mt19937_64 rng(rd());
  std::uniform_int_distribution<uint64_t> dist;

  for (int attempt = 0; attempt < 128; ++attempt) {
    std::ostringstream name;
    name << prefix << std::hex << dist(rng);
    auto candidate = temp_dir / name.str();

    std::error_code ec;
    if (fs::create_directory(candidate, ec)) return candidate;
    if (!ec) continue; // directory already existed
    if (ec != std::make_error_code(std::errc::file_exists)) {
      throw std::system_error(ec, "make_unique_temp_directory: failed to create '" + candidate.string() + "'");
    }
  }
  throw std::runtime_error("make_unique_temp_directory: exhausted attempts to create unique directory");
}

// Simple in-memory sink collecting payloads from a single "contacts" channel
class CollectorSink : public dynamic::IDynamicSink {
public:
  std::vector<std::string> payloads;
  void write(dynamic::EmissionView e) override { payloads.emplace_back(e.payload); }
  void close() override {}
  void flush() override {}
};

// Very naive impl: looks for "num_contacts":<number>
static std::optional<std::size_t> extract_num_contacts(const std::string &json) {
  const std::string key = "\"num_contacts\":";
  auto pos = json.find(key);
  if (pos == std::string::npos) return std::nullopt;
  pos += key.size();
  while (pos < json.size() && std::isspace(static_cast<unsigned char>(json[pos]))) ++pos;
  std::size_t end = pos;
  while (end < json.size() && std::isdigit(static_cast<unsigned char>(json[end]))) ++end;
  if (end == pos) return std::nullopt;
  try {
    return static_cast<std::size_t>(std::stoull(json.substr(pos, end - pos)));
  } catch (...) {
    return std::nullopt;
  }
}

// Create a fresh LMDB database at db_path from all *.cif.gz files in data_dir.
static void build_database_from_directory(const fs::path &data_dir, const fs::path &db_path, int threads) {
  if (fs::exists(db_path)) fs::remove_all(db_path);

  auto db = std::make_shared<LMDBDatabase>(db_path.string());

  auto src = dynamic::sources_factory::from_directory(
      data_dir.string(), ".cif.gz", /*recursive=*/false, /*batch_size=*/50);
  dynamic::StageManager mgr(std::move(src));

  auto task = std::make_shared<analysis::system::ModelPackTask>("db");
  mgr.add_task("createdb", /*deps*/{}, task, /*thread_safe=*/true);

  mgr.connect_sink("db", std::make_shared<dynamic::LmdbSink>(db, 50));

  mgr.compile();
  mgr.run(static_cast<std::size_t>(threads));
}

// Run the compute pipeline over ALL keys in the LMDB at db_path, and collect contacts.
static std::shared_ptr<CollectorSink> run_contacts_pipeline_over_db(const fs::path &db_path, int threads) {
  auto src = dynamic::sources_factory::from_lmdb(db_path.string(), std::string{}, 64);
  dynamic::StageManager mgr(std::move(src));

  mgr.set_auto_builtins(true);
  mgr.get_system_params().is_model = true;
  mgr.get_topology_params().atom_typing_method = AtomTypingMethod::Arpeggio;

  {
    pipeline::compute::ContactsParams p{};
    p.provider = analysis::contacts::ContactProvider::Arpeggio;
    p.type = InteractionType::All;
    p.channel = "contacts";
    p.format = pipeline::compute::ContactsOutputFormat::Json;
    mgr.add_computation(
        "contacts",
        {},
        [label = std::string("contacts"), p]() {
          return std::make_unique<analysis::contacts::ContactsComputation>(label, p);
        },
        /*thread_safe=*/true);
  }

  auto contacts_collector = std::make_shared<CollectorSink>();
  mgr.connect_sink("contacts", contacts_collector);

  mgr.compile();
  mgr.run(static_cast<std::size_t>(threads));

  return contacts_collector;
}

} // namespace

// This test builds a DB from two models then runs two separate parameterized runs.
// Each run uses a fresh DB path, and asserts they complete and produce output.
TEST(ModelDatabasePipeline, MultiItemDbParameterizedRunsDoNotCrash) {
  const int threads = 1;
  LahutaRuntime::ensure_initialized(static_cast<std::size_t>(threads));

  fs::path here(__FILE__);
  fs::path core_dir = here.parent_path().parent_path().parent_path();
  fs::path data_dir = core_dir / "data" / "models";

  ASSERT_TRUE(fs::exists(data_dir)) << "Missing models directory: " << data_dir;
  ASSERT_TRUE(fs::exists(data_dir / "AF-P0CL56-F1-model_v4.cif.gz"));
  ASSERT_TRUE(fs::exists(data_dir / "AF-Q57552-F1-model_v4.cif.gz"));

  fs::path base;
  ASSERT_NO_THROW(base = make_unique_temp_directory("lahuta_test_db_param_"));

  auto cleanup = [&]() {
    std::error_code ec;
    fs::remove_all(base, ec);
  };

  std::vector<std::string> targets = {
      "AF-P0CL56-F1-model_v4.cif.gz",
      "AF-Q57552-F1-model_v4.cif.gz",
  };

  int failures = 0;

  for (std::size_t i = 0; i < targets.size(); ++i) {
    SCOPED_TRACE(testing::Message() << "param run #" << (i + 1));
    fs::path db_path = base / (std::string("run") + std::to_string(i + 1));

    ASSERT_NO_THROW(build_database_from_directory(data_dir, db_path, threads));

    std::shared_ptr<CollectorSink> sink;
    ASSERT_NO_THROW(sink = run_contacts_pipeline_over_db(db_path, threads));
    ASSERT_TRUE(sink);

    // Check we collected something and that the target model appears in at least one payload
    EXPECT_GT(sink->payloads.size(), 0u);

    bool found_target = false;
    for (const auto &payload : sink->payloads) {
      if (payload.find(targets[i]) != std::string::npos) {
        found_target = true;
        (void)extract_num_contacts(payload);
        break;
      }
    }

    if (!found_target) failures += 1;
    EXPECT_TRUE(found_target) << "Did not find target model in outputs: " << targets[i];
  }

  cleanup();

  EXPECT_EQ(failures, 0) << "One of the parameterized runs failed. Could be a regression.";
}

TEST(ModelDatabasePipeline, SingleItemDbPipelineProcessesOneModel) {
  const int threads = 1;
  LahutaRuntime::ensure_initialized(static_cast<std::size_t>(threads));

  fs::path here(__FILE__);
  fs::path core_dir = here.parent_path().parent_path().parent_path();
  fs::path data_dir = core_dir / "data" / "models";
  ASSERT_TRUE(fs::exists(data_dir));

  fs::path target = data_dir / "AF-P0CL56-F1-model_v4.cif.gz";
  ASSERT_TRUE(fs::exists(target));

  fs::path base;
  ASSERT_NO_THROW(base = make_unique_temp_directory("lahuta_test_db_single_"));
  auto cleanup = [&]() {
    std::error_code ec;
    fs::remove_all(base, ec);
  };

  ASSERT_NO_THROW(build_database_from_directory(data_dir, base, threads));

  auto src2 = dynamic::sources_factory::from_vector(std::vector<std::string>{target.string()});
  dynamic::StageManager mgr(std::move(src2));

  mgr.set_auto_builtins(true);
  mgr.get_system_params().is_model = true;
  mgr.get_topology_params().atom_typing_method = AtomTypingMethod::Arpeggio;

  {
    pipeline::compute::ContactsParams p{};
    p.provider = analysis::contacts::ContactProvider::Arpeggio;
    p.type = InteractionType::All;
    p.channel = "contacts";
    p.format = pipeline::compute::ContactsOutputFormat::Json;
    mgr.add_computation(
        "contacts",
        {},
        [label = std::string("contacts"), p]() {
          return std::make_unique<analysis::contacts::ContactsComputation>(label, p);
        },
        /*thread_safe=*/true);
  }

  auto sink = std::make_shared<CollectorSink>();
  mgr.connect_sink("contacts", sink);

  ASSERT_NO_THROW(mgr.compile());
  ASSERT_NO_THROW(mgr.run(static_cast<std::size_t>(threads)));
  ASSERT_FALSE(sink->payloads.empty());
  bool found = false;
  for (const auto &payload : sink->payloads) {
    if (payload.find(target.filename().string()) != std::string::npos) { found = true; break; }
  }
  EXPECT_TRUE(found) << "Single-item pipeline did not produce contacts for the target.";

  cleanup();
}
