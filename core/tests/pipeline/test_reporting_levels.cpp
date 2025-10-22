#include <vector>

#include <gtest/gtest.h>

#include "pipeline/dynamic/manager.hpp"
#include "pipeline/dynamic/sources.hpp"

using namespace lahuta::pipeline::dynamic;

TEST(DynamicPipelineReporting, ReportingLevelsToggleMetrics) {
  auto src = sources_factory::from_vector(std::vector<std::string>{"item_a", "item_b"});
  StageManager mgr(std::move(src));
  mgr.set_auto_builtins(false);
  mgr.compile();

  auto basic_report = mgr.run(1);
  EXPECT_TRUE(basic_report.metrics_enabled);
  EXPECT_GT(basic_report.total_seconds, 0.0);

  mgr.set_reporting_level(StageManager::ReportingLevel::Off);
  EXPECT_EQ(mgr.get_reporting_level(), StageManager::ReportingLevel::Off);
  auto off_report = mgr.run(1);
  EXPECT_FALSE(off_report.metrics_enabled);
  EXPECT_DOUBLE_EQ(off_report.ingest_seconds, 0.0);
  EXPECT_DOUBLE_EQ(off_report.cpu_seconds, 0.0);
  EXPECT_EQ(off_report.items_total, std::size_t{0});

  mgr.set_reporting_level(StageManager::ReportingLevel::Debug);
  EXPECT_EQ(mgr.get_reporting_level(), StageManager::ReportingLevel::Debug);
  auto debug_report = mgr.run(1);
  EXPECT_TRUE(debug_report.metrics_enabled);
}
