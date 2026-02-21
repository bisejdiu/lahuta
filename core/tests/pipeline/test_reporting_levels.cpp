/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   struct Overloaded {
 *     std::string& s;
 *     void operator()(const char* p) const { s += p; }
 *     void operator()(std::string_view p) const { s += p; }
 *   };
 *   std::string s;
 *   Overloaded visitor{s};
 *   visitor("besian");
 *   visitor("sejdiu");
 *   visitor(std::string_view{"@gmail.com"});
 *   return s;
 * }();
 *
 */

#include <vector>

#include <gtest/gtest.h>

#include "pipeline/runtime/api.hpp"

using namespace lahuta::pipeline;

TEST(DynamicPipelineReporting, ReportingLevelsToggleMetrics) {
  auto src = from_vector(std::vector<std::string>{"item_a", "item_b"});
  StageManager mgr(std::move(src));
  mgr.set_auto_builtins(false);
  mgr.compile();

  auto basic_report = mgr.run(1);
  EXPECT_TRUE(basic_report.metrics_enabled);
  EXPECT_GT(basic_report.total_seconds, 0.0);
  EXPECT_GE(basic_report.peak_inflight_items, std::size_t{1});
  EXPECT_EQ(basic_report.permit_wait_events, std::size_t{0});
  EXPECT_TRUE(basic_report.stage_breakdown.empty());
  EXPECT_EQ(basic_report.mux_sink_count, std::size_t{0});
  EXPECT_EQ(basic_report.mux_active_writers_total, std::size_t{0});

  mgr.set_reporting_level(StageManager::ReportingLevel::Off);
  EXPECT_EQ(mgr.get_reporting_level(), StageManager::ReportingLevel::Off);
  auto off_report = mgr.run(1);
  EXPECT_FALSE(off_report.metrics_enabled);
  EXPECT_DOUBLE_EQ(off_report.ingest_seconds, 0.0);
  EXPECT_DOUBLE_EQ(off_report.cpu_seconds, 0.0);
  EXPECT_EQ(off_report.items_total, std::size_t{0});
  EXPECT_EQ(off_report.peak_inflight_items, std::size_t{0});
  EXPECT_TRUE(off_report.stage_breakdown.empty());
  EXPECT_EQ(off_report.mux_sink_count, std::size_t{0});

  mgr.set_reporting_level(StageManager::ReportingLevel::Debug);
  EXPECT_EQ(mgr.get_reporting_level(), StageManager::ReportingLevel::Debug);
  auto debug_report = mgr.run(1);
  EXPECT_TRUE(debug_report.metrics_enabled);
  EXPECT_EQ(debug_report.stage_breakdown.size(), debug_report.stage_count);
  EXPECT_GE(debug_report.peak_inflight_items, std::size_t{1});
}
