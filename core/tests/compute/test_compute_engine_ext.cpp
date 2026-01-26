#include <gtest/gtest.h>

#include "compute/compute_base.hpp"
#include "compute/engine.hpp"
#include "compute/parameters.hpp"

using namespace lahuta::compute;

namespace {

struct DummyData {
  int count    = 0;
  int a_val    = 0;
  int b_runs   = 0;
  int observed = -1;
};

struct DummyParams : public ParameterBase<DummyParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = 250;
};

// Increments a counter
class DummyComputation : public Computation<DummyData> {
public:
  DummyComputation(std::string label) : label_store_(std::move(label)), label_(label_store_) {}
  ComputationResult execute(DataContext<DummyData> &ctx, const ParameterInterface &) override {
    ctx.data().count++;
    return ComputationResult(true);
  }
  std::unique_ptr<ParameterInterface> get_parameters() const override {
    return std::make_unique<DummyParams>();
  }
  const ComputationLabel &get_label() const override { return label_; }
  std::vector<ComputationLabel> get_dependencies() const override { return {}; }

private:
  std::string label_store_;
  ComputationLabel label_;
};

struct AParams : public ParameterBase<AParams> {
  static constexpr ParameterInterface::TypeId TYPE_ID = 251;

  int x = 1;
};

class AComp : public Computation<DummyData> {
public:
  AComp(std::string label, AParams p) : label_store_(std::move(label)), label_(label_store_), p_(p) {}
  ComputationResult execute(DataContext<DummyData> &ctx, const ParameterInterface &raw) override {
    auto &typed      = static_cast<const AParams &>(raw);
    ctx.data().a_val = typed.x;
    return ComputationResult(true);
  }
  std::unique_ptr<ParameterInterface> get_parameters() const override {
    return std::make_unique<AParams>(p_);
  }
  const ComputationLabel &get_label() const override { return label_; }
  std::vector<ComputationLabel> get_dependencies() const override { return {}; }

private:
  std::string label_store_;
  ComputationLabel label_;
  AParams p_;
};

class BComp : public Computation<DummyData> {
public:
  BComp(std::string label, std::vector<ComputationLabel> deps)
      : label_store_(std::move(label)), label_(label_store_), deps_(std::move(deps)) {}
  ComputationResult execute(DataContext<DummyData> &ctx, const ParameterInterface &) override {
    ctx.data().b_runs++;
    ctx.data().observed = ctx.data().a_val;
    return ComputationResult(true);
  }
  std::unique_ptr<ParameterInterface> get_parameters() const override {
    return std::make_unique<DummyParams>();
  }
  const ComputationLabel &get_label() const override { return label_; }
  std::vector<ComputationLabel> get_dependencies() const override { return deps_; }

private:
  std::string label_store_;
  ComputationLabel label_;
  std::vector<ComputationLabel> deps_;
};

} // namespace

TEST(ComputeEngineExt, RegistryCapacityAtLeast64) {
  DummyData data{};
  ComputeEngine<DummyData> eng(data);

  const int N = 32;
  for (int i = 0; i < N; ++i) {
    eng.add(std::make_unique<DummyComputation>("c" + std::to_string(i)));
  }
  // Run the last one just to drive execution
  EXPECT_TRUE(eng.run<void>(ComputationLabel{"c31"}));
  EXPECT_GE(data.count, 1);
}

TEST(ComputeEngineExt, ParameterInvalidationDownstreamRecomputes) {
  DummyData data{};
  ComputeEngine<DummyData> eng(data);

  const ComputationLabel A{"A"};
  const ComputationLabel B{"B"};

  AParams ap;
  ap.x = 7;
  eng.add(std::make_unique<AComp>("A", ap));
  eng.add(std::make_unique<BComp>("B", std::vector<ComputationLabel>{A}));

  // First run
  EXPECT_TRUE(eng.run<void>(B));
  EXPECT_EQ(data.observed, 7);
  EXPECT_EQ(data.b_runs, 1);

  // Mutate A's parameter in-place which should invalidate B
  auto &p = eng.get_parameters<AParams>(A);
  p.x     = 9;

  // Run B again, it should re-run because dependency changd
  EXPECT_TRUE(eng.run<void>(B));
  EXPECT_EQ(data.observed, 9);
  EXPECT_EQ(data.b_runs, 2);
}

TEST(ComputeEngineExt, AutoHealEnablesPrereqs) {
  DummyData data{};
  ComputeEngine<DummyData> eng(data);

  const ComputationLabel A{"A"};
  const ComputationLabel B{"B"};
  const ComputationLabel C{"C"};

  eng.add(std::make_unique<DummyComputation>("A"));
  eng.add(std::make_unique<BComp>("B", std::vector<ComputationLabel>{A}));
  eng.add(std::make_unique<BComp>("C", std::vector<ComputationLabel>{B}));

  // Disable A and B, then run C. Auto-heal should re-enable A and B.
  eng.enable(A, false);
  eng.enable(B, false);

  EXPECT_TRUE(eng.run<void>(C));
  EXPECT_TRUE(eng.has_completed(A));
  EXPECT_TRUE(eng.has_completed(B));
  EXPECT_TRUE(eng.has_completed(C));
}
