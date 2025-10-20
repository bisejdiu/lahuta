#include <regex>
#include <string>

#include <gtest/gtest.h>

#include "version.hpp"

namespace {

constexpr auto VersionPattern = R"(^(\d+)\.(\d+)\.(\d+)(.*)$)";

} // namespace

TEST(VersionMetadata, ComponentsAreSelfConsistent) {
  const std::string version_string{lahuta::version};
  ASSERT_FALSE(version_string.empty()) << "Lahuta version string must not be empty";

  std::cmatch match;
  ASSERT_TRUE(std::regex_match(version_string.c_str(), match, std::regex(VersionPattern)))
      << "Version string '" << version_string << "' does not match expected pattern";

  const int parsed_major = std::stoi(match[1].str());
  const int parsed_minor = std::stoi(match[2].str());
  const int parsed_patch = std::stoi(match[3].str());
  const std::string parsed_suffix = match[4].str();

  EXPECT_EQ(parsed_major, lahuta::version_major);
  EXPECT_EQ(parsed_minor, lahuta::version_minor);
  EXPECT_EQ(parsed_patch, lahuta::version_patch);
  EXPECT_EQ(parsed_suffix, std::string(lahuta::version_suffix));
  EXPECT_EQ(parsed_suffix.empty(), !lahuta::version_is_prerelease);
}
