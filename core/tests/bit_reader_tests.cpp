#include <gtest/gtest.h>

#include <cstddef>
#include <vector>

#include "md/detail/bit_reader.hpp"

using lahuta::md::detail::BitReader;

static std::vector<std::byte> B(std::initializer_list<unsigned> v) {
  std::vector<std::byte> out;
  out.reserve(v.size());
  for (auto x : v) {
    out.push_back(static_cast<std::byte>(x & 0xFFu));
  }
  return out;
}

TEST(BitReader, BulkRefillCountsOnlyRealBytes) {
  // 5 real bytes, tail padding is present but must not be counted as bits.
  const auto payload = B({0xFF, 0x00, 0xAA, 0x55, 0x0F});
  BitReader br(payload.data(), payload.size(), /*tail_zeros=*/8);

  // Read exactly 5 * 8 bits
  for (int i = 0; i < 5; ++i) {
    const auto v = br.read_bits(8);
    (void)v;
  }

  // Any additional bit must throw "Unexpected end of XTC bitstream"
  EXPECT_THROW({ (void)br.read_bits(1); }, std::exception);

  // Introspection reflects logical payload only
  EXPECT_EQ(br.bytes_consumed(), payload.size());
  EXPECT_TRUE(br.at_end());
}

TEST(BitReader, TailBoundaryResidualBitsZeroAndNonZero) {
  // One byte: 0xA0 = 1010 0000
  {
    const auto p = B({0xA0});
    BitReader br(p.data(), p.size(), 8);
    // Read 3 bits, remaining top 5 bits are 00000 , residual zero
    EXPECT_EQ(br.read_bits(3), 0b101u);
    EXPECT_TRUE(br.residual_bits_zero());
  }
  // One byte: 0xAF = 1010 1111
  {
    const auto p = B({0xAF});
    BitReader br(p.data(), p.size(), 8);
    EXPECT_EQ(br.read_bits(3), 0b101u);
    // Remaining top 5 bits are 01111
    EXPECT_FALSE(br.residual_bits_zero());
  }
}

TEST(BitReader, AtEndAndBytesConsumedClampToPayload) {
  const auto p = B({0x12, 0x34, 0x56});
  BitReader br(p.data(), p.size(), 8);

  // Read all 24 bits
  EXPECT_EQ(br.read_bits(12), 0x123u);
  EXPECT_EQ(br.read_bits(12), 0x456u);

  EXPECT_EQ(br.bytes_consumed(), p.size());
  EXPECT_TRUE(br.at_end());
  EXPECT_THROW({ (void)br.read_bits(1); }, std::exception);
}
