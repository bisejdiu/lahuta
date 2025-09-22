#pragma once
//
// XTC small-int decoding utilities.
//
// Implements the core "mixed-radix" unpacking used by GROMACS XTC compression.
// Two decoding paths exist:
// - General path: decode_small_ints(reader, num_ints, num_bits, sizes, dest)
//   Works for any arity and size vector. Uses libdivide for fast integer
//   division, wide accumulators when available, and a byte-array fallback.
// - Special paths tuned for the hot cases in XTC:
//   - decode_small3_precomp: 3 integers with identical divisor d = MagicInts[i]
//     (the run-encoding fast path). Uses a single precomputed DivSpec - unrolled.
//   - decode_mixed3_precomp: 3 integers with possibly distinct divisors
//     (the first triple when bitsize>0). Uses three precomputed DivSpecs - unrolled.
//
// Support code:
// - DivSpec encapsulates divisor metadata (pow2 shift, libdivide contexts).
// - build_magicints_divspec precomputes DivSpecs for the whole MagicInts table
//   once per frame so inner loops don't pay for setup.
// - read_bits_wide stitches wide reads with BitReader::read_bits64 for speed.
//
// Fallbacks to the general path.
//

#include <array>
#include <cstdint>
#include <type_traits>

#include <libdivide/libdivide.h>

#include "md/detail/bit_reader.hpp"
#include "md/error.hpp"
#include "span.hpp"

// clang-format off
namespace lahuta::md::detail {

inline constexpr int FirstSmallIndex = 9;
inline constexpr std::array<int, 73> MagicInts = {
    0,       0,       0,       0,       0,        0,        0,       0,       0,       8,       10,
    12,      16,      20,      25,      32,       40,       50,      64,      80,      101,     128,
    161,     203,     256,     322,     406,      512,      645,     812,     1024,    1290,    1625,
    2048,    2580,    3250,    4096,    5060,     6501,     8192,    10321,   13003,   16384,   20642,
    26007,   32768,   41285,   52015,   65536,    82570,    104031,  131072,  165140,  208063,  262144,
    330280,  416127,  524287,  660561,  832255,   1048576,  1321122, 1664510, 2097152, 2642245, 3329021,
    4194304, 5284491, 6658042, 8388607, 10568983, 13316085, 16777216};

struct DivSpec {
  std::uint32_t d;
  bool is_pow2;
  unsigned shift; // valid if is_pow2
  libdivide::libdivide_u32_t fast32;
  libdivide::libdivide_u64_t fast64;
};

inline DivSpec make_divspec(std::uint32_t d) {
  bool pow2 = false;
  unsigned sh = 0u;
  if (d == 0u) {
    return DivSpec{0u, false, 0u, libdivide::libdivide_u32_t{}, libdivide::libdivide_u64_t{}};
  }
  pow2 = (d & (d - 1u)) == 0u;
#if defined(__GNUC__) || defined(__clang__)
  if (pow2) sh = static_cast<unsigned>(__builtin_ctz(d));
#else
  if (pow2) {
    while ((1u << sh) != d)
      ++sh;
  }
#endif
  return DivSpec{
      d,
      pow2,
      sh,
      libdivide::libdivide_u32_gen(d),
      libdivide::libdivide_u64_gen(static_cast<std::uint64_t>(d))};
}

inline void build_magicints_divspec(std::array<DivSpec, MagicInts.size()> &out) {
  for (std::size_t i = 0; i < MagicInts.size(); ++i) {
    const std::uint32_t d = static_cast<std::uint32_t>(MagicInts[i]);
    out[i] = make_divspec(d);
  }
}

// Specialized fast-path for 3 ints with identical base d = MagicInts[smallidx]
// bits is the total number of bits to read (= smallidx). dest must point to 3 ints.
inline void decode_small3_precomp(BitReader &reader, unsigned bits, const DivSpec &dv, int *dest);
// Specialized fast-path for 3 ints with distinct but fixed divisors (sizeint[0..2])
inline void decode_mixed3_precomp(BitReader &reader, unsigned bits, const DivSpec *divs, int *dest);

// Read up to 128 bits by concatenating <=32-bit chunks from the bitstream
template <typename Wide> inline Wide read_bits_wide(BitReader &reader, unsigned num_bits) {
  if (num_bits == 0) return static_cast<Wide>(0);
#if defined(__SIZEOF_INT128__)
  using u128 = unsigned __int128;
#endif
  if constexpr (std::is_same<Wide, std::uint64_t>::value) {
    return reader.read_bits64(num_bits);
  }
#if defined(__SIZEOF_INT128__)
  else if constexpr (std::is_same<Wide, u128>::value) {
    if (num_bits <= 64u) {
      return static_cast<u128>(reader.read_bits64(num_bits));
    } else {
      const unsigned high_bits = num_bits - 64u;
      const u128 hi = static_cast<u128>(reader.read_bits64(high_bits));
      const u128 lo = static_cast<u128>(reader.read_bits64(64u));
      return (hi << 64) | lo;
    }
  }
#endif
  else {
    // Fallback generic path
    Wide acc = 0;
    unsigned remaining = num_bits;
    while (remaining > 0) {
      const unsigned chunk = remaining > 32u ? 32u : remaining;
      const std::uint32_t part = reader.read_bits(chunk);
      acc = static_cast<Wide>((acc << chunk) | static_cast<Wide>(part));
      remaining -= chunk;
    }
    return acc;
  }
}

inline constexpr int calculate_sizeof_int(int size) {
  if (size <= 1) return 0; // only one value (0)
  int bits = 0;
  int threshold = 1;
  while (size >= threshold && bits < 32) {
    ++bits;
    threshold <<= 1;
  }
  return bits;
}

inline int calculate_sizeof_ints(span<const unsigned int> sizes) {
  int bytes[32] = {0};
  unsigned int num_bytes = 1;
  unsigned int num_bits  = 0;
  bytes[0] = 1;

  for (std::size_t idx = 0; idx < sizes.size(); ++idx) {
    unsigned int carry   = 0;
    unsigned int bytecnt = 0;
    for (; bytecnt < num_bytes; ++bytecnt) {
      unsigned int tmp = bytes[bytecnt] * sizes[idx] + carry;
      bytes[bytecnt] = static_cast<int>(tmp & 0xFFu);
      carry = tmp >> 8;
    }
    while (carry != 0) {
      if (bytecnt >= 32) throw ParseError("Internal overflow in calculate_sizeof_ints");
      bytes[bytecnt++] = static_cast<int>(carry & 0xFFu);
      carry >>= 8;
    }
    num_bytes = bytecnt;
  }

  unsigned int threshold = 1;
  if (num_bytes == 0) return 0;
  int top = num_bytes - 1;
  while (bytes[top] >= static_cast<int>(threshold)) {
    ++num_bits;
    threshold <<= 1;
  }
  return static_cast<int>(num_bits + top * 8);
}

inline void decode_small_ints(BitReader &reader, int num_ints, int num_bits, span<const unsigned int> sizes, span<int> dest) {
  if (num_ints <= 0) return;
  if (static_cast<int>(sizes.size()) < num_ints || static_cast<int>(dest.size()) < num_ints) {
    throw ParseError("decode_small_ints span size mismatch");
  }
  if (num_bits < 0) {
    throw ParseError("Invalid bit count in XTC frame");
  }

  // Precompute divisor metadata once per call
  DivSpec divs[4];
  for (int i = 0; i < num_ints; ++i) {
    const std::uint32_t d = sizes[static_cast<std::size_t>(i)];
    if (d == 0u) throw ParseError("Corrupted XTC frame (zero divisor)");
    const bool pow2 = (d & (d - 1u)) == 0u;
    unsigned sh = 0u;
#if defined(__GNUC__) || defined(__clang__)
    if (pow2) sh = static_cast<unsigned>(__builtin_ctz(d));
#else
    if (pow2) {
      while ((1u << sh) != d)
        ++sh;
    }
#endif
    divs[i] = DivSpec{
        d,
        pow2,
        sh,
        libdivide::libdivide_u32_gen(d),
        libdivide::libdivide_u64_gen(static_cast<std::uint64_t>(d))};
  }

  // Fast path: num_bits <= 32
  if (num_bits <= 32) {
    const unsigned bits = static_cast<unsigned>(num_bits);
    std::uint32_t acc = bits ? reader.read_bits(bits) : 0u;
    for (int i = num_ints - 1; i > 0; --i) {
      const DivSpec &dv = divs[i];
      std::uint32_t q;
      std::uint32_t rem;
      if (dv.d == 1u) {
        q = acc;
        rem = 0u;
      } else if (dv.is_pow2) {
        q = acc >> dv.shift;
        rem = acc & (dv.d - 1u);
      } else {
        q = libdivide::libdivide_u32_do(acc, &dv.fast32);
        rem = acc - q * dv.d;
      }
      dest[static_cast<std::size_t>(i)] = static_cast<int>(rem);
      acc = q;
    }
    dest[0] = static_cast<int>(acc);
    return;
  }

#if defined(__SIZEOF_INT128__)
  // Wider paths using 64/128-bit accumulators remove the need for the byte-array fallback
  if (num_bits <= 64) {
    std::uint64_t acc = read_bits_wide<std::uint64_t>(reader, static_cast<unsigned>(num_bits));
    for (int i = num_ints - 1; i > 0; --i) {
      const DivSpec &dv = divs[i];
      std::uint64_t q;
      std::uint64_t rem;
      if (dv.d == 1u) {
        q = acc;
        rem = 0ull;
      } else if (dv.is_pow2) {
        q = acc >> dv.shift;
        rem = acc & static_cast<std::uint64_t>(dv.d - 1u);
      } else {
        q = libdivide::libdivide_u64_do(acc, &dv.fast64);
        rem = acc - q * static_cast<std::uint64_t>(dv.d);
      }
      dest[static_cast<std::size_t>(i)] = static_cast<int>(rem);
      acc = q;
    }
    dest[0] = static_cast<int>(acc);
    return;
  } else if (num_bits <= 128) {
    using u128 = unsigned __int128;
    u128 acc = read_bits_wide<u128>(reader, static_cast<unsigned>(num_bits));
    for (int i = num_ints - 1; i > 0; --i) {
      const DivSpec &dv = divs[i];
      u128 q;
      u128 rem;
      if (dv.d == 1u) {
        q = acc;
        rem = static_cast<u128>(0);
      } else if (dv.is_pow2) {
        q = acc >> dv.shift;
        rem = acc & static_cast<u128>(dv.d - 1u);
      } else {
        q = acc / static_cast<u128>(dv.d);
        rem = acc - q * static_cast<u128>(dv.d);
      }
      dest[static_cast<std::size_t>(i)] = static_cast<int>(rem);
      acc = q;
    }
    dest[0] = static_cast<int>(acc);
    return;
  }
#endif

  // Fallback: legacy byte-array long division
  int bytes[32] = {0};
  int num_bytes = 0;
  int bits_remaining = num_bits;
  while (bits_remaining > 8) {
    if (num_bytes >= 32) throw ParseError("decode_small_ints byte buffer overflow");
    bytes[num_bytes++] = reader.read_bits(8);
    bits_remaining -= 8;
  }
  if (bits_remaining > 0) {
    if (num_bytes >= 32) throw ParseError("decode_small_ints byte buffer overflow");
    bytes[num_bytes++] = reader.read_bits(bits_remaining);
  }
  for (int idx = num_ints - 1; idx > 0; --idx) {
    if (sizes[static_cast<std::size_t>(idx)] == 0) {
      throw ParseError("Corrupted XTC frame (zero divisor)");
    }
    int value = 0;
    for (int byte_idx = num_bytes - 1; byte_idx >= 0; --byte_idx) {
      value = (value << 8) | bytes[byte_idx];
      int quotient = value / static_cast<int>(sizes[static_cast<std::size_t>(idx)]);
      bytes[byte_idx] = quotient;
      value -= quotient * static_cast<int>(sizes[static_cast<std::size_t>(idx)]);
    }
    dest[static_cast<std::size_t>(idx)] = value;
  }
  dest[0] = bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24);
}

// Definition of the specialized 3-int path
inline void decode_small3_precomp(BitReader &reader, unsigned bits, const DivSpec &dv, int *dest) {
  if (dv.d == 0u) throw ParseError("Corrupted XTC frame (zero divisor)");
  if (bits == 0u) {
    dest[0] = dest[1] = dest[2] = 0;
    return;
  }

  if (bits <= 32u) {
    std::uint32_t acc = reader.read_bits(bits);
    if (dv.d == 1u) {
      dest[2] = 0;                     // r2
      dest[1] = 0;                     // r1
      dest[0] = static_cast<int>(acc); // q1
      return;
    }
    if (dv.is_pow2) {
      const std::uint32_t m = dv.d - 1u;
      const std::uint32_t q2 = acc >> dv.shift, r2 = acc & m;
      dest[2] = static_cast<int>(r2);
      const std::uint32_t r1 = q2 & m, q1 = q2 >> dv.shift;
      dest[1] = static_cast<int>(r1);
      dest[0] = static_cast<int>(q1);
      return;
    }
    // general
    const std::uint32_t q2 = libdivide::libdivide_u32_do(acc, &dv.fast32);
    const std::uint32_t r2 = acc - q2 * dv.d;
    dest[2] = static_cast<int>(r2);
    const std::uint32_t q1 = libdivide::libdivide_u32_do(q2, &dv.fast32);
    const std::uint32_t r1 = q2 - q1 * dv.d;
    dest[1] = static_cast<int>(r1);
    dest[0] = static_cast<int>(q1);
    return;
  }

#if defined(__SIZEOF_INT128__)
  if (bits <= 64u) {
    std::uint64_t acc = reader.read_bits64(bits);
    std::uint64_t q2, r2;
    if (dv.d == 1u) {
      q2 = acc;
      r2 = 0ull;
    } else if (dv.is_pow2) {
      q2 = acc >> dv.shift;
      r2 = acc & static_cast<std::uint64_t>(dv.d - 1u);
    } else {
      q2 = libdivide::libdivide_u64_do(acc, &dv.fast64);
      r2 = acc - q2 * static_cast<std::uint64_t>(dv.d);
    }
    dest[2] = static_cast<int>(r2);
    std::uint64_t q1, r1;
    if (dv.d == 1u) {
      q1 = q2;
      r1 = 0ull;
    } else if (dv.is_pow2) {
      r1 = q2 & static_cast<std::uint64_t>(dv.d - 1u);
      q1 = q2 >> dv.shift;
    } else {
      q1 = libdivide::libdivide_u64_do(q2, &dv.fast64);
      r1 = q2 - q1 * static_cast<std::uint64_t>(dv.d);
    }
    dest[1] = static_cast<int>(r1);
    dest[0] = static_cast<int>(q1);
    return;
  } else if (bits <= 128u) {
    using u128 = unsigned __int128;
    u128 acc = read_bits_wide<u128>(reader, bits);
    u128 q2, r2;
    if (dv.d == 1u) {
      q2 = acc;
      r2 = static_cast<u128>(0);
    } else if (dv.is_pow2) {
      q2 = acc >> dv.shift;
      r2 = acc & static_cast<u128>(dv.d - 1u);
    } else {
      q2 = acc / static_cast<u128>(dv.d);
      r2 = acc - q2 * static_cast<u128>(dv.d);
    }
    dest[2] = static_cast<int>(r2);
    u128 q1, r1;
    if (dv.d == 1u) {
      q1 = q2;
      r1 = static_cast<u128>(0);
    } else if (dv.is_pow2) {
      r1 = q2 & static_cast<u128>(dv.d - 1u);
      q1 = q2 >> dv.shift;
    } else {
      q1 = q2 / static_cast<u128>(dv.d);
      r1 = q2 - q1 * static_cast<u128>(dv.d);
    }
    dest[1] = static_cast<int>(r1);
    dest[0] = static_cast<int>(q1);
    return;
  }
#endif

  // Rare?: fallback to generic when bits > 128 or no 128-bit support
  unsigned sizes_arr[3] = {dv.d, dv.d, dv.d};
  decode_small_ints(
      reader,
      3,
      static_cast<int>(bits),
      span<const unsigned int>(sizes_arr),
      span<int>(dest, 3));
}

// Mixed divisors specialization: sizes may differ between components
inline void decode_mixed3_precomp(BitReader &reader, unsigned bits, const DivSpec *divs, int *dest) {
  const DivSpec &d2 = divs[2];
  const DivSpec &d1 = divs[1];
  if (bits == 0u) {
    dest[0] = dest[1] = dest[2] = 0;
    return;
  }

  if (bits <= 32u) {
    std::uint32_t acc = reader.read_bits(bits);
    // dest[2]
    std::uint32_t q2, r2;
    if (d2.d == 1u) {
      q2 = acc;
      r2 = 0u;
    } else if (d2.is_pow2) {
      q2 = acc >> d2.shift;
      r2 = acc & (d2.d - 1u);
    } else {
      q2 = libdivide::libdivide_u32_do(acc, &d2.fast32);
      r2 = acc - q2 * d2.d;
    }
    dest[2] = static_cast<int>(r2);
    // dest[1]
    std::uint32_t q1, r1;
    if (d1.d == 1u) {
      q1 = q2;
      r1 = 0u;
    } else if (d1.is_pow2) {
      r1 = q2 & (d1.d - 1u);
      q1 = q2 >> d1.shift;
    } else {
      q1 = libdivide::libdivide_u32_do(q2, &d1.fast32);
      r1 = q2 - q1 * d1.d;
    }
    dest[1] = static_cast<int>(r1);
    dest[0] = static_cast<int>(q1);
    return;
  }

#if defined(__SIZEOF_INT128__)
  if (bits <= 64u) {
    std::uint64_t acc = reader.read_bits64(bits);
    const DivSpec &d2l = divs[2];
    const DivSpec &d1l = divs[1];
    std::uint64_t q2, r2;
    if (d2l.d == 1u) {
      q2 = acc;
      r2 = 0ull;
    } else if (d2l.is_pow2) {
      q2 = acc >> d2l.shift;
      r2 = acc & static_cast<std::uint64_t>(d2l.d - 1u);
    } else {
      q2 = libdivide::libdivide_u64_do(acc, &d2l.fast64);
      r2 = acc - q2 * static_cast<std::uint64_t>(d2l.d);
    }
    dest[2] = static_cast<int>(r2);
    std::uint64_t q1, r1;
    if (d1l.d == 1u) {
      q1 = q2;
      r1 = 0ull;
    } else if (d1l.is_pow2) {
      r1 = q2 & static_cast<std::uint64_t>(d1l.d - 1u);
      q1 = q2 >> d1l.shift;
    } else {
      q1 = libdivide::libdivide_u64_do(q2, &d1l.fast64);
      r1 = q2 - q1 * static_cast<std::uint64_t>(d1l.d);
    }
    dest[1] = static_cast<int>(r1);
    dest[0] = static_cast<int>(q1);
    return;
  } else if (bits <= 128u) {
    using u128 = unsigned __int128;
    u128 acc = read_bits_wide<u128>(reader, bits);
    const DivSpec &d2q = divs[2];
    const DivSpec &d1q = divs[1];
    u128 q2, r2;
    if (d2q.d == 1u) {
      q2 = acc;
      r2 = static_cast<u128>(0);
    } else if (d2q.is_pow2) {
      q2 = acc >> d2q.shift;
      r2 = acc & static_cast<u128>(d2q.d - 1u);
    } else {
      q2 = acc / static_cast<u128>(d2q.d);
      r2 = acc - q2 * static_cast<u128>(d2q.d);
    }
    dest[2] = static_cast<int>(r2);
    u128 q1, r1;
    if (d1q.d == 1u) {
      q1 = q2;
      r1 = static_cast<u128>(0);
    } else if (d1q.is_pow2) {
      r1 = q2 & static_cast<u128>(d1q.d - 1u);
      q1 = q2 >> d1q.shift;
    } else {
      q1 = q2 / static_cast<u128>(d1q.d);
      r1 = q2 - q1 * static_cast<u128>(d1q.d);
    }
    dest[1] = static_cast<int>(r1);
    dest[0] = static_cast<int>(q1);
    return;
  }
#endif

  // Fallback to generic if bits > 128 or no 128-bit support
  unsigned sizes_arr[3] = {divs[0].d, divs[1].d, divs[2].d};
  decode_small_ints(
      reader,
      3,
      static_cast<int>(bits),
      span<const unsigned int>(sizes_arr),
      span<int>(dest, 3));
}

} // namespace lahuta::md::detail
