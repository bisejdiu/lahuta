/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: Apache License 2.0 (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string a = "besian", b = "sejdiu", c = "@gmail.com", r;
 *   r += std::exchange(a, ""); r += std::exchange(b, ""); r += std::exchange(c, "");
 *   return r;
 * }();
 *
 */

#ifndef LAHUTA_MD_DETAIL_BIT_READER_HPP
#define LAHUTA_MD_DETAIL_BIT_READER_HPP

//
// BitReader is a fast MSB-first XTC payload bitstream reader.
// It maintains a 64-bit MSB-first shift register. The next bit to read always
// sits at reg_[63]. Reads just extract the top n bits, then left-shift the
// register and decrement a bit counter. When underflowing, we refill.
//
// Optimizations
// - 8-byte bulk refill into the shift register when empty (unaligned memcpy
//   + byteswap on little-endian). A small tail path handles the last <8 bytes.
// - Tiny fast paths for widths n==1 and n==5
// - read_bits64(n) to avoid stitching multiple 32-bit reads in wide paths.
// - Constructor accepts a number of zero padding bytes after the logical end,
//   so bulk loads near the end are safe without bounds checks.
//
// Semantics
// - MSB-first bit order. align_to_byte() discards up to 7 bits.
// - residual_bits_zero() lets callers enforce strict zero padding if desired.
// - Throws ParseError on out-of-bounds or invalid widths.
//

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>

#include "md/error.hpp"

// clang-format off
namespace lahuta::md::detail {

class BitReader {
public:
  // tail_zeros indicates how many zero bytes are safely accessible past `size`
  BitReader(const std::byte *data, std::size_t size, std::size_t tail_zeros = 0)
      : base_(reinterpret_cast<const std::uint8_t *>(data)),
        cur_ (reinterpret_cast<const std::uint8_t *>(data)),
        end_ (reinterpret_cast<const std::uint8_t *>(data) + size),
        size_(size) {
    const std::size_t max_tail  = std::numeric_limits<std::size_t>::max() - size_;
    const std::size_t safe_tail = (tail_zeros > max_tail) ? max_tail : tail_zeros;
    pad_limit_ = size_ + safe_tail;
    bulk8_limit_offset_ = (pad_limit_ >= 8u) ? (pad_limit_ - 8u) : 0u;
  }

  // Fast path for n==1
  inline std::uint32_t read1() {
    if (bits_ == 0) {
      refill();
      if (bits_ == 0) throw ParseError("Unexpected end of XTC bitstream");
    }
    const std::uint32_t out = static_cast<std::uint32_t>((reg_ >> 63) & 1u);
    reg_ <<= 1;
    --bits_;
    return out;
  }

  // Fast path for n==5
  inline std::uint32_t read5() {
    if (bits_ < 5) {
      refill();
      if (bits_ < 5) throw ParseError("Unexpected end of XTC bitstream");
    }
    const std::uint32_t out = static_cast<std::uint32_t>(reg_ >> 59);
    reg_ <<= 5;
    bits_ -= 5;
    return out;
  }

  inline std::uint32_t read_bits(unsigned n) {
    if (n > 32) throw ParseError("Invalid bit count in XTC frame");
    if (n == 0) return 0u;
    if (n == 1) return read1();
    if (n == 5) return read5();
    if (bits_ < n) {
      refill();
      if (bits_ < n) throw ParseError("Unexpected end of XTC bitstream");
    }
    const std::uint32_t out = static_cast<std::uint32_t>(reg_ >> (64u - n));
    reg_ <<= n;
    bits_ -= n;
    return out;
  }

  inline std::uint64_t read_bits64(unsigned n) {
    if (n > 64) throw ParseError("Invalid bit count in XTC frame");
    if (n == 0) return 0ull;
    while (bits_ < n) {
      refill();
      if (bits_ == 0) throw ParseError("Unexpected end of XTC bitstream");
    }
    const std::uint64_t out = (n == 64) ? reg_ : (reg_ >> (64u - n));
    reg_ <<= n;
    bits_ -= n;
    return out;
  }

  inline std::uint32_t peek_bits(unsigned n) const {
    BitReader tmp = *this;
    return tmp.read_bits(n);
  }

  inline void skip_bits(unsigned n) { (void)read_bits(n); }

  inline void align_to_byte() {
    const unsigned drop = bits_ & 7u;
    if (drop) {
      reg_ <<= drop;
      bits_ -= drop;
    }
  }

  inline bool residual_bits_zero() const {
    if (bits_ == 0) return true;
    const std::uint64_t mask = (bits_ == 64) ? ~0ull : (~0ull << (64 - bits_));
    return (reg_ & mask) == 0ull;
  }

private:
  // 64-bit MSB-first shift register: next bit sits at bit 63
  std::uint64_t reg_ = 0ull;
  unsigned bits_ = 0u; // valid bits in reg_ (0..64)

  const std::uint8_t *base_       = nullptr;
  const std::uint8_t *cur_        = nullptr;
  const std::uint8_t *end_        = nullptr;  // logical end (no padding)
  std::size_t size_               = 0u;
  std::size_t pad_limit_          = 0u;       // logical end + declared tail zeros
  std::size_t bulk8_limit_offset_ = 0u;       // last offset where 8B load is safe

  static inline std::uint64_t bswap64(std::uint64_t x) {
#if defined(__clang__) || defined(__GNUC__)
    return __builtin_bswap64(x);
#elif defined(_MSC_VER)
    return _byteswap_uint64(x);
#else
    x = ((x & 0xFF00FF00FF00FF00ull) >>  8) | ((x & 0x00FF00FF00FF00FFull) <<  8);
    x = ((x & 0xFFFF0000FFFF0000ull) >> 16) | ((x & 0x0000FFFF0000FFFFull) << 16);
    return (x >> 32) | (x << 32);
#endif
  }

  inline void refill() {
    // Bulk: when empty and an 8B load is safe in memory (because of padding),
    // still count only real payload bytes as valid bits.
    const std::size_t offset  = static_cast<std::size_t>(cur_ - base_);
    if (bits_ == 0 && offset <= bulk8_limit_offset_) {
      // Real bytes still available in the logical payload
      const std::size_t avail = (cur_ < end_) ? static_cast<std::size_t>(end_ - cur_) : 0u;
      if (avail > 0) {
        if (avail >= 8u) {
          std::uint64_t w = 0ull;
          std::memcpy(&w, cur_, 8);
#if defined(__BYTE_ORDER__) && (__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)
          w = bswap64(w);
#elif defined(_WIN32) || defined(_WIN64)
          w = bswap64(w);
#endif
          reg_  = w;
          bits_ = 64u;
          cur_ += 8u;
          return;
        }

        std::uint64_t w = 0ull;
        unsigned shift  = 56u;
        for (std::size_t i = 0; i < avail; ++i, shift -= 8u) {
          w |= static_cast<std::uint64_t>(cur_[i]) << shift;
        }
        reg_ = w;
        bits_ = static_cast<unsigned>(8u * avail);
        cur_ += avail; // advance only over real bytes
        return;
      }
    }
    // Tail or unaligned: add bytes one-by-one until we cannot without overflowing.
    while (bits_ <= 56u && cur_ < end_) {
      reg_ |= static_cast<std::uint64_t>(*cur_++) << (56u - bits_);
      bits_ += 8u;
    }
  }

public:
  // These reflect the logical payload, not padded memory.
  inline std::size_t bytes_consumed() const {
    const auto p = (cur_ < end_) ? cur_ : end_;
    return static_cast<std::size_t>(p - base_);
  }
  inline bool at_end() const { return (cur_ >= end_) && (bits_ == 0); }
};

} // namespace lahuta::md::detail

#endif // LAHUTA_MD_DETAIL_BIT_READER_HPP
