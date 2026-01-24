
//
// XtcReader reads frame headers via XDR, then decodes coordinates either uncompressed or compressed
//
// Performance-oriented bits
// - Bit IO: a 64-bit MSB-first BitReader does bulk 8B refills and has tiny
//   fast paths for 1- and 5-bit reads. I pad the payload with 8 zero bytes
//   to make near-end unaligned loads safe.
// - Decoding: the two hottest decode steps use precomputed libdivide contexts
//   and unrolled 3-int splitters:
//   - decode_small3_precomp: runs (3 equal divisors = MagicInts[smallidx])
//   - decode_mixed3_precomp: the first triple when bitsize>0 (3 distinct divisors = sizeint[0..2])
// - Safety checks (caps, bounds, invariants) are hoisted outside hot loops
//   where possible, without sacrificing correctness.
//
// Random access helpers build a lightweight index by scanning headers, keeping
// file offsets, steps, and times. All floating-point box elements are validated
// for finiteness and plausible diagonals.
//

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <memory>
#include <optional>
#include <ostream>
#include <streambuf>
#include <utility>
#include <vector>

#include "md/detail/bit_reader.hpp"
#include "md/detail/xtc_smallints.hpp"
#include "md/xtc_reader.hpp"
#include "utils/span.hpp"

// clang-format off
namespace lahuta::md {

namespace {

constexpr std::size_t Dim     = 3;
constexpr float NmToAngstrom  = 10.0f;
constexpr int XtcMagic        = 1995;
constexpr int XtcNewMagic     = 2023;
constexpr std::uint32_t MaxNatomsOldMagic = 298261617u;

[[ nodiscard ]] inline RDGeom::Point3D make_point(const int coord[3], double scale) {
  return RDGeom::Point3D(
      static_cast<double>(coord[0]) * scale,
      static_cast<double>(coord[1]) * scale,
      static_cast<double>(coord[2]) * scale);
}

} // namespace

namespace {

inline bool checked_mul_size_t(std::size_t a, std::size_t b, std::size_t &out) {
  if (a == 0 || b == 0) { out = 0; return true; }
  if (a > (std::numeric_limits<std::size_t>::max() / b)) return false;
  out = a * b;
  return true;
}

class VectorOutputBuffer final : public std::streambuf {
public:
  explicit VectorOutputBuffer(std::vector<char> &target) : target_(target) {}

protected:
  int_type overflow(int_type ch) override {
    if (ch == traits_type::eof()) {
      return traits_type::not_eof(ch);
    }
    target_.push_back(static_cast<char>(ch));
    return ch;
  }

  std::streamsize xsputn(const char *s, std::streamsize count) override {
    if (count > 0) {
      target_.insert(target_.end(), s, s + count);
    }
    return count;
  }

private:
  std::vector<char> &target_;
};

} // namespace

XtcReader::XtcReader(std::string_view filename) : filename_(std::string(filename)) {
  auto file_stream = std::make_shared<std::ifstream>(std::string(filename), std::ios::binary);
  if (!file_stream->good()) {
    throw FileError(std::string("Cannot open XTC file: ") + std::string(filename));
  }
  stream_ptr_ = file_stream;
  stream_     = stream_ptr_.get();
  xdr_        = std::make_unique<internal::XDRReader>(stream_);
}

XtcReader::XtcReader(std::istream &stream) : stream_(&stream) {
  xdr_ = std::make_unique<internal::XDRReader>(stream_);
}

void XtcReader::initialize() {
  if (!read_next_frame()) {
    throw ParseError("Unable to read first XTC frame");
  }
}

bool XtcReader::begin_frame(FrameHeader &header, FramePacket *packet_meta, std::uint32_t &payload_atoms) {
  if (!read_frame_header(header)) {
    return false;
  }

  current_header_ = header;
  natoms_         = header.natoms;
  current_precision_ = 0.0f;

  box_ = RDGeom::Point3D(
      header.box[0] * NmToAngstrom,
      header.box[4] * NmToAngstrom,
      header.box[8] * NmToAngstrom);

  if (header.natoms > caps_.max_atoms) {
    throw ParseError("Atom count exceeds configured cap");
  }

  if (!xdr_->read(&payload_atoms)) {
    return false;
  }
  if (payload_atoms != header.natoms) {
    throw ParseError("Atom count in XTC frame does not match header");
  }

  if (packet_meta) {
    packet_meta->atom_count   = header.natoms;
    packet_meta->step         = header.step;
    packet_meta->timestamp_ps = std::nullopt;
    if (std::isfinite(static_cast<double>(header.time))) {
      packet_meta->timestamp_ps = static_cast<double>(header.time);
    }
  }

  return true;
}

void XtcReader::update_current_index() {
  if (!frame_index_ready_) return;

  auto it = std::find_if(frame_index_.begin(), frame_index_.end(), [&](const FrameIndexEntry &entry) {
    return entry.step == current_header_.step && std::fabs(entry.time - current_header_.time) < 1e-6f;
  });
  if (it != frame_index_.end()) {
    current_index_ = static_cast<std::size_t>(std::distance(frame_index_.begin(), it));
  }
}

bool XtcReader::read_next_frame() {
  FrameHeader header{};
  std::uint32_t payload_atoms = 0;
  if (!begin_frame(header, nullptr, payload_atoms)) {
    has_current_ = false;
    coords_.clear();
    return false;
  }

  coords_.clear();
  coords_.reserve(static_cast<std::size_t>(payload_atoms));

  bool ok = true;
  if (payload_atoms == 0) {
    current_precision_ = -1.0f;
  } else if (payload_atoms <= min_compressed_system_size) {
    ok = read_uncompressed_coords_into(coords_, payload_atoms);
  } else {
    ok = read_compressed_coords_into(coords_, payload_atoms);
  }

  has_current_ = ok;
  if (!ok) {
    coords_.clear();
    return false;
  }

  update_current_index();
  return true;
}

bool XtcReader::read_next_into(RDGeom::POINT3D_VECT &out) {
  FrameHeader header{};
  std::uint32_t payload_atoms = 0;
  if (!begin_frame(header, nullptr, payload_atoms)) {
    has_current_ = false;
    out.clear();
    return false;
  }

  bool ok = true;
  if (payload_atoms == 0) {
    out.clear();
    current_precision_ = -1.0f;
  } else if (payload_atoms <= min_compressed_system_size) {
    ok = read_uncompressed_coords_into(out, payload_atoms);
  } else {
    ok = read_compressed_coords_into(out, payload_atoms);
  }

  has_current_ = ok;
  if (!ok) {
    out.clear();
    return false;
  }

  update_current_index();
  return true;
}

bool XtcReader::read_next_packet(FramePacket &packet) {
  FrameHeader header{};
  std::uint32_t payload_atoms = 0;
  if (!begin_frame(header, &packet, payload_atoms)) {
    has_current_ = false;
    packet.payload.reset();
    return false;
  }

  std::vector<char> storage;
  storage.reserve(64 + static_cast<std::size_t>(header.natoms) * Dim * sizeof(float));
  VectorOutputBuffer buf(storage);
  std::ostream os(&buf);
  internal::XDRWriter writer(&os);

  auto xtc_write_scalar = [&writer](auto value) {
    if (!writer.write(value)) {
      throw ParseError("Failed to serialize XTC frame component");
    }
  };
  auto xtc_write_array = [&writer](const auto *data, std::size_t count) {
    if (!writer.write(data, count)) {
      throw ParseError("Failed to serialize XTC frame array");
    }
  };

  xtc_write_scalar(header.magic);
  xtc_write_scalar(static_cast<int>(header.natoms));
  xtc_write_scalar(static_cast<int>(static_cast<std::int32_t>(header.step)));
  xtc_write_scalar(header.time);
  xtc_write_array(header.box.data(), header.box.size());
  xtc_write_scalar(payload_atoms);

  if (payload_atoms == 0) {
    current_precision_ = -1.0f;
  } else if (payload_atoms <= min_compressed_system_size) {
    std::size_t size3 = 0;
    if (!checked_mul_size_t(static_cast<std::size_t>(payload_atoms), Dim, size3)) {
      throw ParseError("Coordinate array size overflow");
    }
    std::vector<float> raw(size3);
    if (size3 > 0 && !xdr_->read(raw.data(), raw.size())) {
      has_current_ = false;
      return false;
    }
    if (size3 > 0) {
      xtc_write_array(raw.data(), raw.size());
    }
    current_precision_ = -1.0f;
  } else {
    float precision = 0.0f;
    if (!xdr_->read(&precision)) {
      has_current_ = false;
      return false;
    }
    if (precision < 1e-12f) {
      throw ParseError("XTC precision too small");
    }
    if (!(precision > 0.0f) || !std::isfinite(precision)) {
      throw ParseError("Invalid XTC precision value");
    }
    xtc_write_scalar(precision);

    std::array<int, Dim> minint{};
    if (!xdr_->read(minint.data(), minint.size())) {
      has_current_ = false;
      return false;
    }
    xtc_write_array(minint.data(), minint.size());

    std::array<int, Dim> maxint{};
    if (!xdr_->read(maxint.data(), maxint.size())) {
      has_current_ = false;
      return false;
    }
    xtc_write_array(maxint.data(), maxint.size());

    for (std::size_t i = 0; i < Dim; ++i) {
      if (maxint[i] < minint[i]) {
        throw ParseError("XTC frame has invalid bounds");
      }
      auto span64 = static_cast<long long>(maxint[i]) - static_cast<long long>(minint[i]) + 1LL;
      if (span64 <= 0) {
        throw ParseError("XTC frame has invalid bounds");
      }
      if (span64 > std::numeric_limits<unsigned int>::max()) {
        throw ParseError("XTC span too large");
      }
    }

    int smallidx = 0;
    if (!xdr_->read(&smallidx)) {
      has_current_ = false;
      return false;
    }
    if (smallidx < detail::FirstSmallIndex || smallidx >= static_cast<int>(detail::MagicInts.size())) {
      throw ParseError("Invalid small index in XTC frame");
    }
    xtc_write_scalar(smallidx);

    std::uint64_t byte_count64 = 0;
    if (current_header_.magic == XtcNewMagic) {
      if (!xdr_->read(&byte_count64)) {
        has_current_ = false;
        return false;
      }
      xtc_write_scalar(byte_count64);
    } else {
      unsigned int byte_count32 = 0;
      if (!xdr_->read(&byte_count32)) {
        has_current_ = false;
        return false;
      }
      xtc_write_scalar(byte_count32);
      byte_count64 = byte_count32;
    }

    if (byte_count64 > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
      throw ParseError("XTC frame bitstream is too large to load");
    }
    std::size_t byte_count = static_cast<std::size_t>(byte_count64);

    std::uint64_t cap64 = static_cast<std::uint64_t>(header.natoms)
        * static_cast<std::uint64_t>(caps_.max_bitstream_bytes_per_atom) + 64ull;
    if (byte_count64 > cap64) {
      throw ParseError("Unreasonable XTC payload size for atom count");
    }

    std::vector<char> payload(byte_count);
    if (byte_count > 0 && !xdr_->read(payload.data(), payload.size())) {
      has_current_ = false;
      return false;
    }
    if (!payload.empty()) {
      xtc_write_array(payload.data(), payload.size());
    }

    current_precision_ = precision;
  }

  has_current_ = true;
  update_current_index();

  auto owned = std::make_shared<std::vector<char>>(std::move(storage));
  packet.payload = std::const_pointer_cast<const std::vector<char>>(owned);
  return true;
}

std::vector<RDGeom::Point3D> XtcReader::move_current_coords() {
  if (!has_current_) {
    throw ParseError("No current XTC frame available");
  }
  return std::move(coords_);
}

bool XtcReader::read_frame_header(FrameHeader &header) {
  int magic = 0;
  if (!xdr_->read(&magic)) {
    return false;
  }
  if (magic != XtcMagic && magic != XtcNewMagic) {
    throw ParseError("Invalid XTC magic number");
  }
  header.magic = magic;

  int natoms = 0;
  if (!xdr_->read(&natoms)) return false;
  if (natoms < 0) throw ParseError("Negative atom count in XTC frame");
  header.natoms = static_cast<std::uint32_t>(natoms);

  if (header.natoms > caps_.max_atoms) {
    throw ParseError("Atom count exceeds configured cap");
  }

  if (header.magic == XtcMagic && header.natoms > MaxNatomsOldMagic) {
    throw ParseError("XTC frame exceeds size limit for legacy magic number");
  }

  int step = 0;
  if (!xdr_->read(&step)) {
    return false;
  }
  header.step = static_cast<std::int64_t>(step);

  if (!xdr_->read(&header.time)) {
    return false;
  }
  if (!std::isfinite(header.time)) {
    throw ParseError("Non-finite XTC time in header");
  }

  if (!xdr_->read(header.box.data(), header.box.size())) {
    return false;
  }

  // Validate box finite and reasonable diagonal lengths
  for (float v : header.box) {
    if (!std::isfinite(v)) throw ParseError("Non-finite XTC box element in header");
  }
  // Diagonals in nm at [0], [4], [8]; convert to A for cap check
  const float diag_x = header.box[0] * NmToAngstrom;
  const float diag_y = header.box[4] * NmToAngstrom;
  const float diag_z = header.box[8] * NmToAngstrom;
  auto pos_and_lt_cap = [&](float ang) {
    return (ang > 0.0f) && (ang < caps_.max_box_length_angstrom);
  };
  if (!(pos_and_lt_cap(diag_x) && pos_and_lt_cap(diag_y) && pos_and_lt_cap(diag_z))) {
    throw ParseError("Unreasonable XTC box diagonal in header");
  }

  return true;
}

bool XtcReader::read_uncompressed_coords_into(RDGeom::POINT3D_VECT &dest, std::uint32_t atom_count) {
  if (atom_count > caps_.max_atoms) throw ParseError("Atom count exceeds configured cap");
  std::size_t size3 = 0;
  if (!checked_mul_size_t(static_cast<std::size_t>(atom_count), Dim, size3)) {
    throw ParseError("Coordinate array size overflow");
  }
  std::vector<float> raw(size3);
  if (size3 > 0 && !xdr_->read(raw.data(), raw.size())) {
    return false;
  }

  dest.resize(atom_count);
  for (std::size_t idx = 0, atom = 0; atom < atom_count; ++atom, idx += Dim) {
    dest[atom] = RDGeom::Point3D(
        static_cast<double>(raw[idx])     * NmToAngstrom,
        static_cast<double>(raw[idx + 1]) * NmToAngstrom,
        static_cast<double>(raw[idx + 2]) * NmToAngstrom);
  }
  current_precision_ = -1.0f;
  return true;
}

bool XtcReader::read_compressed_coords_into(RDGeom::POINT3D_VECT &dest, std::uint32_t atom_count) {
  float precision = 0.0f;
  if (!xdr_->read(&precision)) {
    return false;
  }
  if (precision < 1e-12f) throw ParseError("XTC precision too small");
  if (!(precision > 0.0f) || !std::isfinite(precision)) throw ParseError("Invalid XTC precision value");

  current_precision_ = precision;
  const double scale = static_cast<double>(1.0f / precision) * NmToAngstrom;

  std::array<int, Dim> minint{};
  std::array<int, Dim> maxint{};
  if (!xdr_->read(minint.data(), minint.size())) return false;
  if (!xdr_->read(maxint.data(), maxint.size())) return false;

  std::array<unsigned int, Dim> sizeint{};
  for (std::size_t i = 0; i < Dim; ++i) {
    if (maxint[i] < minint[i]) {
      throw ParseError("XTC frame has invalid bounds");
    }
    auto span64 = static_cast<long long>(maxint[i]) - static_cast<long long>(minint[i]) + 1LL;
    if (span64 <= 0) throw ParseError("XTC frame has invalid bounds");
    if (span64 > std::numeric_limits<unsigned int>::max()) throw ParseError("XTC span too large");
    sizeint[i] = static_cast<unsigned int>(span64);
  }

  std::array<int, Dim> bitsizeint{};
  unsigned int bitsize = 0;
  if ((sizeint[0] | sizeint[1] | sizeint[2]) > 0xFFFFFFu) {
    bitsizeint[0] = detail::calculate_sizeof_int(static_cast<int>(sizeint[0]));
    bitsizeint[1] = detail::calculate_sizeof_int(static_cast<int>(sizeint[1]));
    bitsizeint[2] = detail::calculate_sizeof_int(static_cast<int>(sizeint[2]));
    bitsize = 0;
  } else {
    bitsize = static_cast<unsigned int>(detail::calculate_sizeof_ints(span<const unsigned int>(sizeint)));
  }

  // Precompute DivSpecs for per-dimension sizes once per frame (used when bitsize > 0)
  detail::DivSpec size_divspec[Dim] = {
      detail::make_divspec(sizeint[0]),
      detail::make_divspec(sizeint[1]),
      detail::make_divspec(sizeint[2])};

  int smallidx = 0;
  if (!xdr_->read(&smallidx)) return false;
  if (smallidx < detail::FirstSmallIndex || smallidx >= static_cast<int>(detail::MagicInts.size())) {
    throw ParseError("Invalid small index in XTC frame");
  }

  int smaller  = detail::MagicInts[std::max(detail::FirstSmallIndex, smallidx - 1)] / 2;
  int smallnum = detail::MagicInts[smallidx] / 2;
  std::array<unsigned int, Dim> sizesmall{
      static_cast<unsigned int>(detail::MagicInts[smallidx]),
      static_cast<unsigned int>(detail::MagicInts[smallidx]),
      static_cast<unsigned int>(detail::MagicInts[smallidx])};

  std::uint64_t byte_count64 = 0;
  if (current_header_.magic == XtcNewMagic) {
    if (!xdr_->read(&byte_count64)) {
      return false;
    }
  } else {
    unsigned int byte_count32 = 0;
    if (!xdr_->read(&byte_count32)) {
      return false;
    }
    byte_count64 = byte_count32;
  }

  if (byte_count64 > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
    throw ParseError("XTC frame bitstream is too large to load");
  }
  std::size_t byte_count = static_cast<std::size_t>(byte_count64);

  // Config cap: natoms * factor + 64 (avoid overflow)
  std::uint64_t cap64 = static_cast<std::uint64_t>(atom_count) * static_cast<std::uint64_t>(caps_.max_bitstream_bytes_per_atom) + 64ull;
  if (byte_count64 > cap64) throw ParseError("Unreasonable XTC payload size for atom count");

  // Pad with 8 zero bytes to enable safe unaligned 64-bit loads in the bit reader
  bitbuf_.resize(byte_count + 8);
  if (byte_count > 0 && !xdr_->read(reinterpret_cast<char *>(bitbuf_.data()), byte_count)) {
    return false;
  }

  detail::BitReader reader(bitbuf_.data(), byte_count, 8);
  if (atom_count > caps_.max_atoms) throw ParseError("Atom count exceeds configured cap");
  std::size_t size3 = 0;
  if (!checked_mul_size_t(static_cast<std::size_t>(atom_count), Dim, size3)) {
    throw ParseError("Coordinate array size overflow");
  }

  coordbuf_.resize(size3);
  dest.resize(atom_count);

  std::array<int, 3> prevcoord{};
  int run = 0; // Persist across atoms which is what GROMACS also does
  std::size_t out = 0;

  // Precompute fast-div contexts for MagicInts table once per frame
  std::array<detail::DivSpec, detail::MagicInts.size()> magic_divspec{};
  detail::build_magicints_divspec(magic_divspec);

  while (out < atom_count) {
    int *thiscoord = coordbuf_.data() + out * Dim;

    if (bitsize == 0) {
      thiscoord[0] = reader.read_bits(bitsizeint[0]);
      thiscoord[1] = reader.read_bits(bitsizeint[1]);
      thiscoord[2] = reader.read_bits(bitsizeint[2]);
    } else {
      // Fast mixed-radix specialization with precomputed divisors
      detail::decode_mixed3_precomp(
          reader,
          static_cast<unsigned>(bitsize),
          size_divspec,
          thiscoord);
    }

    thiscoord[0] += minint[0];
    thiscoord[1] += minint[1];
    thiscoord[2] += minint[2];

    prevcoord[0] = thiscoord[0];
    prevcoord[1] = thiscoord[1];
    prevcoord[2] = thiscoord[2];

    int flag = reader.read_bits(1);
    int is_smaller = 0;
    if (flag == 1) {
      run        = reader.read_bits(5);
      is_smaller = run % 3;
      run       -= is_smaller;
      --is_smaller;
    }

    if (run > 0) {
      // Write the first coordinate (prevcoord) and then the subsequent ones
      dest[out++] = make_point(prevcoord.data(), scale);

      for (int k = 0; k < run; k += 3) {
        thiscoord += 3;
        // Fast specialized path: 3 ints, identical base MagicInts[smallidx]
        detail::decode_small3_precomp(
            reader,
            static_cast<unsigned>(smallidx),
            magic_divspec[static_cast<std::size_t>(smallidx)],
            thiscoord);

        thiscoord[0] += prevcoord[0] - smallnum;
        thiscoord[1] += prevcoord[1] - smallnum;
        thiscoord[2] += prevcoord[2] - smallnum;

        prevcoord[0] = thiscoord[0];
        prevcoord[1] = thiscoord[1];
        prevcoord[2] = thiscoord[2];

        if (out >= atom_count) throw ParseError("XTC frame produced too many coordinates");
        dest[out++] = make_point(thiscoord, scale);
      }
    } else {
      dest[out++] = make_point(thiscoord, scale);
    }

    smallidx += is_smaller;
    smallidx  = std::clamp(smallidx, detail::FirstSmallIndex - 1, static_cast<int>(detail::MagicInts.size()) - 1);

    if (is_smaller < 0) {
      smallnum = smaller;
      smaller  = smallidx > detail::FirstSmallIndex ? detail::MagicInts[smallidx - 1] / 2 : 0;
    } else if (is_smaller > 0) {
      smaller  = smallnum;
      smallnum = detail::MagicInts[smallidx] / 2;
    }
  }

  if (out != atom_count) throw ParseError("XTC frame decoded incorrect atom count");

  // Allow residual zero padding bits in last byte?
  if (strict_ && !reader.residual_bits_zero()) {
    throw ParseError("Non-zero padding bits at end of XTC bitstream");
  }
  reader.align_to_byte();
  return true;
}

bool XtcReader::is_seekable() const {
  if (!stream_) return false;
  auto *buf = stream_->rdbuf();
  if (!buf) return false;
  // Try a no-op seek via the buffer, no state change
  auto cur = stream_->tellg();
  if (cur < 0) return false;
  auto test = buf->pubseekoff(0, std::ios::cur, std::ios::in);
  if (test == std::streampos(-1)) return false;
  return true;
}

std::streampos XtcReader::file_size_bytes() const {
  auto *ifs = dynamic_cast<std::ifstream *>(stream_);
  if (!ifs) return std::streampos(-1);
  auto cur = ifs->tellg();
  ifs->seekg(0, std::ios::end);
  auto end = ifs->tellg();
  ifs->clear();
  ifs->seekg(cur);
  return end;
}

std::optional<std::pair<std::int64_t, float>> XtcReader::at_header_start() {
  // Probe for header pattern (magic, natoms, step, time, box) and validate natoms and box pattern
  auto *ifs = dynamic_cast<std::ifstream *>(stream_);
  if (!ifs) return std::nullopt;
  auto off = ifs->tellg();
  int magic = 0, natoms_i = 0, step_i = 0;
  std::array<float, 10> f_all{};
  internal::XDRReader probe(stream_);
  // Attempt reads, on any failure, move to off + 4 and return nullopt.
  auto advance4 = [&]() {
    ifs->clear();
    ifs->seekg(off + static_cast<std::streamoff>(sizeof(internal::XDRReader::block_type)));
  };
  if (!probe.read(&magic))    { advance4(); return std::nullopt; }
  if (!probe.read(&natoms_i)) { advance4(); return std::nullopt; }
  if (!probe.read(&step_i))   { advance4(); return std::nullopt; }
  // Read time + 9 box floats
  for (int i = 0; i < 10; ++i) {
    if (!probe.read(&f_all[static_cast<std::size_t>(i)])) {
      advance4();
      return std::nullopt;
    }
  }
  const float time = f_all[0];
  // Validate
  const bool magic_ok  = (magic == XtcMagic || magic == XtcNewMagic);
  const bool natoms_ok = (natoms_i >= 0);

  constexpr float eps = 1e-12f;
  auto finite = [](float f){ return std::isfinite(f); };
  bool box_finite = std::all_of(f_all.begin()+1, f_all.end(), finite);

  // Diagonals are box[0], box[4], box[8] --> f_all[1], f_all[5], f_all[9]
  bool diag_reasonable = (std::fabs(f_all[1]) > eps) &&
                         (std::fabs(f_all[5]) > eps) &&
                         (std::fabs(f_all[9]) > eps);

  // Peek lsize and check equals natoms
  unsigned int lsize = 0;
  if (!probe.read(&lsize)) { advance4(); return std::nullopt; }
  bool lsize_ok = (static_cast<int>(lsize) == natoms_i);

  const bool box_ok = box_finite && diag_reasonable;
  // Reset to off + 4 for linear scan continuation
  ifs->clear();
  ifs->seekg(off + static_cast<std::streamoff>(sizeof(internal::XDRReader::block_type)));
  if (!magic_ok || !natoms_ok || !box_ok || !lsize_ok) return std::nullopt;
  const std::int64_t step = static_cast<std::int64_t>(step_i);
  return std::make_optional(std::make_pair(step, time));
}

std::optional<std::streampos> XtcReader::find_next_frame_start() {
  auto *ifs = dynamic_cast<std::ifstream *>(stream_);
  if (!ifs) return std::nullopt;
  auto fsize = file_size_bytes();
  if (fsize <= 0) return std::nullopt;
  while (true) {
    auto pos = ifs->tellg();
    if (pos < 0) return std::nullopt;

    // Stop scanning if we are very close to the end // no full header can fit
    constexpr std::streamoff MinHeaderBytes = 4 + 4 + 4 + 4 + 9 * 4; // 52
    if (pos >= fsize - MinHeaderBytes) return std::nullopt;

    if (auto hdr = at_header_start(); hdr.has_value()) {
      // Header was at pos (we advanced by 4, reset to pos)
      ifs->clear();
      ifs->seekg(pos);
      return pos;
    }
    // Keep scanning, at_header_start already advanced by 4 bytes
    if (ifs->fail()) {
      ifs->clear();
      // Move forward by 4 bytes to try to recover
      ifs->seekg(pos + static_cast<std::streamoff>(sizeof(internal::XDRReader::block_type)));
    }
  }
}

namespace {

void skip_frame_payload(internal::XDRReader &reader, int magic, std::uint32_t lsize, unsigned int min_compressed_system_size) {
  if (lsize <= min_compressed_system_size) {
    std::vector<float> tmp(static_cast<std::size_t>(lsize) * 3);
    if (!reader.read(tmp.data(), tmp.size())) {
      throw ParseError("XTC truncated while skipping coordinates");
    }
    return;
  }
  float precision = 0.0f;
  if (!reader.read(&precision)) {
    throw ParseError("XTC missing precision during skip");
  }
  std::array<int, 3> minint{};
  std::array<int, 3> maxint{};
  int smallidx;
  if (!reader.read(minint.data(), minint.size()) ||
      !reader.read(maxint.data(), maxint.size()) ||
      !reader.read(&smallidx)) {
    throw ParseError("XTC truncated while skipping bounds");
  }
  std::uint64_t byte_count64 = 0;
  if (magic == XtcNewMagic) {
    if (!reader.read(&byte_count64)) {
      throw ParseError("XTC truncated while skipping payload size");
    }
  } else {
    unsigned int byte_count32 = 0;
    if (!reader.read(&byte_count32)) {
      throw ParseError("XTC truncated while skipping payload size");
    }
    byte_count64 = byte_count32;
  }
  if (byte_count64 > 0) {
    if (byte_count64 > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
      throw ParseError("XTC payload too large to skip");
    }
    std::vector<char> skip(static_cast<std::size_t>(byte_count64));
    if (!reader.read(skip.data(), skip.size())) {
      throw ParseError("XTC truncated while skipping payload");
    }
  }
}

} // namespace

bool XtcReader::ensure_index() {
  if (frame_index_ready_) return true;
  if (filename_.empty())  return false;

  std::ifstream idx_stream(filename_, std::ios::binary);
  if (!idx_stream.good()) return false;

  internal::XDRReader idx_reader(&idx_stream);

  frame_index_.clear();

  try {
    while (true) {
      auto header_pos = idx_stream.tellg();
      if (header_pos < 0) {
        break;
      }
      int magic = 0;
      if (!idx_reader.read(&magic)) break; // EOF
      if (magic != XtcMagic && magic != XtcNewMagic) {
        break; // Stop indexing at inconsistency and keep partial index
      }
      int natoms = 0;
      int step_i = 0;
      float time = 0.0f;
      if (!idx_reader.read(&natoms) || !idx_reader.read(&step_i) || !idx_reader.read(&time)) {
        break;
      }
      std::array<float, 9> box_tmp{};
      if (!idx_reader.read(box_tmp.data(), box_tmp.size())) {
        break;
      }
      frame_index_.push_back(FrameIndexEntry{header_pos, static_cast<std::int64_t>(step_i), time});

      std::uint32_t lsize = 0;
      if (!idx_reader.read(&lsize)) {
        break;
      }
      if (lsize != static_cast<std::uint32_t>(natoms)) {
        break; // stop indexing on mismatch to avoid desync
      }
      skip_frame_payload(idx_reader, magic, lsize, min_compressed_system_size);
    }
  } catch (const ParseError &) {} // tolerate partial index, swallow

  if (frame_index_.empty()) return false;

  frame_index_ready_ = true;
  if (has_current_) {
    auto it = std::find_if(frame_index_.begin(), frame_index_.end(), [&](const FrameIndexEntry &entry) {
      return entry.step == current_header_.step && std::fabs(entry.time - current_header_.time) < 1e-6f;
    });
    if (it != frame_index_.end()) {
      current_index_ = static_cast<std::size_t>(std::distance(frame_index_.begin(), it));
    } else {
      current_index_ = 0;
    }
  } else {
    current_index_ = 0;
  }

  return true;
}

bool XtcReader::seek_to_index(std::size_t idx) {
  if (!frame_index_ready_ || idx >= frame_index_.size()) {
    return false;
  }
  auto *ifs = dynamic_cast<std::ifstream *>(stream_);
  if (!ifs) return false;

  ifs->clear();
  ifs->seekg(frame_index_[idx].offset);
  if (!ifs->good()) return false;

  xdr_ = std::make_unique<internal::XDRReader>(stream_);
  if (!read_next_frame()) return false;

  current_index_ = idx;
  return true;
}

bool XtcReader::seek_step(std::int64_t target_step) {
  if (!is_seekable())  return false;
  if (!ensure_index()) return false;

  if (frame_index_ready_ && current_index_ < frame_index_.size()
      && frame_index_[current_index_].step == target_step) {
    return true;
  }

  auto it = std::lower_bound(
      frame_index_.begin(),
      frame_index_.end(),
      target_step,
      [](const FrameIndexEntry &entry, std::int64_t value) { return entry.step < value; });

  if (it == frame_index_.end() || it->step != target_step) return false;
  return seek_to_index(static_cast<std::size_t>(std::distance(frame_index_.begin(), it)));
}

bool XtcReader::seek_time(float target_time, bool forward_only) {
  if (!is_seekable())  return false;
  if (!ensure_index()) return false;

  if (frame_index_ready_ && current_index_ < frame_index_.size()
      && std::fabs(frame_index_[current_index_].time - target_time) < 1e-6f) {
    return true;
  }

  auto it = std::lower_bound(
      frame_index_.begin(),
      frame_index_.end(),
      target_time,
      [](const FrameIndexEntry &entry, float value) { return entry.time < value; });

  std::size_t idx;
  if (it == frame_index_.end()) {
    idx = frame_index_.size() - 1;
  } else {
    idx = static_cast<std::size_t>(std::distance(frame_index_.begin(), it));
  }

  if (forward_only && frame_index_ready_) {
    if (idx < current_index_) return false;
  }

  return seek_to_index(idx);
}

} // namespace lahuta::md
