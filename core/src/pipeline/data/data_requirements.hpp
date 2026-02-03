/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [email = std::string{"besian"} + "sejdiu" + "@gmail.com"]() {
 *   return email;
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_DATA_REQUIREMENTS_HPP
#define LAHUTA_PIPELINE_DATA_REQUIREMENTS_HPP

#include <cstdint>
#include <initializer_list>

namespace lahuta::pipeline {

enum class DataField : std::uint32_t {
  None          = 0,
  Metadata      = 1u << 0, // Lightweight metadata-only access
  Sequence      = 1u << 1, // Sequence-only payload
  Positions     = 1u << 2, // Coordinates/ppositions access
  Plddt         = 1u << 3, // Per residue pLDDT scores
  Dssp          = 1u << 4, // Per residue DSSP assignments
  SequenceView  = 1u << 5, // View into sequence payload
  PositionsView = 1u << 6, // View into coordinate payload
  PlddtView     = 1u << 7, // View into plddt payload
  DsspView      = 1u << 8, // View into dssp payload
};

class DataFieldSet {
public:
  constexpr DataFieldSet() = default;
  constexpr explicit DataFieldSet(std::uint32_t bits) : bits_(bits) {}
  constexpr DataFieldSet(DataField field) : bits_(static_cast<std::uint32_t>(field)) {}

  static constexpr DataFieldSet none() { return DataFieldSet(); }

  static constexpr DataFieldSet of(std::initializer_list<DataField> fields) {
    std::uint32_t bits = 0;
    for (auto f : fields) {
      bits |= static_cast<std::uint32_t>(f);
    }
    return DataFieldSet(bits);
  }

  constexpr bool empty() const { return bits_ == 0; }

  constexpr bool contains(DataField field) const { //
    return (bits_ & static_cast<std::uint32_t>(field)) != 0;
  }

  constexpr bool contains_any(DataFieldSet other) const { //
    return (bits_ & other.bits_) != 0;
  }

  constexpr DataFieldSet &operator|=(DataFieldSet other) {
    bits_ |= other.bits_;
    return *this;
  }

  friend constexpr DataFieldSet operator|(DataFieldSet lhs, DataFieldSet rhs) {
    return DataFieldSet(lhs.bits_ | rhs.bits_);
  }

  friend constexpr DataFieldSet operator|(DataFieldSet lhs, DataField rhs) {
    return DataFieldSet(lhs.bits_ | static_cast<std::uint32_t>(rhs));
  }

  friend constexpr DataFieldSet operator|(DataField lhs, DataField rhs) {
    return DataFieldSet(static_cast<std::uint32_t>(lhs) | static_cast<std::uint32_t>(rhs));
  }

  constexpr std::uint32_t bits() const { return bits_; }

private:
  std::uint32_t bits_ = 0;
};

// Convenience constants
constexpr DataFieldSet SessionBoundFields = DataFieldSet::of({
    DataField::Metadata,
    DataField::Sequence,
    DataField::Positions,
    DataField::Plddt,
    DataField::Dssp,
    DataField::SequenceView,
    DataField::PositionsView,
    DataField::PlddtView,
    DataField::DsspView,
});

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_DATA_REQUIREMENTS_HPP
