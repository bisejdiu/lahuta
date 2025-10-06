#pragma once

#include <string>
#include <string_view>

#include <rdkit/GraphMol/MonomerInfo.h>

// TODO: Here we duplicate code that exist elsewhere in the project.
namespace lahuta::bonds {

struct ResidueKey {
  std::string resname;
  std::string chain_id;
  std::string icode;
  std::string alt_loc;
  int resnum = 0;

  ResidueKey() = default;

  explicit ResidueKey(const RDKit::AtomPDBResidueInfo &info)
      : resname(info.getResidueName()), chain_id(info.getChainId()), icode(info.getInsertionCode()),
        alt_loc(info.getAltLoc()), resnum(info.getResidueNumber()) {}

  bool operator==(const ResidueKey &other) const {
    return resname == other.resname && chain_id == other.chain_id && icode == other.icode
           && alt_loc == other.alt_loc && resnum == other.resnum;
  }

  bool operator!=(const ResidueKey &other) const { return !(*this == other); }
};

namespace residue_key_utils {

inline std::string trim_copy(std::string_view value) {
  const auto first = value.find_first_not_of(' ');
  if (first == std::string_view::npos) return {};
  const auto last = value.find_last_not_of(' ');
  return std::string(value.substr(first, last - first + 1));
}

inline std::string build_instance_key(const ResidueKey &key) {
  std::string result;
  result.reserve(key.resname.size() + key.chain_id.size() + 24);
  result.append(key.resname);
  result.push_back('|');
  result.append(key.chain_id);
  result.push_back('|');
  result.append(std::to_string(key.resnum));
  result.push_back('|');
  result.append(key.icode);
  result.push_back('|');
  result.append(key.alt_loc);
  return result;
}

} // namespace residue_key_utils

} // namespace lahuta::bonds

namespace std {
template <> struct hash<lahuta::bonds::ResidueKey> {
  size_t operator()(const lahuta::bonds::ResidueKey &key) const {
    size_t h1 = hash<string>{}(key.resname);
    size_t h2 = hash<string>{}(key.chain_id);
    size_t h3 = hash<string>{}(key.icode);
    size_t h4 = hash<string>{}(key.alt_loc);
    size_t h5 = hash<int>{}(key.resnum);

    return h1 ^ (h2 << 1) ^ (h3 << 2) ^ (h4 << 3) ^ (h5 << 4);
  }
};
} // namespace std
