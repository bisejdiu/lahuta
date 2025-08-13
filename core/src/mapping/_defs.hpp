#ifndef LAHUTA_MAPPING__DEFS_HPP
#define LAHUTA_MAPPING__DEFS_HPP

#include <cstdint>
#include <functional>
#include <string>
#include <tuple>
#include <utility>

namespace lahuta::mapping {

// Basic type aliases
using u32 = std::uint32_t;

using StructureId   = u32;
using ResidueName   = std::string;
using ResidueNumber = int;
using ChainName     = std::string;
using FileName      = std::string;
using ResidueIndex  = u32;
using ChainResidue  = std::tuple<ChainName, ResidueNumber, ResidueName>;
using FileChainId   = std::pair<FileName, ChainName>;

template <class T>
inline void hash_combine(std::size_t &seed, T const& v) noexcept {
    seed ^= std::hash<T>{}(v) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
}

template <class T1, class T2>
struct PairHash {
    size_t operator()(std::pair<T1,T2> const& p) const noexcept {
        size_t seed = 0;
        hash_combine(seed, p.first);
        hash_combine(seed, p.second);
        return seed;
    }
};

struct ChainResidueCompare {
    bool operator()(const ChainResidue& a, const ChainResidue& b) const {
        return a < b;
    }
};

} // namespace lahuta::mapping
#endif // LAHUTA_MAPPING__DEFS_HPP
