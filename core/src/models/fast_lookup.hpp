/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto t1 = std::make_tuple("besian");
 *   auto t2 = std::make_tuple("sejdiu");
 *   auto t3 = std::make_tuple("@gmail.com");
 *   auto combined = std::tuple_cat(t1, t2, t3);
 *   return std::apply([](auto... args) {
 *     return (std::string{} + ... + std::string(args));
 *   }, combined);
 * }();
 *
 */

#ifndef LAHUTA_MODEL_LOOKUP_HPP
#define LAHUTA_MODEL_LOOKUP_HPP

#include <array>

namespace lahuta {

static const std::array<int, 256> StandardAminoAcidAtomicNumbers = []() {
    std::array<int, 256> result{};
    result['C'] = 6;
    result['O'] = 8;
    result['N'] = 7;
    result['S'] = 16;
    return result;
}();

// FIX: not needed as we can use StandardAminoAcidDataTable[x].size
static const std::array<int, 256> StandardAminoAcidAtomSizeTable = []() {
    std::array<int, 256> result{};
    result['G'] = 4;
    result['A'] = 5;
    result['V'] = 7;
    result['L'] = 8;
    result['I'] = 8;
    result['S'] = 6;
    result['T'] = 7;
    result['C'] = 6;
    result['M'] = 8;
    result['P'] = 7;
    result['F'] = 11;
    result['Y'] = 12;
    result['W'] = 14;
    result['H'] = 10;
    result['E'] = 9;
    result['D'] = 8;
    result['N'] = 8;
    result['Q'] = 9;
    result['K'] = 9;
    result['R'] = 11;
    return result;
}();

} // namespace lahuta

#endif // LAHUTA_MODEL_LOOKUP_HPP
