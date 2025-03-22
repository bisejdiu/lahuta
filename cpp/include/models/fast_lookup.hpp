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
