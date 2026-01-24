/*
 * Aims to provide a compile-time substitution matrix. Most of the challenge comes from the non-ideal input
 * file format. Performance-wise, there isn't much to gain from compile-time evaluation, but it's a fun
 * exercise.
 */
#ifndef LAHUTA_SUBSTITUTION_MATRIX_H
#define LAHUTA_SUBSTITUTION_MATRIX_H

#include <array>
#include <climits>

#include "BaseMatrix.h"
#include "fseek/internal/constexpr_utils.hpp"
#include "fseek/ops.hpp"

namespace lahuta {

namespace {
constexpr size_t skipNonWhitespace(const unsigned char *data, size_t size, size_t pos) {
  size_t counter = pos;
  while (counter < size && !isspace(data[counter])) {
    counter++;
  }
  return counter - pos;
}

constexpr size_t skipWhitespace(const unsigned char *data, size_t size, size_t pos) {
  size_t counter = pos;
  while (counter < size && isspace(data[counter]) && data[counter] != '\n') {
    counter++;
  }
  return counter - pos;
}

constexpr size_t getWordsOfLine(
    const unsigned char *data, size_t size, size_t startPos, const unsigned char **words, size_t maxElement) {
  size_t pos = startPos;
  size_t elementCounter = 0;

  while (pos < size && data[pos] != '\n') {
    pos += skipWhitespace(data, size, pos);
    if (pos >= size || data[pos] == '\n') break;

    words[elementCounter] = &data[pos];
    elementCounter++;
    if (elementCounter >= maxElement) return elementCounter;

    pos += skipNonWhitespace(data, size, pos);
  }

  if (elementCounter < maxElement && pos < size) {
    words[elementCounter] = &data[pos];
  }
  return elementCounter;
}

/// NOTE: unsafe! Only for use in the context of this file
constexpr bool compare_str(const unsigned char *str1, const char *str2) {
  size_t i = 0;
  while (str2[i] != '\0') {
    if (str1[i] != str2[i]) return false;
    i++;
  }
  return true;
}

constexpr double strtod(const unsigned char *str) {
  double result = 0.0;
  bool negative = false;

  if (*str == '-') {
    negative = true;
    ++str;
  }

  // integer part
  while (*str >= '0' && *str <= '9') {
    result = result * 10.0 + (*str - '0');
    ++str;
  }

  // decimal part
  if (*str == '.') {
    ++str;
    double fraction = 0.1;
    while (*str >= '0' && *str <= '9') {
      result += (*str - '0') * fraction;
      fraction *= 0.1;
      ++str;
    }
  }

  return negative ? -result : result;
}

constexpr std::array<unsigned char, UCHAR_MAX> make_filled_char_array() {
  std::array<unsigned char, UCHAR_MAX> arr{};
  for (auto &elem : arr) {
    elem = static_cast<unsigned char>(UCHAR_MAX);
  }
  return arr;
}

constexpr std::array<char, UCHAR_MAX> make_filled_uchar_array() {
  std::array<char, UCHAR_MAX> arr{};
  for (auto &elem : arr) {
    elem = static_cast<char>(UCHAR_MAX);
  }
  return arr;
}
} // namespace

struct __SubMatrix__ {
  static constexpr int alphabetSize = 21;
  static constexpr int containsX = 1;
  constexpr static const double ANY_BACK = 1E-5;

  constexpr __SubMatrix__(const CompileTimeMatrix *matrix, float bitFactor, float scoreBias)
      : aa2num(make_filled_char_array()), num2aa(make_filled_uchar_array()), bitFactor(bitFactor),
        scoreBias(scoreBias) {

    setAaMappingDetectAlphSize(matrix->data, matrix->size);
    readProbMatrix(matrix->data, matrix->size);
    setupLetterMapping();
    generateSubMatrix();
  }

  using Matrix = std::array<std::array<short, alphabetSize>, alphabetSize>;
  using FloatMatrix = std::array<std::array<float, alphabetSize>, alphabetSize>;
  using DoubleMatrix = std::array<std::array<double, alphabetSize>, alphabetSize>;
  using ProbArray = std::array<double, alphabetSize>;

  constexpr void generateSubMatrix() noexcept {
    ProbArray localPBack{};
    computeBackground(localPBack);

    // Precompute matrix R for amino acid pseudocounts
    for (int i = 0; i < alphabetSize; i++) {
      for (int j = 0; j < alphabetSize; j++) {
        subMatrixPseudoCounts[i][j] = probMatrix[i][j] / pBack[j];
      }
    }

    DoubleMatrix tempSubMatrix{};
    for (int i = 0; i < alphabetSize; i++) {
      for (int j = 0; j < alphabetSize; j++) {
        double v = probMatrix[i][j] / (pBack[i] * pBack[j]);
        tempSubMatrix[i][j] = log2_approx(v);
      }
    }

    // Convert to short data type matrix with bit scaling
    for (int i = 0; i < alphabetSize; i++) {
      for (int j = 0; j < alphabetSize; j++) {
        double pValNBitScale = (bitFactor * tempSubMatrix[i][j] + scoreBias);
        subMatrix[i][j] = (pValNBitScale < 0.0) ? pValNBitScale - 0.5 : pValNBitScale + 0.5;
      }
    }
  }

  constexpr void computeBackground(ProbArray &local_pback) noexcept {
    for (int i = 0; i < alphabetSize; i++) {
      local_pback[i] = 0;
      for (int j = 0; j < alphabetSize; j++) {
        local_pback[i] += probMatrix[i][j];
      }
    }
    if (containsX) {
      local_pback[alphabetSize - 1] = ANY_BACK;
    }
  }

  constexpr std::pair<int, bool> setAaMappingDetectAlphSize(const unsigned char *matrixData, size_t size) {
    const unsigned char *words[256] = {};
    size_t pos = 0;
    size_t lineStart = 0;

    while (pos < size) {
      if (matrixData[pos] == '#') {
        while (pos < size && matrixData[pos] != '\n')
          pos++;
        pos++;
        lineStart = pos;
        continue;
      }

      size_t wordCnt = getWordsOfLine(matrixData, size, lineStart, words, 256);

      if (wordCnt > 1) {
        for (size_t i = 0; i < wordCnt; i++) {
          if (!isalpha(words[i][0])) {
            return {-1, false};
          }
          int aa = __toupper__(words[i][0]);
          aa2num[aa] = static_cast<unsigned char>(i);
          num2aa[i] = aa;
          auto v = static_cast<unsigned char>(i);
        }
        return {static_cast<int>(wordCnt), containsX};
      }

      while (pos < size && matrixData[pos] != '\n')
        pos++;
      pos++;
      lineStart = pos;
    }

    return {-1, false};
  }

  constexpr void readProbMatrix(const unsigned char *matrixData, size_t size) {
    const unsigned char *words[256] = {};
    size_t pos = 0;
    size_t lineStart = 0;

    bool probMatrixStart = false;
    bool hasLambda = false;
    bool hasBackground = false;

    while (pos < size) {
      if (matrixData[pos] == '#') {
        if (compare_str(&matrixData[pos], "# Background (precomputed optional):")) {
          size_t wordCnt = getWordsOfLine(matrixData, size, pos, words, 256);
          for (size_t i = 4; i < wordCnt; i++) {
            auto v = words[i];
            pBack[i - 4] = strtod(words[i]);
          }
          hasBackground = true;
        } else if (compare_str(&matrixData[pos], "# Lambda     (precomputed optional):")) {
          size_t wordCnt = getWordsOfLine(matrixData, size, pos, words, 256);
          lambda = strtod(words[4]);
          hasLambda = true;
        }

        while (pos < size && matrixData[pos] != '\n')
          pos++;
        pos++;
        lineStart = pos;
        continue;
      }

      size_t wordCnt = getWordsOfLine(matrixData, size, lineStart, words, 256);

      if (wordCnt > 1) {
        if (!probMatrixStart) {
          probMatrixStart = true;
        } else {
          unsigned char upperChar = __toupper__(words[0][0]);
          size_t aa = aa2num[upperChar];
          for (size_t i = 0; i < alphabetSize; i++) {
            probMatrix[aa][i] = strtod(words[i + 1]);
          }
        }
      }

      while (pos < size && matrixData[pos] != '\n')
        pos++;
      pos++;
      lineStart = pos;
    }

    bool xIsPositive = false;
    if (containsX) {
      int xIndex = aa2num['X'];
      for (int j = 0; j < alphabetSize; j++) {
        if ((probMatrix[xIndex][j] > 0) || (probMatrix[j][xIndex] > 0)) {
          xIsPositive = true;
          break;
        }
      }
    }

    if (!containsX) return; // here for visibility

    if (!hasLambda || !hasBackground) {
      // NOTE: mmseqs invokes estimateLambdaAndBackground
      pBack[aa2num['X']] = ANY_BACK;
    }

    if (!xIsPositive) {
      for (int i = 0; i < alphabetSize - 1; i++) {
        pBack[i] = pBack[i] * (1.0 - pBack[aa2num['X']]);
      }
    }

    // Reconstruct Probability Sab=(Pab/Pa*Pb) -> Pab = exp(Sab) * Pa * Pb
    for (int i = 0; i < alphabetSize; i++) {
      for (int j = 0; j < alphabetSize; j++) {
        probMatrix[i][j] = exp_approx(lambda * probMatrix[i][j]) * pBack[i] * pBack[j];
      }
    }
  }

  constexpr void setupLetterMapping() {
    for (int letter = 0; letter < UCHAR_MAX; letter++) {
      char upperLetter = __toupper__(static_cast<char>(letter));
      switch (upperLetter) {
        case 'A':
        case 'T':
        case 'G':
        case 'C':
        case 'D':
        case 'E':
        case 'F':
        case 'H':
        case 'I':
        case 'K':
        case 'L':
        case 'M':
        case 'N':
        case 'P':
        case 'Q':
        case 'R':
        case 'S':
        case 'V':
        case 'W':
        case 'Y':
        case 'X':
          this->aa2num[static_cast<int>(letter)] = this->aa2num[static_cast<int>(upperLetter)];
          break;
        case 'J':
          this->aa2num[static_cast<int>(letter)] = this->aa2num[(int)'L'];
          break;
        case 'U':
        case 'O':
          this->aa2num[static_cast<int>(letter)] = this->aa2num[(int)'X'];
          break;
        case 'Z':
          this->aa2num[static_cast<int>(letter)] = this->aa2num[(int)'E'];
          break;
        case 'B':
          this->aa2num[static_cast<int>(letter)] = this->aa2num[(int)'D'];
          break;
        default:
          this->aa2num[static_cast<int>(letter)] = this->aa2num[(int)'X'];
          break;
      }
    }
  }

  std::array<unsigned char, UCHAR_MAX> aa2num{};
  std::array<char, UCHAR_MAX> num2aa{};

  Matrix subMatrix{};
  FloatMatrix subMatrixPseudoCounts{};
  DoubleMatrix probMatrix{};
  ProbArray pBack{};

  double lambda{};

  const float bitFactor;
  const float scoreBias;
};

inline BaseMatrix *create_base_matrix(const __SubMatrix__ &source) {
  BaseMatrix *ptr = new BaseMatrix;
  ptr->alphabetSize = source.alphabetSize;
  ptr->allocatedAlphabetSize = source.alphabetSize;
  ptr->lambda = source.lambda;
  ptr->scoreBias = source.scoreBias;

  ptr->initMatrixMemory(ptr->alphabetSize);

  ptr->num2aa = new char[255];
  ptr->aa2num = new unsigned char[UCHAR_MAX];
  for (int i = 0; i < UCHAR_MAX; ++i) {
    ptr->aa2num[i] = UCHAR_MAX;
  }

  std::memcpy(ptr->aa2num, source.aa2num.data(), UCHAR_MAX);
  std::memcpy(ptr->num2aa, source.num2aa.data(), UCHAR_MAX);

  for (int i = 0; i < ptr->alphabetSize; ++i) {
    ptr->subMatrix[i] = new short[ptr->alphabetSize];
    std::memcpy(ptr->subMatrix[i], source.subMatrix[i].data(), ptr->alphabetSize * sizeof(short));
  }

  for (int i = 0; i < ptr->alphabetSize; ++i) {
    ptr->subMatrixPseudoCounts[i] = new float[ptr->alphabetSize];
    std::memcpy(
        ptr->subMatrixPseudoCounts[i],
        source.subMatrixPseudoCounts[i].data(),
        ptr->alphabetSize * sizeof(float));
  }

  for (int i = 0; i < ptr->alphabetSize; ++i) {
    ptr->probMatrix[i] = new double[ptr->alphabetSize];
    std::memcpy(ptr->probMatrix[i], source.probMatrix[i].data(), ptr->alphabetSize * sizeof(double));
  }

  std::memcpy(ptr->pBack, source.pBack.data(), ptr->alphabetSize * sizeof(double));

  return ptr;
}

inline constexpr __SubMatrix__ CTSM8 = __SubMatrix__(&SubMatrix3Di, 8.0f, -0.2f);
inline constexpr __SubMatrix__ CTSM2 = __SubMatrix__(&SubMatrix3Di, 2.0f, -0.2f);

} // namespace lahuta

#endif // LAHUTA_SUBSTITUTION_MATRIX_H
