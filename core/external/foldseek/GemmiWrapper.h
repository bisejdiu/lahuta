//
// Created by Martin Steinegger on 6/7/21.
//

#ifndef FOLDSEEK_GEMMIWRAPPER_H
#define FOLDSEEK_GEMMIWRAPPER_H
#include <vector>
#include <unordered_map>
#include <string>
#include "gemmi/model.hpp"
#include "structureto3di.h"

class GemmiWrapper {
public:
    enum class Format {
        Detect = 0,
        Pdb,
        Mmcif,
        Mmjson,
        ChemComp,
        Foldcomp,
        Unknown
    };

    GemmiWrapper();
    ~GemmiWrapper() {
        if (fixupBuffer) {
            delete fixupBuffer;
        }
    }

    bool load(const std::string& filename, Format format = Format::Detect);

    std::pair<size_t, size_t> nextChain();

    std::vector<Alphabet3Di::Vec3> ca; // C-alpha coordinates
    std::vector<float> ca_bfactor; // C-alpha B-factors
    std::vector<Alphabet3Di::Vec3> n; // N coordinates
    std::vector<Alphabet3Di::Vec3> c; // C coordinates
    std::vector<Alphabet3Di::Vec3> cb; // C-beta coordinates
    std::vector<char> ami; // Amino acid sequence
    std::vector<char> seq3di;
    std::vector<std::string> names;
    std::vector<std::string> chainNames;
    std::vector<unsigned int> modelIndices;
    unsigned int modelCount = 0;
    std::vector<std::pair<size_t ,size_t>> chain;
    std::vector<int> taxIds;
    std::string title;

    char* fixupBuffer;
    size_t fixupBufferSize;

    gemmi::Structure st;

private:
    std::unordered_map<std::string,char> threeAA2oneAA;
    int modelIt;
    int chainIt;

    void updateStructure(void * structure, const std::string & filename, std::unordered_map<std::string, int>& entity_to_tax_id);
};


#endif //FOLDSEEK_GEMMIWRAPPER_H
