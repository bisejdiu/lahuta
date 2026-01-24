#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "fseek/ops.hpp"
#include "fseek/prefilter.hpp"

namespace py = pybind11;

// clang-format off
namespace {
// Uses a std::vector<char> as a buffer
static std::string result_to_buffer_wrapper(const Matcher::result_t &result, bool addBacktrace, bool compress, bool addOrfPosition) {
    std::vector<char> buff(1024);
    size_t len = Matcher::resultToBuffer(buff.data(), result, addBacktrace, compress, addOrfPosition);
    return std::string(buff.data(), len);
}
} // namespace


namespace lahuta::bindings {

void bind_foldseek(py::module_ &m) {

    py::enum_<AlignType>(m, "AlignType")
        .value("AA_3Di",  AlignType::AA_3Di)
        .value("AA",      AlignType::AA)
        .value("_3Di",    AlignType::_3Di);

    py::enum_<TMScoreThrMode> (m, "TMScoreThrMode")
        .value("alignment", TMScoreThrMode::alignment)
        .value("query",     TMScoreThrMode::query)
        .value("target",    TMScoreThrMode::target)
        .value("min",       TMScoreThrMode::min);

    py::enum_<SeqType>(m, "SeqType")
        .value("AminoAcid",   SeqType::AminoAcid)
        .value("Nucleotide",  SeqType::Nucleotide)
        .value("HMM",         SeqType::HMM);

    py::class_<FoldSeekOps>(m, "FoldSeekOps")
        .def(py::init<>())
        .def_readwrite("alignType",               &FoldSeekOps::alignType)
        .def_readwrite("tmScoreThr",              &FoldSeekOps::tmScoreThr)
        .def_readwrite("tmScoreThrMode",          &FoldSeekOps::tmScoreThrMode)
        .def_readwrite("exactTMscore",            &FoldSeekOps::exactTMscore)
        .def_readwrite("lddtThr",                 &FoldSeekOps::lddtThr)
        .def_readwrite("sortByStructureBits",     &FoldSeekOps::sortByStructureBits)
        .def_readwrite("alignmentMode",           &FoldSeekOps::alignmentMode)
        .def_readwrite("alignmentOutputMode",     &FoldSeekOps::alignmentOutputMode)
        .def_readwrite("wrappedScoring",          &FoldSeekOps::wrappedScoring)
        .def_readwrite("maxSeqLen",               &FoldSeekOps::maxSeqLen)
        .def_readwrite("compBiasCorrection",      &FoldSeekOps::compBiasCorrection)
        .def_readwrite("compBiasCorrectionScale", &FoldSeekOps::compBiasCorrectionScale)
        .def_readwrite("scoreBias",               &FoldSeekOps::scoreBias)
        .def_readwrite("realign",                 &FoldSeekOps::realign)
        .def_readwrite("correlationScoreWeight",  &FoldSeekOps::correlationScoreWeight)
        .def_readwrite("addBacktrace",            &FoldSeekOps::addBacktrace)
        .def_readwrite("covThr",                  &FoldSeekOps::covThr)
        .def_readwrite("covMode",                 &FoldSeekOps::covMode)
        .def_readwrite("evalThr",                 &FoldSeekOps::evalThr)
        .def_readwrite("seqIdThr",                &FoldSeekOps::seqIdThr)
        .def_readwrite("seqIdMode",               &FoldSeekOps::seqIdMode)
        .def_readwrite("alnLenThr",               &FoldSeekOps::alnLenThr)
        .def_readwrite("chainNameMode",           &FoldSeekOps::chainNameMode)
        .def_readwrite("maskBfactorThreshold",    &FoldSeekOps::maskBfactorThreshold)
        .def_readwrite("altAlignment",            &FoldSeekOps::altAlignment)
        .def_readwrite("inputFormat",             &FoldSeekOps::inputFormat)
        .def_readwrite("gapOpen",                 &FoldSeekOps::gapOpen)
        .def_readwrite("gapExtend",               &FoldSeekOps::gapExtend);

    py::class_<PrefilterOptions>(m, "PrefilterOptions")
        .def(py::init<>())
        .def_readwrite("use_prefilter",         &PrefilterOptions::use_prefilter)
        .def_readwrite("alphabetSize",          &PrefilterOptions::alphabetSize)
        .def_readwrite("maskMode",              &PrefilterOptions::maskMode)
        .def_readwrite("maskLowerCaseMode",     &PrefilterOptions::maskLowerCaseMode)
        .def_readwrite("maskProb",              &PrefilterOptions::maskProb)
        .def_readwrite("kmerSize",              &PrefilterOptions::kmerSize)
        .def_readwrite("kmerThr",               &PrefilterOptions::kmerThr)
        .def_readwrite("spacedKmer",            &PrefilterOptions::spacedKmer)
        .def_readwrite("spacedKmerPattern",     &PrefilterOptions::spacedKmerPattern)
        .def_readwrite("takeOnlyBestKmer",      &PrefilterOptions::takeOnlyBestKmer)
        .def_readwrite("querySeqType",          &PrefilterOptions::querySeqType)
        .def_readwrite("targetSeqType",         &PrefilterOptions::targetSeqType)
        .def_readwrite("targetSearchMode",      &PrefilterOptions::targetSearchMode)
        .def_readwrite("sensitivity",           &PrefilterOptions::sensitivity)
        .def_readwrite("maxSeqLen",             &PrefilterOptions::maxSeqLen)
        .def_readwrite("diagonalScoring",       &PrefilterOptions::diagonalScoring)
        .def_readwrite("minDiagScoreThr",       &PrefilterOptions::minDiagScoreThr)
        .def_readwrite("aaBiasCorrection",      &PrefilterOptions::aaBiasCorrection)
        .def_readwrite("aaBiasCorrectionScale", &PrefilterOptions::aaBiasCorrectionScale)
        .def_readwrite("covThr",                &PrefilterOptions::covThr)
        .def_readwrite("covMode",               &PrefilterOptions::covMode)
        .def_readwrite("maxResListLen",         &PrefilterOptions::maxResListLen);

    py::class_<Matcher::result_t>(m, "MatcherResult")
        .def(py::init<>())
        .def(py::init<unsigned int, int, float, float, float, double,
                      unsigned int, int, int, unsigned int, int, int,
                      unsigned int, int, int, int, int, std::string>(),
             py::arg("dbkey"),    py::arg("score"), py::arg("qcov"),       py::arg("dbcov"),
             py::arg("seqId"),    py::arg("eval"),  py::arg("alnLength"),  py::arg("qStartPos"),
             py::arg("qEndPos"),  py::arg("qLen"),  py::arg("dbStartPos"), py::arg("dbEndPos"),
             py::arg("dbLen"),    py::arg("queryOrfStartPos"), py::arg("queryOrfEndPos"),
             py::arg("dbOrfStartPos"), py::arg("dbOrfEndPos"), py::arg("backtrace"))
        .def(py::init<unsigned int, int, float, float, float, double,
                      unsigned int, int, int, unsigned int, int, int,
                      unsigned int, std::string>(),
             py::arg("dbkey"),    py::arg("score"), py::arg("qcov"),       py::arg("dbcov"),
             py::arg("seqId"),    py::arg("eval"),  py::arg("alnLength"),  py::arg("qStartPos"),
             py::arg("qEndPos"),  py::arg("qLen"),  py::arg("dbStartPos"), py::arg("dbEndPos"),
             py::arg("dbLen"),    py::arg("backtrace"))

        .def_readwrite("dbKey",             &Matcher::result_t::dbKey)
        .def_readwrite("score",             &Matcher::result_t::score)
        .def_readwrite("qcov",              &Matcher::result_t::qcov)
        .def_readwrite("dbcov",             &Matcher::result_t::dbcov)
        .def_readwrite("seqId",             &Matcher::result_t::seqId)
        .def_readwrite("eval",              &Matcher::result_t::eval)
        .def_readwrite("alnLength",         &Matcher::result_t::alnLength)
        .def_readwrite("qStartPos",         &Matcher::result_t::qStartPos)
        .def_readwrite("qEndPos",           &Matcher::result_t::qEndPos)
        .def_readwrite("qLen",              &Matcher::result_t::qLen)
        .def_readwrite("dbStartPos",        &Matcher::result_t::dbStartPos)
        .def_readwrite("dbEndPos",          &Matcher::result_t::dbEndPos)
        .def_readwrite("dbLen",             &Matcher::result_t::dbLen)
        .def_readwrite("queryOrfStartPos",  &Matcher::result_t::queryOrfStartPos)
        .def_readwrite("queryOrfEndPos",    &Matcher::result_t::queryOrfEndPos)
        .def_readwrite("dbOrfStartPos",     &Matcher::result_t::dbOrfStartPos)
        .def_readwrite("dbOrfEndPos",       &Matcher::result_t::dbOrfEndPos)
        .def_readwrite("backtrace",         &Matcher::result_t::backtrace)
        .def_static("protein2nucl", &Matcher::result_t::protein2nucl, py::arg("backtrace"), py::arg("newBacktrace"));

    py::class_<Matcher>(m, "Matcher")
        .def_static("compareHits",          &Matcher::compareHits,          py::arg("first"), py::arg("second"))
        .def_static("compressAlignment",    &Matcher::compressAlignment,    py::arg("bt"))
        .def_static("uncompressAlignment",  &Matcher::uncompressAlignment,  py::arg("cbt"))
        .def_static("result_to_string",     &Matcher::result_to_string,     py::arg("r"),  py::arg("compress_backtrace") = true)
        .def_static("results_to_string",    &Matcher::results_to_string,    py::arg("vr"), py::arg("compress_backtrace") = true)
        .def_static("computeAlnLength",     &Matcher::computeAlnLength,
                    py::arg("qStart"),  py::arg("qEnd"),
                    py::arg("dbStart"), py::arg("dbEnd"))
        .def_static("resultToBuffer", &result_to_buffer_wrapper,
                    py::arg("result"),          py::arg("addBacktrace"),
                    py::arg("compress") = true, py::arg("addOrfPosition") = false);
}

} // namespace lahuta::bindings
