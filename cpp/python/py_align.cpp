#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "aligner.hpp"
#include "array.hpp"
#include "ops.hpp"
#include "prefilter.hpp"
#include "array.hpp"
#include "procs/search_and_align.hpp"
#include "py_align.hpp"

// clang-format off
namespace py = pybind11;
using namespace lahuta;

// A wrapper that uses a std::vector<char> as a buffer
static std::string resultToBufferWrapper(const Matcher::result_t &result, bool addBacktrace, bool compress, bool addOrfPosition) {
    std::vector<char> buff(1024);
    size_t len = Matcher::resultToBuffer(buff.data(), result, addBacktrace, compress, addOrfPosition);
    return std::string(buff.data(), len);
}


void bind_align(py::module &_lahuta) {

    py::enum_<AlignType>      align_type      (_lahuta, "AlignType");
    py::enum_<TMScoreThrMode> tmscore_thr_mode(_lahuta, "TMScoreThrMode");
    py::enum_<SeqType>        seq_type        (_lahuta, "SeqType");

    py::class_<FoldSeekOps>      foldseekopts     (_lahuta, "FoldSeekOps");
    py::class_<PrefilterOptions> prefilteropts    (_lahuta, "PrefilterOptions");
    py::class_<ProcessingConfig> processing_config(_lahuta, "ProcessingConfig");

    py::class_<Matcher> matcher(_lahuta, "Matcher");
    py::class_<Matcher::result_t> matcher_result(_lahuta, "MatcherResult");

    py::class_<SeqData,            std::shared_ptr<SeqData>>            seq_data    (_lahuta, "SeqData");
    py::class_<LahutaAlignerBase,  std::shared_ptr<LahutaAlignerBase>>  aligner_base(_lahuta, "LahutaAlignerBase");
    py::class_<LahutaProcessor,    std::shared_ptr<LahutaProcessor>>    processor   (_lahuta, "LahutaProcessor");
    py::class_<AlignerResults,     std::shared_ptr<AlignerResults>>     alig_result (_lahuta, "AlignerResults");

    py::class_<LahutaAligner, LahutaAlignerBase,  std::shared_ptr<LahutaAligner>>  aligner(_lahuta, "LahutaAligner");


    align_type
        .value("AA_3Di",  AlignType::AA_3Di)
        .value("AA",      AlignType::AA)
        .value("_3Di",    AlignType::_3Di);

    tmscore_thr_mode
        .value("alignment", TMScoreThrMode::alignment)
        .value("query",     TMScoreThrMode::query)
        .value("target",    TMScoreThrMode::target)
        .value("min",       TMScoreThrMode::min);

    seq_type
        .value("AminoAcid",   SeqType::AminoAcid)
        .value("Nucleotide",  SeqType::Nucleotide)
        .value("HMM",         SeqType::HMM);

    foldseekopts
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

    prefilteropts
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

    processing_config
        .def(py::init<>())
        .def_readwrite("query_chunk_size",  &ProcessingConfig::query_chunk_size)
        .def_readwrite("target_chunk_size", &ProcessingConfig::target_chunk_size)
        .def_readwrite("allow_self_ops",    &ProcessingConfig::allow_self_ops);

    aligner
        .def(py::init<FoldSeekOps, PrefilterOptions, unsigned int>(),
             py::arg("ops") = FoldSeekOps(), py::arg("pf_ops") = PrefilterOptions(), py::arg("n_threads") = 0)
        .def("run",         &LahutaAligner::run, py::arg("query_files"), py::arg("target_files"))
        .def("get_results", &LahutaAligner::get_results);

    processor
        .def(py::init<std::shared_ptr<LahutaAlignerBase>, const ProcessingConfig &>(), py::arg("aligner"), py::arg("config"))
        .def("process",
             py::overload_cast<const std::vector<std::string>&, const std::vector<std::string>&>(&LahutaProcessor::process),
             py::arg("query_files"), py::arg("target_files"))
        .def("process",
             py::overload_cast<const std::vector<std::string>&>(&LahutaProcessor::process),
             py::arg("query_files"));

    alig_result
        .def(py::init<>())
        .def_property("query",
             [](AlignerResults &self) { return self.query; },
             [](AlignerResults &self, std::shared_ptr<SeqData> q) { self.query = q; },
             py::keep_alive<1, 2>()) // Keep AlignerResults alive as long as returned query is alive
        .def_property("target",
             [](AlignerResults &self) { return self.target; },
             [](AlignerResults &self, std::shared_ptr<SeqData> t) { self.target = t; },
             py::keep_alive<1, 2>()) // Keep AlignerResults alive as long as returned target is alive
        .def_readwrite("results", &AlignerResults::results);

    matcher_result
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

    matcher
        .def_static("compareHits",          &Matcher::compareHits,          py::arg("first"), py::arg("second"))
        .def_static("compressAlignment",    &Matcher::compressAlignment,    py::arg("bt"))
        .def_static("uncompressAlignment",  &Matcher::uncompressAlignment,  py::arg("cbt"))
        .def_static("result_to_string",     &Matcher::result_to_string,     py::arg("r"),  py::arg("compress_backtrace") = true)
        .def_static("results_to_string",    &Matcher::results_to_string,    py::arg("vr"), py::arg("compress_backtrace") = true)
        .def_static("computeAlnLength",     &Matcher::computeAlnLength,
                    py::arg("qStart"),  py::arg("qEnd"),
                    py::arg("dbStart"), py::arg("dbEnd"))
        .def_static("resultToBuffer", &resultToBufferWrapper,
                    py::arg("result"),          py::arg("addBacktrace"),
                    py::arg("compress") = true, py::arg("addOrfPosition") = false);

    seq_data
        .def(py::init<>())
        .def("size", (int (SeqData::*)() const) &SeqData::size, "Returns the number of residues in the sequence")
        .def_readwrite("Seq3Di",      &SeqData::Seq3Di)
        .def_readwrite("SeqAA",       &SeqData::SeqAA)
        .def_readwrite("CaData",      &SeqData::CaData)
        .def_readwrite("file_name",   &SeqData::file_name)
        .def_readwrite("chain_name",  &SeqData::chain_name)
        .def_property_readonly("x", [](SeqData &self) {
            size_t n = self.SeqAA.size();
            auto vec = std::vector<float>(self.CaData.begin(), self.CaData.begin() + n);
            return float_array(vec);
        })
        .def_property_readonly("y", [](SeqData &self) {
            size_t n = self.SeqAA.size();
            auto vec = std::vector<float>(self.CaData.begin() + n, self.CaData.begin() + 2 * n);
            return float_array(vec);
        })
        .def_property_readonly("z", [](SeqData &self) {
            size_t n = self.SeqAA.size();
            auto vec = std::vector<float>(self.CaData.begin() + 2 * n, self.CaData.begin() + 3 * n);
            return float_array(vec);
        })
        .def("__lt__", [](const SeqData &a, const SeqData &b) { return a < b; });
        /*.def("map_3di", &SeqData::map_3di, py::arg("matrix"), py::arg("ops"),*/
        /*     "Maps the 3Di sequence to a Sequence object")*/
        /*.def("map_aa", &SeqData::map_aa, py::arg("matrix"), py::arg("ops"),*/
        /*     "Maps the amino acid sequence to a Sequence object")*/
        /*.def("build_sequence", &SeqData::build_sequence, py::arg("matrix"), py::arg("ops"),*/
        /*     "Builds a Sequence object")*/

}

