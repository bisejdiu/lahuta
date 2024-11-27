// Written by Martin Steinegger martin.steinegger@snu.ac.kr
//
// Represents parameters of MMseqs2
//
#ifndef MMSEQS_PARAMETERS
#define MMSEQS_PARAMETERS
#include <string>
#include <vector>
#include <map>
#include <typeinfo>
#include <cstddef>
#include <utility>
#include <cstdint>

class Parameters {
public:

    static const int DBTYPE_AMINO_ACIDS = 0;
    static const int DBTYPE_NUCLEOTIDES = 1;
    static const int DBTYPE_HMM_PROFILE = 2;
    //static const int DBTYPE_PROFILE_STATE_SEQ = 3;
    //static const int DBTYPE_PROFILE_STATE_PROFILE = 4;
    static const int DBTYPE_ALIGNMENT_RES = 5;
    static const int DBTYPE_CLUSTER_RES = 6;
    static const int DBTYPE_PREFILTER_RES = 7;
    static const int DBTYPE_TAXONOMICAL_RESULT = 8;
    static const int DBTYPE_INDEX_DB = 9;
    static const int DBTYPE_CA3M_DB = 10;
    static const int DBTYPE_MSA_DB = 11;
    static const int DBTYPE_GENERIC_DB = 12;
    static const int DBTYPE_OMIT_FILE = 13;
    static const int DBTYPE_PREFILTER_REV_RES = 14;
    static const int DBTYPE_OFFSETDB = 15;
    static const int DBTYPE_DIRECTORY = 16; // needed for verification
    static const int DBTYPE_FLATFILE = 17; // needed for verification
    static const int DBTYPE_SEQTAXDB = 18; // needed for verification
    static const int DBTYPE_STDIN = 19; // needed for verification
    static const int DBTYPE_URI = 20; // needed for verification

    static const unsigned int DBTYPE_EXTENDED_COMPRESSED = 1;
    static const unsigned int DBTYPE_EXTENDED_INDEX_NEED_SRC = 2;
    static const unsigned int DBTYPE_EXTENDED_CONTEXT_PSEUDO_COUNTS = 4;

    static const int SEARCH_TYPE_AUTO = 0;
    static const int SEARCH_TYPE_PROTEIN = 1;
    static const int SEARCH_TYPE_TRANSLATED = 2;
    static const int SEARCH_TYPE_NUCLEOTIDES = 3;
    static const int SEARCH_TYPE_TRANS_NUCL_ALN = 4;
    // flag
    static const int SEARCH_MODE_FLAG_QUERY_AMINOACID = 1;
    static const int SEARCH_MODE_FLAG_TARGET_AMINOACID = 2;
    static const int SEARCH_MODE_FLAG_QUERY_TRANSLATED = 4;
    static const int SEARCH_MODE_FLAG_TARGET_TRANSLATED = 8;
    static const int SEARCH_MODE_FLAG_QUERY_PROFILE = 16;
    static const int SEARCH_MODE_FLAG_TARGET_PROFILE = 32;
    static const int SEARCH_MODE_FLAG_QUERY_NUCLEOTIDE = 64;
    static const int SEARCH_MODE_FLAG_TARGET_NUCLEOTIDE = 128;

    static const unsigned int ALIGNMENT_MODE_FAST_AUTO = 0;
    static const unsigned int ALIGNMENT_MODE_SCORE_ONLY = 1;
    static const unsigned int ALIGNMENT_MODE_SCORE_COV = 2;
    static const unsigned int ALIGNMENT_MODE_SCORE_COV_SEQID = 3;
    static const unsigned int ALIGNMENT_MODE_UNGAPPED = 4;

    static const unsigned int ALIGNMENT_OUTPUT_ALIGNMENT = 0;
    static const unsigned int ALIGNMENT_OUTPUT_CLUSTER = 1;

    static const unsigned int EXPAND_TRANSFER_EVALUE = 0;
    static const unsigned int EXPAND_RESCORE_BACKTRACE = 1;

    static const unsigned int PCMODE_SUBSTITUTION_SCORE = 0;
    static const unsigned int PCMODE_CONTEXT_SPECIFIC = 1;


    static const unsigned int WRITER_ASCII_MODE = 0;
    static const unsigned int WRITER_COMPRESSED_MODE = 1;
    static const unsigned int WRITER_LEXICOGRAPHIC_MODE = 2;

    // convertalis alignment
    static const int FORMAT_ALIGNMENT_BLAST_TAB = 0;
    static const int FORMAT_ALIGNMENT_SAM = 1;
    static const int FORMAT_ALIGNMENT_BLAST_WITH_LEN = 2;
    static const int FORMAT_ALIGNMENT_HTML = 3;
    static const int FORMAT_ALIGNMENT_BLAST_TAB_WITH_HEADERS = 4;

    // result2msa
    static const int FORMAT_MSA_CA3M = 0;
    static const int FORMAT_MSA_CA3M_CONSENSUS = 1;
    static const int FORMAT_MSA_FASTADB = 2;
    static const int FORMAT_MSA_FASTADB_SUMMARY = 3;
    static const int FORMAT_MSA_STOCKHOLM_FLAT = 4;
    static const int FORMAT_MSA_A3M = 5;
    static const int FORMAT_MSA_A3M_ALN_INFO = 6;
    // outfmt
    static const int OUTFMT_QUERY = 0;
    static const int OUTFMT_TARGET = 1;
    static const int OUTFMT_EVALUE = 2;
    static const int OUTFMT_GAPOPEN = 3;
    static const int OUTFMT_PIDENT = 4;
    static const int OUTFMT_NIDENT = 5;
    static const int OUTFMT_QSTART = 6;
    static const int OUTFMT_QEND = 7;
    static const int OUTFMT_QLEN = 8;
    static const int OUTFMT_TSTART = 9;
    static const int OUTFMT_TEND = 10;
    static const int OUTFMT_TLEN = 11;
    static const int OUTFMT_ALNLEN = 12;
    static const int OUTFMT_RAW = 13;
    static const int OUTFMT_BITS = 14;
    static const int OUTFMT_CIGAR = 15;
    static const int OUTFMT_QSEQ = 16;
    static const int OUTFMT_TSEQ = 17;
    static const int OUTFMT_QHEADER = 18;
    static const int OUTFMT_THEADER = 19;
    static const int OUTFMT_QALN = 20;
    static const int OUTFMT_TALN = 21;
    static const int OUTFMT_QFRAME = 22;
    static const int OUTFMT_TFRAME = 23;
    static const int OUTFMT_MISMATCH = 24;
    static const int OUTFMT_QCOV = 25;
    static const int OUTFMT_TCOV = 26;
    static const int OUTFMT_EMPTY = 27;
    static const int OUTFMT_QSET = 28;
    static const int OUTFMT_QSETID = 29;
    static const int OUTFMT_TSET = 30;
    static const int OUTFMT_TSETID = 31;
    static const int OUTFMT_TAXID = 32;
    static const int OUTFMT_TAXNAME = 33;
    static const int OUTFMT_TAXLIN = 34;
    static const int OUTFMT_QORFSTART = 35;
    static const int OUTFMT_QORFEND = 36;
    static const int OUTFMT_TORFSTART = 37;
    static const int OUTFMT_TORFEND = 38;
    static const int OUTFMT_FIDENT = 39;
    static const int OUTFMT_PPOS = 40;

    static const int INDEX_SUBSET_NORMAL = 0;
    static const int INDEX_SUBSET_NO_HEADERS = 1;
    static const int INDEX_SUBSET_NO_PREFILTER = 2;
    static const int INDEX_SUBSET_NO_ALIGNMENT = 4;


    // static std::vector<int> getOutputFormat(int formatMode, const std::string &outformat, bool &needSequences, bool &needBacktrace, bool &needFullHeaders,
    //                                         bool &needLookup, bool &needSource, bool &needTaxonomyMapping, bool &needTaxonomy);

    // clustering
    static const int SET_COVER = 0;
    static const int CONNECTED_COMPONENT = 1;
    static const int GREEDY = 2;
    static const int GREEDY_MEM = 3;

    // clustering
    static const int APC_ALIGNMENTSCORE=1;
    static const int APC_SEQID=2;
    // split mode
    static const int TARGET_DB_SPLIT = 0;
    static const int QUERY_DB_SPLIT = 1;
    static const int DETECT_BEST_DB_SPLIT = 2;

    // taxonomy output
    static const int TAXONOMY_OUTPUT_LCA = 0;
    static const int TAXONOMY_OUTPUT_ALIGNMENT = 1;
    static const int TAXONOMY_OUTPUT_BOTH = 2;

    // aggregate taxonomy
    static const int AGG_TAX_UNIFORM = 0;
    static const int AGG_TAX_MINUS_LOG_EVAL = 1;
    static const int AGG_TAX_SCORE = 2;

    // pairaln dummy mode
    static const int PAIRALN_DUMMY_MODE_OFF = 0;
    static const int PAIRALN_DUMMY_MODE_ON = 1;

    // pairaln mode
    static const int PAIRALN_MODE_ALL_PER_SPECIES = 0;
    static const int PAIRALN_MODE_COVER_ALL_CHAINS = 1;

    // taxonomy search strategy
    static const int TAXONOMY_SINGLE_SEARCH = 1;
    static const int TAXONOMY_2BLCA = 2;
    static const int TAXONOMY_APPROX_2BLCA = 3;
    static const int TAXONOMY_TOP_HIT = 4;

    static const int PARSE_VARIADIC = 1;
    static const int PARSE_REST = 2;
    static const int PARSE_ALLOW_EMPTY = 4;

    // preload mode
    static const int PRELOAD_MODE_AUTO = 0;
    static const int PRELOAD_MODE_FREAD = 1;
    static const int PRELOAD_MODE_MMAP = 2;
    static const int PRELOAD_MODE_MMAP_TOUCH = 3;

    static std::string getSplitModeName(int splitMode) {
        switch (splitMode) {
            case 0: return "Target";
            case 1: return "Query";
            case 2: return "Auto";
            default: return "Error";
        }
    };

    // split
    static const int AUTO_SPLIT_DETECTION = 0;

    static const int MAX_SEQ_LEN = 65535;

    // extractalignedregion
    static const int EXTRACT_QUERY  = 1;
    static const int EXTRACT_TARGET = 2;

    static const int CLUST_HASH_DEFAULT_ALPH_SIZE = 3;
    static const int CLUST_HASH_DEFAULT_MIN_SEQ_ID = 99;
    static const int CLUST_LINEAR_DEFAULT_ALPH_SIZE = 13;
    static const int CLUST_LINEAR_DEFAULT_K = 0;
    static const int CLUST_LINEAR_KMER_PER_SEQ = 0;

    // cov mode
    static const int COV_MODE_BIDIRECTIONAL  = 0;
    static const int COV_MODE_TARGET = 1;
    static const int COV_MODE_QUERY = 2;
    static const int COV_MODE_LENGTH_QUERY = 3;
    static const int COV_MODE_LENGTH_TARGET = 4;
    static const int COV_MODE_LENGTH_SHORTER = 5;

    // seq. id mode
    static const int SEQ_ID_ALN_LEN  = 0;
    static const int SEQ_ID_SHORT = 1;
    static const int SEQ_ID_LONG = 2;

    // seq. split mode
    static const int SEQUENCE_SPLIT_MODE_HARD = 0;
    static const int SEQUENCE_SPLIT_MODE_SOFT = 1;

    // rescorediagonal
    static const int RESCORE_MODE_HAMMING = 0;
    static const int RESCORE_MODE_SUBSTITUTION = 1;
    static const int RESCORE_MODE_ALIGNMENT = 2;
    static const int RESCORE_MODE_END_TO_END_ALIGNMENT = 3;
    static const int RESCORE_MODE_WINDOW_QUALITY_ALIGNMENT = 4;

    // combinepvalperset
    static const int AGGREGATION_MODE_MULTIHIT = 0;
    static const int AGGREGATION_MODE_MIN_PVAL = 1;
    static const int AGGREGATION_MODE_PRODUCT = 2;
    static const int AGGREGATION_MODE_TRUNCATED_PRODUCT = 3;

    // header type
    static const int HEADER_TYPE_UNICLUST = 1;
    static const int HEADER_TYPE_METACLUST = 2;

    // createsubdb, filtertaxseqdb type
    static const int SUBDB_MODE_HARD = 0;
    static const int SUBDB_MODE_SOFT = 1;

    static const int ID_MODE_KEYS = 0;
    static const int ID_MODE_LOOKUP = 1;

    // prefilter mode
    static const int PREF_MODE_KMER = 0;
    static const int PREF_MODE_UNGAPPED = 1;
    static const int PREF_MODE_EXHAUSTIVE = 2;

    // unpackdb
    static const int UNPACK_NAME_KEY = 0;
    static const int UNPACK_NAME_ACCESSION = 1;

    // result direction
    static const int PARAM_RESULT_DIRECTION_QUERY  = 0;
    static const int PARAM_RESULT_DIRECTION_TARGET = 1;

    // MultiParam<NuclAA<std::string>> scoringMatrixFile;       // path to scoring matrix
    // MultiParam<NuclAA<std::string>> seedScoringMatrixFile;   // seed sub. matrix

    // static std::vector<std::string> findMissingTaxDbFiles(const std::string &filename);
    // static void printTaxDbError(const std::string &filename, const std::vector<std::string>& missingFiles);

    static const uint32_t DBTYPE_MASK = 0x0000FFFF; // 0x0000FFFF = 65535

    // isEqualDbtype works by comparing the lower 16 bits of the dbtype
    static bool isEqualDbtype(const int type1, const int type2) {
        return ((type1 & DBTYPE_MASK) == (type2 & DBTYPE_MASK));
    }


    template <typename T>
    typename std::enable_if<
    std::is_enum<T>::value, void>::type
    static _isEqualDbtype(T value) {
        int intValue = static_cast<int>(value);
        return isEqualDbtype(intValue, intValue);
    };

    static const char* getDbTypeName(int dbtype) {
        switch (dbtype & DBTYPE_MASK) {
            case DBTYPE_AMINO_ACIDS: return "Aminoacid";
            case DBTYPE_NUCLEOTIDES: return "Nucleotide";
            case DBTYPE_HMM_PROFILE: return "Profile";
            case DBTYPE_ALIGNMENT_RES: return "Alignment";
            case DBTYPE_CLUSTER_RES: return "Clustering";
            case DBTYPE_PREFILTER_RES: return "Prefilter";
            case DBTYPE_TAXONOMICAL_RESULT: return "Taxonomy";
            case DBTYPE_INDEX_DB: return "Index";
            case DBTYPE_CA3M_DB: return "CA3M";
            case DBTYPE_MSA_DB: return "MSA";
            case DBTYPE_GENERIC_DB: return "Generic";
            case DBTYPE_PREFILTER_REV_RES: return "Bi-directional prefilter";
            case DBTYPE_OFFSETDB: return "Offsetted headers";
            case DBTYPE_DIRECTORY: return "Directory";
            case DBTYPE_FLATFILE: return "Flatfile";
            case DBTYPE_STDIN: return "stdin";
            case DBTYPE_URI: return "uri";

            default: return "Unknown";
        }
    }

protected:
    Parameters();
    static Parameters* instance;

private:
    Parameters(Parameters const&);
    void operator=(Parameters const&);

};

#endif
