#ifndef LAHUTA_DB_WRITER_HPP
#define LAHUTA_DB_WRITER_HPP

#include "db/data_type.hpp"
#include "lmdb/lmdb++.h"
#include "models/parser.hpp"

// clang-format off
namespace lahuta {

class LMDBWriter {
public:
    LMDBWriter(lmdb::env &env, lmdb::dbi &dbi)
        : m_env(env), m_dbi(dbi), m_txn(nullptr) {}

    ~LMDBWriter() {
        if (m_txn) {
            m_txn.abort();
        }
    }

    void begin_txn() {
        if (!m_txn) {
            m_txn = lmdb::txn::begin(m_env.handle());
        }
    }

    void commit_txn() {
        if (m_txn) {
            m_txn.commit();
            m_txn = nullptr;
        }
    }

    void abort_txn() {
        if (m_txn) {
            m_txn.abort();
            m_txn = nullptr;
        }
    }

    bool store(const std::string &key, const ModelParserResult &data, bool commit = true) {
        size_t sequence_size  = data.sequence.size();
        size_t taxonomy_size  = data.metadata.ncbi_taxonomy_id.size();
        size_t organism_size  = data.metadata.organism_scientific.size();
        size_t coords_size    = data.coords.size() * 3 * sizeof(float); // 3 floats per point
        const uint32_t plddt_len = static_cast<uint32_t>(data.plddt_per_residue.size());
        const size_t plddt_bytes = static_cast<size_t>(plddt_len) * sizeof(pLDDTCategory);
        size_t total_size     = sizeof(SerializedModelData) + sequence_size + taxonomy_size + organism_size + coords_size +
                                sizeof(uint32_t) + plddt_bytes;

        // Allocate a buffer and serialize the header, sequence, and coordinates
        std::unique_ptr<char[]> buffer(new char[total_size]);
        auto* serialized = reinterpret_cast<SerializedModelData*>(buffer.get());
        serialized->sequence_length = static_cast<uint32_t>(sequence_size);
        serialized->num_points = static_cast<uint32_t>(data.coords.size());
        serialized->ncbi_taxonomy_id_length = static_cast<uint32_t>(taxonomy_size);
        serialized->organism_scientific_length = static_cast<uint32_t>(organism_size);

        char* seq_ptr = buffer.get() + sizeof(SerializedModelData);
        if (sequence_size) std::memcpy(seq_ptr, data.sequence.data(), sequence_size);

        char* taxonomy_ptr = seq_ptr + sequence_size;
        if (taxonomy_size) std::memcpy(taxonomy_ptr, data.metadata.ncbi_taxonomy_id.data(), taxonomy_size);

        char* organism_ptr = taxonomy_ptr + taxonomy_size;
        if (organism_size) std::memcpy(organism_ptr, data.metadata.organism_scientific.data(), organism_size);

        // convert Point3D to float
        float* float_coords = reinterpret_cast<float*>(organism_ptr + organism_size);
        for (size_t i = 0; i < data.coords.size(); ++i) {
            float_coords[i * 3]     = static_cast<float>(data.coords[i].x);
            float_coords[i * 3 + 1] = static_cast<float>(data.coords[i].y);
            float_coords[i * 3 + 2] = static_cast<float>(data.coords[i].z);
        }

        char* plddt_ptr = reinterpret_cast<char*>(float_coords + data.coords.size() * 3);
        std::memcpy(plddt_ptr, &plddt_len, sizeof(plddt_len));
        plddt_ptr += sizeof(plddt_len);
        if (plddt_bytes) {
            std::memcpy(plddt_ptr, data.plddt_per_residue.data(), plddt_bytes);
        }

        bool created_txn = false;
        if (!m_txn) {
            m_txn = lmdb::txn::begin(m_env.handle());
            created_txn = true;
        }

        std::string_view value(buffer.get(), total_size);
        bool success = m_dbi.put(m_txn.handle(), key, value);

        // if requested or if we created a txn just for this operation
        if (commit || created_txn) {
            m_txn.commit();
            m_txn = nullptr;
        }

        return success;
    }

    bool put_raw(const std::string& key, std::string_view   val, bool commit_now /* default=false */) {
        if (!m_txn) m_txn = lmdb::txn::begin(m_env.handle());

        bool ok = m_dbi.put(m_txn.handle(), key, val);
        if (commit_now) commit_txn();
        return ok;
    }

private:
    lmdb::env &m_env;
    lmdb::dbi &m_dbi;
    lmdb::txn m_txn;
};

} // namespace lahuta

#endif // LAHUTA_DB_WRITER_HPP
