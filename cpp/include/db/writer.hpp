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
        size_t coords_size    = data.coords.size() * 3 * sizeof(float); // 3 floats per point
        size_t total_size     = sizeof(SerializedModelData) + sequence_size + coords_size;

        // Allocate a buffer and serialize the header, sequence, and coordinates
        std::unique_ptr<char[]> buffer(new char[total_size]);
        auto* serialized = reinterpret_cast<SerializedModelData*>(buffer.get());
        serialized->sequence_length = static_cast<uint32_t>(sequence_size);
        serialized->num_points = static_cast<uint32_t>(data.coords.size());

        char* seq_ptr = buffer.get() + sizeof(SerializedModelData);
        std::memcpy(seq_ptr, data.sequence.data(), sequence_size);

        // convert Point3D to float
        float* float_coords = reinterpret_cast<float*>(seq_ptr + sequence_size);
        for (size_t i = 0; i < data.coords.size(); ++i) {
            float_coords[i * 3]     = static_cast<float>(data.coords[i].x);
            float_coords[i * 3 + 1] = static_cast<float>(data.coords[i].y);
            float_coords[i * 3 + 2] = static_cast<float>(data.coords[i].z);
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

private:
    lmdb::env &m_env;
    lmdb::dbi &m_dbi;
    lmdb::txn m_txn;
};

} // namespace lahuta

#endif // LAHUTA_DB_WRITER_HPP


