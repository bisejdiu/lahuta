#ifndef LAHUTA_PIPELINE_DSL_PRIMITIVES_SOURCES_HPP
#define LAHUTA_PIPELINE_DSL_PRIMITIVES_SOURCES_HPP

#include <string>
#include <vector>

#include "pipeline/ingest/db_keys.hpp"
#include "pipeline/ingest/directory.hpp"
#include "pipeline/ingest/file_list.hpp"
#include "pipeline/ingest/file_reader.hpp"
#include "pipeline/ingest/in_memory.hpp"

// clang-format off
namespace lahuta::pipeline {

template<typename>   struct is_source : std::false_type {};
template<typename S> struct source_t { S src; };
template<typename S> struct is_source<source_t<S>> : std::true_type {};
template<typename Src> source_t(Src) -> source_t<Src>;

inline auto directory(std::string root, std::string ext, bool recursive=true) {
    return source_t{ Directory{ std::move(root), std::move(ext), recursive } };
}

inline auto vector_of(std::vector<std::string> v) {
    return source_t{ InMemory<std::string>{ std::move(v) } };
}
inline auto db_keys(LMDBDatabase& db, std::size_t batch) {
    return source_t{ DBKeys{ db, batch } };
}
inline auto file_list(std::string file) {
    return source_t{ FileList{ std::move(file) } };
}

template<typename T>
inline auto test_data(std::vector<T> data, std::size_t batch_size = 1024) {
    return source_t{ InMemory<T>{ std::move(data), batch_size } };
}

template<typename FormatTag, typename T>
inline auto read_serialized_file(std::string file_path, std::size_t batch_size = 1024) {
    return source_t{ FileReader<FormatTag, T>{ file_path, batch_size } };
}

/// recursively find the source in a pipeline
template<typename P>
auto& find_source(P& p) {
  if constexpr (is_source<std::decay_t<P>>::value)
    return p.src;
  else
    return find_source(p.lhs);
}

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_DSL_PRIMITIVES_SOURCES_HPP
