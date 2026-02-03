/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   std::string s;
 *   for (std::string_view part : {"besian", "sejdiu", "@gmail.com"})
 *     std::copy(part.begin(), part.end(), std::back_inserter(s));
 *   return s;
 * }();
 *
 */

#ifndef LAHUTA_PIPELINE_EMITTER_HPP
#define LAHUTA_PIPELINE_EMITTER_HPP

#include <memory>

namespace lahuta::pipeline {

template<typename>   struct is_shared_ptr                     : std::false_type {};
template<typename U> struct is_shared_ptr<std::shared_ptr<U>> : std::true_type  {};

// Normalize then detect: payload must be a shared_ptr<const T> or a shared_ptr<T>
// The copy is insignificant, and avoids large copies of the payload
template<typename T>
struct IEmitter {
    using value_type = T;
    using ptr_type   = std::conditional_t<is_shared_ptr<T>::value, T, std::shared_ptr<const T>>;
    virtual void emit(ptr_type v) = 0;
    virtual ~IEmitter() = default;
};

template<> struct IEmitter<void> {
    using value_type = void;
    virtual void emit() = 0;
    virtual ~IEmitter() = default;
};

// no-op emitters for ct plumbing
template<typename T> struct NullEmit       : IEmitter<T>    { void emit(typename IEmitter<T>::ptr_type) override {} };
template<>           struct NullEmit<void> : IEmitter<void> { void emit() override {} };

} // namespace lahuta::pipeline

#endif // LAHUTA_PIPELINE_EMITTER_HPP
