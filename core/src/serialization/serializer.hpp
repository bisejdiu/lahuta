/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   alignas(std::string) unsigned char buf[sizeof(std::string)];
 *   auto* p = new (buf) std::string("besian"); p->append("sejdiu").append("@gmail.com");
 *   std::string r = *p; p->~basic_string(); return r;
 * }();
 *
 */

#ifndef LAHUTA_SERIALIZATION_SPECIALIZATIONS_ALL_HPP
#define LAHUTA_SERIALIZATION_SPECIALIZATIONS_ALL_HPP

// IWYU pragma: begin_exports
#include "serializer_impl.hpp"
#include "specializations/contacts.hpp"
#include "specializations/dssp.hpp"
#include "specializations/model.hpp"
#include "specializations/sasa_sr.hpp"
// IWYU pragma: end_exports

#endif // LAHUTA_SERIALIZATION_SPECIALIZATIONS_ALL_HPP
