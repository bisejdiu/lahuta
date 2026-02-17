/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   constexpr std::size_t align = alignof(std::string); alignas(align) char buf[sizeof(std::string)];
 *   auto* p = new (buf) std::string{"besian"}; p->append("sejdiu").append("@gmail.com");
 *   std::string r = *p; p->~basic_string(); return r;
 * }();
 *
 */

#ifndef LAHUTA_SERIALIZATION_CONTACTS_SERIALIZER_HPP
#define LAHUTA_SERIALIZATION_CONTACTS_SERIALIZER_HPP

#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string_view>
#include <vector>

#include "analysis/contacts/provider.hpp"
#include "analysis/contacts/records.hpp"
#include "entities/formatter.hpp"
#include "entities/interaction_types.hpp"
#include "entities/resolver.hpp"
#include "serialization/formats.hpp"
#include "serialization/json.hpp"
#include "serialization/serializer_impl.hpp"

namespace lahuta {

struct JsonNumber {
  std::string_view repr;
};

template <>
struct JsonWritable<JsonNumber> : std::true_type {
  static void write(std::string &out, JsonNumber v) {
    out.append(v.repr.data(), v.repr.size());
  }
};

} // namespace lahuta

namespace serialization {
using namespace lahuta;

using ContactsRes = lahuta::analysis::ContactsRecord;

namespace detail {

inline std::string_view provider_short_code(analysis::ContactProvider provider) noexcept {
  switch (provider) {
    case analysis::ContactProvider::MolStar:     return "ms";
    case analysis::ContactProvider::Arpeggio:    return "ap";
    case analysis::ContactProvider::GetContacts: return "gc";
  }
  return {};
}

inline std::optional<analysis::ContactProvider> provider_from_short_code(std::string_view code) noexcept {
  if (code == "ms") return analysis::ContactProvider::MolStar;
  if (code == "ap") return analysis::ContactProvider::Arpeggio;
  if (code == "gc") return analysis::ContactProvider::GetContacts;
  return std::nullopt;
}

inline std::string format_compact_distance(float distance_sq) {
  const double rounded = std::round(static_cast<double>(distance_sq) * 100.0) / 100.0;
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(2) << rounded;
  std::string out = oss.str();
  if (auto dot = out.find('.'); dot != std::string::npos) {
    while (!out.empty() && out.back() == '0') out.pop_back();
    if (!out.empty() && out.back() == '.') out.pop_back();
  }
  return out;
}

inline uint32_t hash_label(std::string_view label) noexcept {
  uint32_t h = 2166136261u;
  for (unsigned char c : label) {
    h ^= c;
    h *= 16777619u;
  }
  return h;
}

inline uint32_t parse_label_index(std::string_view label, bool &parsed) noexcept {
  parsed = false;
  if (label.empty()) return 0;
  std::size_t i = (label.front() == '(') ? 1u : 0u;
  const std::size_t start = i;
  uint64_t value = 0;
  while (i < label.size()) {
    const char c = label[i];
    if (c < '0' || c > '9') break;
    value = (value * 10u) + static_cast<uint64_t>(c - '0');
    if (value > std::numeric_limits<uint32_t>::max()) {
      parsed = false;
      return 0;
    }
    ++i;
  }
  if (i == start) return 0;
  parsed = true;
  return static_cast<uint32_t>(value);
}

inline EntityID compact_entity_id(std::string_view label) noexcept {
  bool parsed = false;
  uint32_t index = parse_label_index(label, parsed);
  if (!parsed) index = hash_label(label);
  const Kind kind = (!label.empty() && label.front() == '(') ? Kind::Group : Kind::Atom;
  return EntityID::make(kind, index);
}

inline std::string_view sv_from_sajson(const sajson::string &str) noexcept {
  return {str.data(), str.length()};
}

inline std::string require_string(const sajson::value &v, std::string_view key) {
  if (v.get_type() != sajson::TYPE_STRING) {
    throw std::runtime_error("expected string for key \"" + std::string(key) + '"');
  }
  return v.as_string();
}

inline std::string optional_string(const sajson::value &v, std::string_view key) {
  if (v.get_type() == sajson::TYPE_NULL) return {};
  return require_string(v, key);
}

inline bool optional_bool(const sajson::value &v, bool default_value, std::string_view key) {
  if (v.get_type() == sajson::TYPE_NULL) return default_value;
  if (v.get_type() == sajson::TYPE_TRUE) return true;
  if (v.get_type() == sajson::TYPE_FALSE) return false;
  throw std::runtime_error("expected boolean for key \"" + std::string(key) + '"');
}

inline double require_number(const sajson::value &v, std::string_view key) {
  if (v.get_type() == sajson::TYPE_INTEGER) return static_cast<double>(v.get_integer_value());
  if (v.get_type() == sajson::TYPE_DOUBLE)  return v.get_double_value();
  throw std::runtime_error("expected number for key \"" + std::string(key) + '"');
}

inline double optional_number(const sajson::value &v, double default_value, std::string_view key) {
  if (v.get_type() == sajson::TYPE_NULL) return default_value;
  return require_number(v, key);
}

} // namespace detail

template <>
struct Serializer<fmt::json, ContactsRes> {
  using Record = ContactsRes;

  static std::string serialize(const ContactsRes &v) {
    JsonBuilder builder;

    const std::string contact_type_str = interaction_type_set_to_string(v.contact_types, "|");

    builder.key("file_path").value(v.file_path);
    if (v.trajectory_file) {
      builder.key("trajectory_file").value(*v.trajectory_file);
    }
    // clang-format off
    builder.key("success")     .value(v.success)
           .key("provider")    .value(std::string(contact_provider_name(v.provider)))
           .key("contact_type").value(contact_type_str)
           .key("num_contacts").value(v.num_contacts)
           .key("frame_index") .value(v.frame_index);

    if (!v.topology) {
      Logger::get_logger()->warn("ContactsRes serialization: topology is null, cannot serialize contacts.");
      return builder.str();
    }

    builder.key("contacts").begin_array();
    EntityResolver resolver(*v.topology);
    for (const auto &contact : v.contacts) {
      auto pair = resolver.resolve(contact);
      builder.begin_object()
             .key("lhs")     .value(ContactTableFormatter::format_entity_compact(*v.topology, pair.first))
             .key("rhs")     .value(ContactTableFormatter::format_entity_compact(*v.topology, pair.second))
             .key("distance_sq").value(contact.distance)
             .key("type")    .value(interaction_type_to_string(contact.type))
             .end_object();
    }
    builder.end_array();
    return builder.str();
  }

  static ContactsRes deserialize(const std::string &s) {
    ContactsRes out;
    JsonReader r{s};

    out.file_path = r.get<std::string>("file_path");
    if (auto traj = r.get_or<std::string>("trajectory_file", ""); !traj.empty()) {
      out.trajectory_file = traj;
    }
    out.success   = r.get<bool>("success");

    // clang-format off
    std::string provider_str = r.get<std::string>("provider");
    if      (provider_str == "molstar")     out.provider = analysis::ContactProvider::MolStar;
    else if (provider_str == "arpeggio")    out.provider = analysis::ContactProvider::Arpeggio;
    else if (provider_str == "getcontacts") out.provider = analysis::ContactProvider::GetContacts;
    else throw std::runtime_error("Unknown contact provider: " + provider_str);
    // clang-format on

    std::string contact_type_str = r.get<std::string>("contact_type");
    if (auto parsed = parse_interaction_type_sequence(contact_type_str, '|')) {
      out.contact_types = *parsed;
    } else {
      throw std::runtime_error("ContactsRes JSON deserialize: unknown contact type list: " +
                               contact_type_str);
    }

    out.num_contacts = r.get<size_t>("num_contacts");
    out.frame_index  = r.get_or<std::size_t>("frame_index", 0);
    out.topology     = nullptr; // we cannot restore the topology from the serialized contact data

    return out;
  }
};

template <>
struct Serializer<fmt::json_compact, ContactsRes> {
  using Record = ContactsRes;

  static std::string serialize(const ContactsRes &v) {
    JsonBuilder builder;

    builder.key("f").value(v.file_path);
    if (v.trajectory_file) {
      builder.key("tf").value(*v.trajectory_file);
    }
    if (!v.success) {
      builder.key("ok").value(false);
    }
    builder.key("p").value(detail::provider_short_code(v.provider));
    if (!v.contact_types.is_all()) {
      builder.key("ct").value(interaction_type_set_to_string(v.contact_types, "|"));
    }
    if (v.frame_index != 0) {
      builder.key("fi").value(v.frame_index);
    }

    if (!v.topology) {
      Logger::get_logger()->warn(
          "ContactsRes compact serialization: topology is null, cannot serialize contacts.");
      return builder.str();
    }

    auto sorted = sort_interactions(v.contacts);
    auto spans  = slice_by_type(sorted);

    builder.key("c").begin_object();
    EntityResolver resolver(*v.topology);
    for (auto s : spans) {
      if (s.empty()) continue;
      const InteractionType type = s[0].type;
      const std::string_view code = interaction_type_to_short_code(type);
      if (code.empty()) {
        throw std::runtime_error("ContactsRes compact serializer: unknown interaction type");
      }
      builder.key(code).begin_array();
      for (const auto &contact : s) {
        auto pair = resolver.resolve(contact);
        std::string lhs = ContactTableFormatter::format_entity_compact(*v.topology, pair.first, ",");
        std::string rhs = ContactTableFormatter::format_entity_compact(*v.topology, pair.second, ",");
        std::string dist = detail::format_compact_distance(contact.distance);
        builder.begin_array()
            .value(lhs)
            .value(rhs)
            .value(JsonNumber{dist})
            .end_array();
      }
      builder.end_array();
    }
    builder.end_object();

    return builder.str();
  }

  static ContactsRes deserialize(const std::string &s) {
    ContactsRes out{};

    auto doc = sajson::parse(sajson::dynamic_allocation(),
                             sajson::mutable_string_view(s.size(), const_cast<char *>(s.data())));
    if (!doc.is_valid()) {
      throw std::runtime_error(doc.get_error_message_as_cstring());
    }
    const auto root = doc.get_root();
    if (root.get_type() != sajson::TYPE_OBJECT) {
      throw std::runtime_error("expected JSON object");
    }

    auto get_key = [&](std::string_view key) -> sajson::value {
      const sajson::string k{key.data(), key.size()};
      return root.get_value_of_key(k);
    };

    out.file_path = detail::require_string(get_key("f"), "f");
    out.success   = detail::optional_bool(get_key("ok"), true, "ok");

    const std::string provider_code = detail::require_string(get_key("p"), "p");
    if (auto provider = detail::provider_from_short_code(provider_code)) {
      out.provider = *provider;
    } else {
      throw std::runtime_error("Unknown contact provider code: " + provider_code);
    }

    const std::string ct_value = detail::optional_string(get_key("ct"), "ct");
    if (ct_value.empty()) {
      out.contact_types = InteractionTypeSet::all();
    } else if (auto parsed = parse_interaction_type_sequence(ct_value, '|')) {
      out.contact_types = *parsed;
    } else {
      throw std::runtime_error("ContactsRes compact deserialize: unknown contact type list: " +
                               ct_value);
    }

    out.frame_index = static_cast<std::size_t>(detail::optional_number(get_key("fi"), 0.0, "fi"));

    const std::string traj = detail::optional_string(get_key("tf"), "tf");
    if (!traj.empty()) {
      out.trajectory_file = traj;
    }

    out.topology = nullptr;

    const sajson::value contacts_value = get_key("c");
    if (contacts_value.get_type() == sajson::TYPE_OBJECT) {
      std::vector<Contact> contacts;
      std::size_t total = 0;
      const std::size_t groups = contacts_value.get_length();
      for (std::size_t i = 0; i < groups; ++i) {
        const auto arr = contacts_value.get_object_value(i);
        if (arr.get_type() == sajson::TYPE_ARRAY) {
          total += arr.get_length();
        }
      }
      contacts.reserve(total);

      for (std::size_t i = 0; i < groups; ++i) {
        const auto key = contacts_value.get_object_key(i);
        const std::string_view type_code = detail::sv_from_sajson(key);
        auto type = short_code_to_interaction_type(type_code);
        if (!type) {
          throw std::runtime_error("ContactsRes compact deserialize: unknown type code '" +
                                   std::string(type_code) + "'");
        }

        const auto arr = contacts_value.get_object_value(i);
        if (arr.get_type() != sajson::TYPE_ARRAY) {
          throw std::runtime_error("ContactsRes compact deserialize: expected array for type '" +
                                   std::string(type_code) + "'");
        }
        const std::size_t arr_size = arr.get_length();
        for (std::size_t j = 0; j < arr_size; ++j) {
          const auto entry = arr.get_array_element(j);
          if (entry.get_type() != sajson::TYPE_ARRAY) {
            throw std::runtime_error("ContactsRes compact deserialize: expected array entry");
          }
          if (entry.get_length() != 3) {
            throw std::runtime_error("ContactsRes compact deserialize: expected 3-element entry");
          }

          const auto lhs_value = entry.get_array_element(0);
          const auto rhs_value = entry.get_array_element(1);
          const auto dist_value = entry.get_array_element(2);
          const std::string lhs_label = detail::require_string(lhs_value, "lhs");
          const std::string rhs_label = detail::require_string(rhs_value, "rhs");
          const float dist = static_cast<float>(detail::require_number(dist_value, "distance_sq"));

          const EntityID lhs_id = detail::compact_entity_id(lhs_label);
          const EntityID rhs_id = detail::compact_entity_id(rhs_label);
          contacts.emplace_back(lhs_id, rhs_id, dist, *type);
        }
      }

      out.contacts = ContactSet(std::move(contacts), true);
      out.num_contacts = out.contacts.size();
    } else if (contacts_value.get_type() == sajson::TYPE_NULL) {
      out.num_contacts = 0;
    } else {
      throw std::runtime_error("ContactsRes compact deserialize: expected object for key \"c\"");
    }

    return out;
  }
};

template <>
struct Serializer<fmt::text, ContactsRes> {
  using Record = ContactsRes;

  static std::string serialize(const ContactsRes &v) {
    std::ostringstream oss;

    const std::string contact_type_str = interaction_type_set_to_string(v.contact_types, "|");

    oss << (v.success ? "1" : "0") << " " << v.file_path << " ";
    if (v.trajectory_file) {
      oss << "[traj:" << *v.trajectory_file << "] ";
    }
    oss << contact_provider_name(v.provider) << " " << contact_type_str << " " << v.num_contacts << " "
        << v.frame_index << "\n";

    if (!v.topology) return oss.str();

    EntityResolver resolver(*v.topology);
    for (const auto &contact : v.contacts) {
      auto pair = resolver.resolve(contact);
      oss << ContactTableFormatter::format_entity_compact(*v.topology, pair.first) << " -&- "
          << ContactTableFormatter::format_entity_compact(*v.topology, pair.second) << " " << std::fixed
          << std::setprecision(3) << contact.distance << " A^2 (" << interaction_type_to_string(contact.type)
          << ")";
      oss << "\n";
    }

    return oss.str();
  }

  static ContactsRes deserialize(const std::string &data) {
    // We should be able to re-create contacts from plain text but we need to be careful
    // with non-atom entities.
    throw std::runtime_error("Text format deserialization not implemented for ContactsRes");
  }
};

template <>
struct Serializer<fmt::binary, ContactsRes> {
  using Record = ContactsRes;

  static std::string serialize(const ContactsRes &v) {
    constexpr uint8_t VERSION = 1;

    if (!v.topology)
      throw std::runtime_error("ContactsRes binary serializer: topology is required for serialization");

    EntityResolver resolver(*v.topology);
    std::vector<std::string> lhs_names;
    std::vector<std::string> rhs_names;
    lhs_names.reserve(v.contacts.size());
    rhs_names.reserve(v.contacts.size());

    for (const auto &contact : v.contacts) {
      auto pair = resolver.resolve(contact);
      lhs_names.emplace_back(ContactTableFormatter::format_entity_compact(*v.topology, pair.first));
      rhs_names.emplace_back(ContactTableFormatter::format_entity_compact(*v.topology, pair.second));
    }

    const auto checked_len = [](std::size_t len) -> uint32_t {
      if (len > std::numeric_limits<uint32_t>::max()) {
        throw std::runtime_error("ContactsRes binary serializer: length exceeds 32-bit range");
      }
      return static_cast<uint32_t>(len);
    };

    const bool filter_all                = v.contact_types.is_all();
    std::vector<InteractionType> filters = filter_all ? std::vector<InteractionType>{}
                                                      : v.contact_types.members();

    std::size_t size = sizeof(uint8_t) * 4 /*version, success, provider, filter_mode*/ +
                       sizeof(uint32_t) /*filter_count*/ +
                       sizeof(uint32_t) * filters.size() /*filter codes*/ + sizeof(uint32_t) /*path_len*/ +
                       sizeof(uint64_t) /*frame_index*/ + v.file_path.size() +
                       v.contacts.size() * (sizeof(uint64_t) * 2 + sizeof(float) + sizeof(uint32_t)) +
                       sizeof(uint32_t) * 2 * v.contacts.size(); /* name lengths */

    for (std::size_t i = 0; i < lhs_names.size(); ++i) {
      size += lhs_names[i].size() + rhs_names[i].size();
    }

    std::string buffer;
    buffer.reserve(size);

    const auto append_pod = [&buffer](auto value) {
      using T = decltype(value);
      buffer.append(reinterpret_cast<const char *>(&value), sizeof(T));
    };

    buffer.push_back(static_cast<char>(VERSION));
    buffer.push_back(static_cast<char>(v.success ? 1 : 0));
    buffer.push_back(static_cast<char>(static_cast<uint8_t>(v.provider)));

    buffer.push_back(static_cast<char>(filter_all ? 0 : 1));
    const uint32_t filter_count = static_cast<uint32_t>(filters.size());
    append_pod(filter_count);
    for (auto type : filters) {
      const uint32_t type_raw = static_cast<uint32_t>(type);
      append_pod(type_raw);
    }

    const uint32_t path_len = checked_len(v.file_path.size());
    append_pod(path_len);
    buffer.append(v.file_path.data(), v.file_path.size());

    const uint64_t frame_index = static_cast<uint64_t>(v.frame_index);
    append_pod(frame_index);

    const uint32_t num_contacts = checked_len(v.contacts.size());
    append_pod(num_contacts);

    for (std::size_t idx = 0; idx < v.contacts.size(); ++idx) {
      const auto &contact = v.contacts[idx];
      append_pod(contact.lhs.raw);
      append_pod(contact.rhs.raw);
      append_pod(contact.distance);
      const uint32_t type_raw = static_cast<uint32_t>(contact.type);
      append_pod(type_raw);

      const uint32_t lhs_len = checked_len(lhs_names[idx].size());
      const uint32_t rhs_len = checked_len(rhs_names[idx].size());
      append_pod(lhs_len);
      buffer.append(lhs_names[idx].data(), lhs_names[idx].size());
      append_pod(rhs_len);
      buffer.append(rhs_names[idx].data(), rhs_names[idx].size());
    }

    const std::string traj  = v.trajectory_file.value_or("");
    const uint32_t traj_len = checked_len(traj.size());
    append_pod(traj_len);
    if (traj_len > 0) {
      buffer.append(traj.data(), traj.size());
    }

    return buffer;
  }

  static ContactsRes deserialize(const char *data, std::size_t size) {
    constexpr uint8_t VERSION = 1;

    auto require = [size](std::size_t offset, std::size_t need) {
      if (offset + need > size) {
        throw std::runtime_error("ContactsRes binary deserialize: truncated payload");
      }
    };

    auto read_u8 = [&](std::size_t &offset) -> uint8_t {
      require(offset, sizeof(uint8_t));
      uint8_t value  = static_cast<uint8_t>(data[offset]);
      offset        += sizeof(uint8_t);
      return value;
    };

    auto read_u32 = [&](std::size_t &offset) -> uint32_t {
      require(offset, sizeof(uint32_t));
      uint32_t value;
      std::memcpy(&value, data + offset, sizeof(uint32_t));
      offset += sizeof(uint32_t);
      return value;
    };

    auto read_u64 = [&](std::size_t &offset) -> uint64_t {
      require(offset, sizeof(uint64_t));
      uint64_t value;
      std::memcpy(&value, data + offset, sizeof(uint64_t));
      offset += sizeof(uint64_t);
      return value;
    };

    auto read_f32 = [&](std::size_t &offset) -> float {
      require(offset, sizeof(float));
      float value;
      std::memcpy(&value, data + offset, sizeof(float));
      offset += sizeof(float);
      return value;
    };

    auto skip_bytes = [&](std::size_t &offset, std::size_t count) {
      require(offset, count);
      offset += count;
    };

    ContactsRes result{};
    std::size_t offset = 0;

    const uint8_t version = read_u8(offset);
    if (version != VERSION) {
      throw std::runtime_error("ContactsRes binary deserialize: unsupported version");
    }

    result.success         = read_u8(offset) != 0;
    const uint8_t provider = read_u8(offset);
    result.provider        = static_cast<analysis::ContactProvider>(provider);

    const uint8_t filter_mode   = read_u8(offset);
    const uint32_t filter_count = read_u32(offset);
    InteractionTypeSet filters;
    if (filter_mode == 0) {
      filters = InteractionTypeSet::all();
      // make sure we skip any stale codes if present
      for (uint32_t i = 0; i < filter_count; ++i) {
        (void)read_u32(offset);
      }
    } else {
      for (uint32_t i = 0; i < filter_count; ++i) {
        const uint32_t type_raw = read_u32(offset);
        InteractionType type;
        type.category = static_cast<Category>(type_raw & 0xFFFFu);
        type.flavor   = static_cast<Flavor>((type_raw >> 16) & 0xFFFFu);
        filters.add(type);
      }
    }
    result.contact_types = filters;

    const uint32_t path_len = read_u32(offset);
    require(offset, path_len);
    result.file_path.assign(data + offset, path_len);
    offset += path_len;

    result.frame_index = static_cast<std::size_t>(read_u64(offset));

    const uint32_t num_contacts = read_u32(offset);
    result.num_contacts         = num_contacts;

    std::vector<Contact> contacts;
    contacts.reserve(result.num_contacts);

    for (std::size_t i = 0; i < result.num_contacts; ++i) {
      EntityID lhs;
      EntityID rhs;
      lhs.raw                 = read_u64(offset);
      rhs.raw                 = read_u64(offset);
      const float distance    = read_f32(offset);
      const uint32_t type_raw = read_u32(offset);

      InteractionType type;
      type.category = static_cast<Category>(type_raw & 0xFFFFu);
      type.flavor   = static_cast<Flavor>((type_raw >> 16) & 0xFFFFu);

      contacts.emplace_back(lhs, rhs, distance, type);

      const uint32_t lhs_len = read_u32(offset);
      skip_bytes(offset, lhs_len);
      const uint32_t rhs_len = read_u32(offset);
      skip_bytes(offset, rhs_len);
    }

    result.contacts = ContactSet(std::move(contacts), true);
    result.topology = nullptr;

    const uint32_t traj_len = read_u32(offset);
    if (traj_len > 0) {
      require(offset, traj_len);
      result.trajectory_file  = std::string(data + offset, traj_len);
      offset                 += traj_len;
    }

    return result;
  }

  static ContactsRes deserialize(const std::string &buf) { return deserialize(buf.data(), buf.size()); }
};

} // namespace serialization

#endif // LAHUTA_SERIALIZATION_CONTACTS_SERIALIZER_HPP
