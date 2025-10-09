#ifndef LAHUTA_SERIALIZATION_CONTACTS_SERIALIZER_HPP
#define LAHUTA_SERIALIZATION_CONTACTS_SERIALIZER_HPP

#include <cstring>
#include <limits>
#include <sstream>
#include <vector>

#include "analysis/contacts/provider.hpp"
#include "analysis/contacts/records.hpp"
#include "entities/formatter.hpp"
#include "entities/interaction_types.hpp"
#include "entities/resolver.hpp"
#include "serialization/formats.hpp"
#include "serialization/json.hpp"
#include "serialization/serializer_impl.hpp"

// clang-format off
namespace serialization {
using namespace lahuta;

using ContactsRes = lahuta::analysis::contacts::ContactsRecord;

template<>
struct Serializer<fmt::json, ContactsRes> {
  using Record = ContactsRes;

  static std::string serialize(const ContactsRes &v) {
    JsonBuilder builder;

    std::string contact_type_str = (v.contact_type == InteractionType::All)
      ? "All"
      : interaction_type_to_string(v.contact_type);

    builder.key("file_path")   .value(v.file_path)
           .key("success")     .value(v.success)
           .key("provider")    .value(v.provider == analysis::contacts::ContactProvider::Arpeggio ? "arpeggio" : "molstar")
           .key("contact_type").value(contact_type_str)
           .key("num_contacts").value(v.num_contacts);

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
             .key("distance").value(contact.distance)
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
    out.success   = r.get<bool>("success");

    std::string provider_str = r.get<std::string>("provider");
    out.provider = (provider_str == "arpeggio") ? analysis::contacts::ContactProvider::Arpeggio : analysis::contacts::ContactProvider::MolStar;

    std::string contact_type_str = r.get<std::string>("contact_type");
    out.contact_type = (contact_type_str == "All")
      ? InteractionType::All
      : get_interaction_type(contact_type_str);

    out.num_contacts = r.get<size_t>("num_contacts");
    out.topology = nullptr; // we cannot restore the topology from the serialized contact data

    return out;
  }
};

template<>
struct Serializer<fmt::text, ContactsRes> {
  using Record = ContactsRes;

  static std::string serialize(const ContactsRes &v) {
    std::ostringstream oss;

    std::string contact_type_str = (v.contact_type == InteractionType::All)
      ? "All"
      : interaction_type_to_string(v.contact_type);

    oss << (v.success ? "1" : "0") << " "
        << v.file_path << " "
        << (v.provider == analysis::contacts::ContactProvider::Arpeggio ? "arpeggio" : "molstar") << " "
        << contact_type_str << " "
        << v.num_contacts << "\n";

    if (!v.topology) return oss.str();

    EntityResolver resolver(*v.topology);
    for (const auto &contact : v.contacts) {
      auto pair = resolver.resolve(contact);
      oss << ContactTableFormatter::format_entity_compact(*v.topology, pair.first)
          << " -&- " << ContactTableFormatter::format_entity_compact(*v.topology, pair.second)
          << " " << std::fixed << std::setprecision(3) << contact.distance
          << " (" << interaction_type_to_string(contact.type) << ")";
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

template<>
struct Serializer<fmt::binary, ContactsRes> {
  using Record = ContactsRes;

  static std::string serialize(const ContactsRes &v) {
    constexpr uint8_t VERSION = 1;

    if (!v.topology) throw std::runtime_error("ContactsRes binary serializer: topology is required for serialization");

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

    std::size_t size = sizeof(uint8_t)  * 3  /*version, success, provider*/ +
                       sizeof(uint32_t) * 3  /*contact_type, path_len, num_contacts*/ +
                       v.file_path.size() +
                       v.contacts.size() * (sizeof(uint64_t) * 2 + sizeof(float) + sizeof(uint32_t)) +
                       sizeof(uint32_t) * 2 * v.contacts.size(); /* name lengths */

    for (std::size_t i = 0; i < lhs_names.size(); ++i) {
      size += lhs_names[i].size() + rhs_names[i].size();
    }

    std::string buffer;
    buffer.reserve(size);

    const auto append_pod = [&buffer](auto value) {
      using T = decltype(value);
      buffer.append(reinterpret_cast<const char*>(&value), sizeof(T));
    };

    buffer.push_back(static_cast<char>(VERSION));
    buffer.push_back(static_cast<char>(v.success ? 1 : 0));
    buffer.push_back(static_cast<char>(static_cast<uint8_t>(v.provider)));

    const uint32_t contact_type = static_cast<uint32_t>(v.contact_type);
    append_pod(contact_type);

    const uint32_t path_len = checked_len(v.file_path.size());
    append_pod(path_len);
    buffer.append(v.file_path.data(), v.file_path.size());

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

    return buffer;
  }

  static ContactsRes deserialize(const char* data, std::size_t size) {
    constexpr uint8_t VERSION = 1;

    auto require = [size](std::size_t offset, std::size_t need) {
      if (offset + need > size) {
        throw std::runtime_error("ContactsRes binary deserialize: truncated payload");
      }
    };

    auto read_u8 = [&](std::size_t &offset) -> uint8_t {
      require(offset, sizeof(uint8_t));
      uint8_t value = static_cast<uint8_t>(data[offset]);
      offset += sizeof(uint8_t);
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

    result.success = read_u8(offset) != 0;
    const uint8_t provider = read_u8(offset);
    result.provider = static_cast<analysis::contacts::ContactProvider>(provider);

    const uint32_t contact_type = read_u32(offset);
    result.contact_type.category = static_cast<Category>(contact_type & 0xFFFFu);
    result.contact_type.flavor   = static_cast<Flavor>((contact_type >> 16) & 0xFFFFu);

    const uint32_t path_len = read_u32(offset);
    require(offset, path_len);
    result.file_path.assign(data + offset, path_len);
    offset += path_len;

    const uint32_t num_contacts = read_u32(offset);
    result.num_contacts = num_contacts;

    std::vector<Contact> contacts;
    contacts.reserve(result.num_contacts);

    for (std::size_t i = 0; i < result.num_contacts; ++i) {
      EntityID lhs;
      EntityID rhs;
      lhs.raw = read_u64(offset);
      rhs.raw = read_u64(offset);
      const float distance    = read_f32(offset);
      const uint32_t type_raw = read_u32(offset);

      InteractionType type;
      type.category = static_cast<Category>(type_raw & 0xFFFFu);
      type.flavor = static_cast<Flavor>((type_raw >> 16) & 0xFFFFu);

      contacts.emplace_back(lhs, rhs, distance, type);

      const uint32_t lhs_len = read_u32(offset);
      skip_bytes(offset, lhs_len);
      const uint32_t rhs_len = read_u32(offset);
      skip_bytes(offset, rhs_len);
    }

    result.contacts = ContactSet(std::move(contacts), true);
    result.topology = nullptr;

    return result;
  }

  static ContactsRes deserialize(const std::string &buf) { return deserialize(buf.data(), buf.size()); }
};

} // namespace serialization

#endif // LAHUTA_SERIALIZATION_CONTACTS_SERIALIZER_HPP
