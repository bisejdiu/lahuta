#ifndef LAHUTA_SERIALIZATION_CONTACTS_SERIALIZER_HPP
#define LAHUTA_SERIALIZATION_CONTACTS_SERIALIZER_HPP

#include "entities/formatter.hpp"
#include "entities/interaction_types.hpp"
#include "serialization/formats.hpp"
#include "serialization/json.hpp"
#include "serialization/serializer_impl.hpp"
#include "tasks/contacts_task.hpp"
#include <sstream>

// clang-format off
namespace serialization {
using namespace lahuta;

using ContactsRes = tasks::ContactsTask::result_type;

template<>
struct Serializer<fmt::json, ContactsRes> {
  using Record = ContactsRes;

  static std::string serialize(const ContactsRes& v) {
    JsonBuilder builder;

    std::string contact_type_str = (v.contact_type == InteractionType::None)
      ? "All"
      : interaction_type_to_string(v.contact_type);

    builder.key("file_path")   .value(v.file_path)
           .key("success")     .value(v.success)
           .key("provider")    .value(v.provider == tasks::ContactProvider::Arpeggio ? "arpeggio" : "molstar")
           .key("contact_type").value(contact_type_str)
           .key("num_contacts").value(v.num_contacts);

    if (!v.topology) {
      Logger::get_logger()->warn("ContactsRes serialization: topology is null, cannot serialize contacts.");
      return builder.str();
    }

    builder.key("contacts").begin_array();
    for (const auto& contact : v.contacts) {
      builder.begin_object()
             .key("lhs")     .value(ContactTableFormatter::format_entity_compact(*v.topology, contact.lhs))
             .key("rhs")     .value(ContactTableFormatter::format_entity_compact(*v.topology, contact.rhs))
             .key("distance").value(contact.distance)
             .key("type")    .value(interaction_type_to_string(contact.type))
             .end_object();
    }
    builder.end_array();
    return builder.str();
  }

  static ContactsRes deserialize(const std::string& s) {
    ContactsRes out;
    JsonReader r{s};

    out.file_path = r.get<std::string>("file_path");
    out.success   = r.get<bool>("success");

    std::string provider_str = r.get<std::string>("provider");
    out.provider = (provider_str == "arpeggio") ? tasks::ContactProvider::Arpeggio : tasks::ContactProvider::MolStar;

    std::string contact_type_str = r.get<std::string>("contact_type");
    out.contact_type = (contact_type_str == "All")
      ? InteractionType::None
      : get_interaction_type(contact_type_str);

    out.num_contacts = r.get<size_t>("num_contacts");
    out.topology = nullptr; // we cannot restore the topology from the serialized contact data

    return out;
  }
};

template<>
struct Serializer<fmt::text, ContactsRes> {
  using Record = ContactsRes;

  static std::string serialize(const ContactsRes& v) {
    std::ostringstream oss;

    std::string contact_type_str = (v.contact_type == InteractionType::None)
      ? "All"
      : interaction_type_to_string(v.contact_type);

    oss << (v.success ? "1" : "0") << " "
        << v.file_path << " "
        << (v.provider == tasks::ContactProvider::Arpeggio ? "arpeggio" : "molstar") << " "
        << contact_type_str << " "
        << v.num_contacts << "\n";

    if (!v.topology) return oss.str();

    for (const auto& contact : v.contacts) {
      oss << ContactTableFormatter::format_entity_compact(*v.topology, contact.lhs)
          << " -&- " << ContactTableFormatter::format_entity_compact(*v.topology, contact.rhs)
          << " " << std::fixed << std::setprecision(3) << contact.distance
          << " (" << interaction_type_to_string(contact.type) << ")";
      oss << "\n";
    }

    return oss.str();
  }

  static ContactsRes deserialize(const std::string& data) {
    // We should be able to re-create contacts from plain text but we need to be careful
    // with non-atom entities.
    throw std::runtime_error("Text format deserialization not implemented for ContactsRes");
  }
};

template<>
struct Serializer<fmt::binary, ContactsRes> {
  using Record = ContactsRes;

  static std::string serialize(const ContactsRes& v) {
    std::string buffer;

    size_t size = sizeof(bool) + sizeof(uint8_t) * 2 + sizeof(size_t) * 2 +
                  v.file_path.size() +
                  v.num_contacts * (sizeof(uint64_t) * 2 + sizeof(float) + sizeof(uint8_t));

    buffer.reserve(size);

    // Success flag
    buffer.append(reinterpret_cast<const char*>(&v.success), sizeof(v.success));

    uint8_t provider = static_cast<uint8_t>(v.provider);
    buffer.append(reinterpret_cast<const char*>(&provider), sizeof(provider));

    uint8_t contact_type = static_cast<uint8_t>(v.contact_type);
    buffer.append(reinterpret_cast<const char*>(&contact_type), sizeof(contact_type));

    size_t path_len = v.file_path.size();
    buffer.append(reinterpret_cast<const char*>(&path_len), sizeof(path_len));
    buffer.append(v.file_path);

    buffer.append(reinterpret_cast<const char*>(&v.num_contacts), sizeof(v.num_contacts));

    for (const auto& contact : v.contacts) {
      buffer.append(reinterpret_cast<const char*>(&contact.lhs.raw),  sizeof(contact.lhs.raw));
      buffer.append(reinterpret_cast<const char*>(&contact.rhs.raw),  sizeof(contact.rhs.raw));
      buffer.append(reinterpret_cast<const char*>(&contact.distance), sizeof(contact.distance));
      uint8_t type_raw = static_cast<uint8_t>(contact.type);
      buffer.append(reinterpret_cast<const char*>(&type_raw), sizeof(type_raw));
    }

    return buffer;
  }

  // FIX: many implicit assumptions, no explicit error handling
  static ContactsRes deserialize(const char* data, std::size_t size) {
    ContactsRes result;
    std::size_t offset = 0;

    std::memcpy(&result.success, data + offset, sizeof(result.success));
    offset += sizeof(result.success);

    uint8_t provider;
    std::memcpy(&provider, data + offset, sizeof(provider));
    result.provider = static_cast<tasks::ContactProvider>(provider);
    offset += sizeof(provider);

    // Contact type
    uint8_t contact_type;
    std::memcpy(&contact_type, data + offset, sizeof(contact_type));
    result.contact_type.category = static_cast<Category>(contact_type & 0x0F);
    result.contact_type.flavor   = static_cast<Flavor>((contact_type >> 4) & 0x0F);
    offset += sizeof(contact_type);

    // File path
    std::size_t path_len;
    std::memcpy(&path_len, data + offset, sizeof(path_len));
    offset += sizeof(path_len);
    result.file_path.assign(data + offset, path_len);
    offset += path_len;

    // Contact count
    std::memcpy(&result.num_contacts, data + offset, sizeof(result.num_contacts));
    offset += sizeof(result.num_contacts);

    // Contacts
    std::vector<Contact> contacts;
    contacts.reserve(result.num_contacts);

    for (size_t i = 0; i < result.num_contacts; ++i) {
      EntityID lhs, rhs;
      float distance;
      uint8_t type_raw;

      std::memcpy(&lhs.raw, data  + offset, sizeof(lhs.raw));
      offset += sizeof(lhs.raw);
      std::memcpy(&rhs.raw, data  + offset, sizeof(rhs.raw));
      offset += sizeof(rhs.raw);
      std::memcpy(&distance, data + offset, sizeof(distance));
      offset += sizeof(distance);
      std::memcpy(&type_raw, data + offset, sizeof(type_raw));
      offset += sizeof(type_raw);

      InteractionType type;
      type.category = static_cast<Category>(type_raw & 0x0F);
      type.flavor   = static_cast<Flavor>((type_raw >> 4) & 0x0F);

      contacts.emplace_back(lhs, rhs, distance, type);
    }

    result.contacts = ContactSet(std::move(contacts), true);
    result.topology = nullptr; // Topology not restored from binary format

    return result;
  }

  static ContactsRes deserialize(const std::string& buf) {
    return deserialize(buf.data(), buf.size());
  }
};

} // namespace serialization

#endif // LAHUTA_SERIALIZATION_CONTACTS_SERIALIZER_HPP
