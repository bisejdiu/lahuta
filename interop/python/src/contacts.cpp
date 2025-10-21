#include <algorithm>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

#include "analysis/contacts/provider.hpp"
#include "compute/topology_snapshot.hpp"
#include "contacts/arpeggio/provider.hpp"
#include "contacts/engine.hpp"
#include "contacts/getcontacts/provider.hpp"
#include "contacts/molstar/provider.hpp"
#include "entities/contact.hpp"
#include "entities/entity_id.hpp"
#include "entities/interaction_types.hpp"
#include "interactions.hpp"
#include "numpy_utils.hpp"
#include "topology.hpp"

namespace py = pybind11;

// clang-format off
namespace lahuta::bindings {

namespace {

std::optional<InteractionTypeSet> parse_optional_python_interactions(py::handle obj) {
  if (!obj || obj.is_none()) return std::nullopt;
  InteractionTypeSet set;
  add_python_interactions(set, obj);
  if (set.empty()) {
    throw py::value_error("Interaction type selection cannot be empty");
  }
  return set;
}

constexpr uint8_t ContactsBinaryVersion = 1;

struct ContactsRecordColumns {
  bool success = false;
  std::string file_path;
  analysis::contacts::ContactProvider provider = analysis::contacts::ContactProvider::MolStar;
  InteractionTypeSet contact_types = InteractionTypeSet::all();
  std::uint64_t frame_index = 0;
  std::vector<std::uint64_t> lhs_ids;
  std::vector<std::uint64_t> rhs_ids;
  std::vector<float> distances;
  std::vector<std::uint32_t> type_codes;
  std::vector<std::string> lhs_names;
  std::vector<std::string> rhs_names;
};

inline void ensure_available(std::size_t offset, std::size_t need, std::size_t size) {
  if (offset + need > size) {
    throw std::runtime_error("decode_contacts_binary: truncated payload");
  }
}

ContactsRecordColumns decode_contacts_binary_payload(const char *data, std::size_t size) {
  ContactsRecordColumns decoded;
  std::size_t offset = 0;

  auto read_u8 = [&](void) -> uint8_t {
    ensure_available(offset, sizeof(uint8_t), size);
    uint8_t value = static_cast<uint8_t>(data[offset]);
    offset += sizeof(uint8_t);
    return value;
  };

  auto read_u32 = [&](void) -> uint32_t {
    ensure_available(offset, sizeof(uint32_t), size);
    uint32_t value;
    std::memcpy(&value, data + offset, sizeof(uint32_t));
    offset += sizeof(uint32_t);
    return value;
  };

  auto read_u64 = [&](void) -> uint64_t {
    ensure_available(offset, sizeof(uint64_t), size);
    uint64_t value;
    std::memcpy(&value, data + offset, sizeof(uint64_t));
    offset += sizeof(uint64_t);
    return value;
  };

  auto read_f32 = [&](void) -> float {
    ensure_available(offset, sizeof(float), size);
    float value;
    std::memcpy(&value, data + offset, sizeof(float));
    offset += sizeof(float);
    return value;
  };

  auto read_string = [&](std::string &out) {
    const uint32_t len = read_u32();
    ensure_available(offset, len, size);
    out.assign(data + offset, len);
    offset += len;
  };

  const uint8_t version = read_u8();
  if (version != ContactsBinaryVersion) {
    throw std::runtime_error("decode_contacts_binary: unsupported version");
  }

  decoded.success = read_u8() != 0;
  decoded.provider = static_cast<analysis::contacts::ContactProvider>(read_u8());

  const uint8_t filter_mode = read_u8();
  const uint32_t filter_count = read_u32();
  InteractionTypeSet filters;
  if (filter_mode == 0) {
    filters = InteractionTypeSet::all();
    for (uint32_t i = 0; i < filter_count; ++i) {
      (void)read_u32();
    }
  } else {
    for (uint32_t i = 0; i < filter_count; ++i) {
      const uint32_t type_raw = read_u32();
      InteractionType type{};
      type.category = static_cast<Category>(type_raw & 0xFFFFu);
      type.flavor   = static_cast<Flavor>((type_raw >> 16) & 0xFFFFu);
      filters.add(type);
    }
  }
  decoded.contact_types = filters;

  std::string file_path;
  read_string(file_path);
  decoded.file_path = std::move(file_path);

  decoded.frame_index = read_u64();

  const uint32_t num_contacts = read_u32();
  decoded.lhs_ids   .reserve(num_contacts);
  decoded.rhs_ids   .reserve(num_contacts);
  decoded.distances .reserve(num_contacts);
  decoded.type_codes.reserve(num_contacts);
  decoded.lhs_names .reserve(num_contacts);
  decoded.rhs_names .reserve(num_contacts);

  for (uint32_t i = 0; i < num_contacts; ++i) {
    const uint64_t lhs_raw = read_u64();
    const uint64_t rhs_raw = read_u64();
    const float dist = read_f32();
    const uint32_t type_code = read_u32();

    decoded.lhs_ids   .push_back(lhs_raw);
    decoded.rhs_ids   .push_back(rhs_raw);
    decoded.distances .push_back(dist);
    decoded.type_codes.push_back(type_code);

    std::string lhs_name;
    std::string rhs_name;
    read_string(lhs_name);
    read_string(rhs_name);
    decoded.lhs_names.emplace_back(std::move(lhs_name));
    decoded.rhs_names.emplace_back(std::move(rhs_name));
  }

  return decoded;
}

py::dict make_contacts_numpy(ContactsRecordColumns decoded) {
  using namespace lahuta::numpy;

  std::vector<std::string> type_strings;
  type_strings.reserve(decoded.type_codes.size());
  for (auto code : decoded.type_codes) {
    InteractionType type{};
    type.category = static_cast<Category>(code & 0xFFFFu);
    type.flavor   = static_cast<Flavor>((code >> 16) & 0xFFFFu);
    type_strings.emplace_back(interaction_type_to_string(type));
  }

  const std::size_t num_contacts = decoded.lhs_ids.size();

  py::dict contacts;
  contacts["lhs"]      = string_array_1d(decoded.lhs_names);
  contacts["rhs"]      = string_array_1d(decoded.rhs_names);
  contacts["distance"] = move_to_numpy_owning<float>(std::move(decoded.distances));
  contacts["type"]     = string_array_1d(type_strings);

  py::dict out;
  out["success"]      = decoded.success;
  out["file_path"]    = decoded.file_path;
  out["provider"]     = std::string(contact_provider_name(decoded.provider));
  out["contact_type"] = interaction_type_set_to_string(decoded.contact_types, "|");
  out["num_contacts"] = static_cast<uint32_t>(num_contacts);
  out["frame_index"]  = py::int_(decoded.frame_index);
  out["contacts"]     = std::move(contacts);
  return out;
}

py::dict decode_contacts_to_dict_direct(ContactsRecordColumns decoded) {
  const std::size_t num_contacts = decoded.lhs_ids.size();

  // Build contacts list - one dict per contact (matches JSON structure)
  py::list contacts_list(num_contacts);

  for (std::size_t i = 0; i < num_contacts; ++i) {
    // Decode interaction type
    InteractionType type{};
    type.category = static_cast<Category>(decoded.type_codes[i] & 0xFFFFu);
    type.flavor   = static_cast<Flavor>((decoded.type_codes[i] >> 16) & 0xFFFFu);

    // Build contact dict
    py::dict contact;
    contact["lhs"]      = py::str(decoded.lhs_names[i]);
    contact["rhs"]      = py::str(decoded.rhs_names[i]);
    contact["distance"] = py::float_(decoded.distances[i]);
    contact["type"]     = py::str(interaction_type_to_string(type));

    contacts_list[i] = std::move(contact);
  }

  // Build result dict with metadata
  py::dict result;
  result["success"]      = decoded.success;
  result["file_path"]    = decoded.file_path;
  result["provider"]     = std::string(contact_provider_name(decoded.provider));
  result["contact_type"] = interaction_type_set_to_string(decoded.contact_types, "|");
  result["num_contacts"] = static_cast<uint32_t>(num_contacts);
  result["frame_index"]  = py::int_(decoded.frame_index);
  result["contacts"]     = std::move(contacts_list);

  return result;
}

py::dict decode_contacts_to_dict_columnar(ContactsRecordColumns decoded) {
  const std::size_t num_contacts = decoded.lhs_ids.size();

  // Pre-allocate columnar lists
  py::list lhs_list(num_contacts);
  py::list rhs_list(num_contacts);
  py::list dist_list(num_contacts);
  py::list type_list(num_contacts);

  for (std::size_t i = 0; i < num_contacts; ++i) {
    InteractionType type{};
    type.category = static_cast<Category>(decoded.type_codes[i] & 0xFFFFu);
    type.flavor   = static_cast<Flavor>((decoded.type_codes[i] >> 16) & 0xFFFFu);

    lhs_list[i]  = py::str(decoded.lhs_names[i]);
    rhs_list[i]  = py::str(decoded.rhs_names[i]);
    dist_list[i] = py::float_(decoded.distances[i]);
    type_list[i] = py::str(interaction_type_to_string(type));
  }

  py::dict contacts;
  contacts["lhs"]      = std::move(lhs_list);
  contacts["rhs"]      = std::move(rhs_list);
  contacts["distance"] = std::move(dist_list);
  contacts["type"]     = std::move(type_list);

  py::dict result;
  result["success"]      = decoded.success;
  result["file_path"]    = decoded.file_path;
  result["provider"]     = std::string(contact_provider_name(decoded.provider));
  result["contact_type"] = interaction_type_set_to_string(decoded.contact_types, "|");
  result["num_contacts"] = static_cast<uint32_t>(num_contacts);
  result["frame_index"]  = py::int_(decoded.frame_index);
  result["contacts"]     = std::move(contacts);

  return result;
}

py::dict decode_contacts_binary(py::bytes payload) {
  char* data = nullptr;
  Py_ssize_t len = 0;
  if (PyBytes_AsStringAndSize(payload.ptr(), &data, &len) != 0) {
    throw py::value_error("decode_contacts_binary: expected bytes");
  }
  if (len <= 0) {
    py::dict empty;
    empty["success"]      = false;
    empty["file_path"]    = "";
    empty["provider"]     = "unknown";
    empty["contact_type"] = "";
    empty["num_contacts"] = 0;
    empty["frame_index"]  = 0;
    empty["contacts"]     = py::dict();
    return empty;
  }

  ContactsRecordColumns decoded = decode_contacts_binary_payload(data, static_cast<std::size_t>(len));
  return make_contacts_numpy(std::move(decoded));
}

py::dict decode_contacts_binary_direct(py::bytes payload) {
  char* data = nullptr;
  Py_ssize_t len = 0;
  if (PyBytes_AsStringAndSize(payload.ptr(), &data, &len) != 0) {
    throw py::value_error("decode_contacts_binary_direct: expected bytes");
  }

  if (len <= 0) {
    py::dict empty;
    empty["success"]      = false;
    empty["file_path"]    = "";
    empty["provider"]     = "unknown";
    empty["contact_type"] = "";
    empty["num_contacts"] = 0;
    empty["frame_index"]  = 0;
    empty["contacts"]     = py::list();
    return empty;
  }

  ContactsRecordColumns decoded = decode_contacts_binary_payload(data, static_cast<std::size_t>(len));
  return decode_contacts_to_dict_direct(std::move(decoded));
}

py::dict decode_contacts_binary_columnar(py::bytes payload) {
  char* data = nullptr;
  Py_ssize_t len = 0;
  if (PyBytes_AsStringAndSize(payload.ptr(), &data, &len) != 0) {
    throw py::value_error("decode_contacts_binary_columnar: expected bytes");
  }

  if (len <= 0) {
    py::dict empty;
    empty["success"]      = false;
    empty["file_path"]    = "";
    empty["provider"]     = "unknown";
    empty["contact_type"] = "";
    empty["num_contacts"] = 0;
    empty["frame_index"]  = 0;
    py::dict contacts;
    contacts["lhs"]      = py::list();
    contacts["rhs"]      = py::list();
    contacts["distance"] = py::list();
    contacts["type"]     = py::list();
    empty["contacts"]    = std::move(contacts);
    return empty;
  }

  ContactsRecordColumns decoded = decode_contacts_binary_payload(
    data, static_cast<std::size_t>(len)
  );

  return decode_contacts_to_dict_columnar(std::move(decoded));
}

py::list decode_contacts_batch_parallel(py::list payloads, bool columnar) {
  const std::size_t num_payloads = payloads.size();

  if (num_payloads == 0) {
    return py::list();
  }

  py::list results(num_payloads);

  struct PayloadInfo {
    const char* data;
    std::size_t size;
  };

  std::vector<PayloadInfo> payload_infos(num_payloads);

  for (std::size_t i = 0; i < num_payloads; ++i) {
    py::bytes payload = py::cast<py::bytes>(payloads[i]);
    char* data = nullptr;
    Py_ssize_t len = 0;

    if (PyBytes_AsStringAndSize(payload.ptr(), &data, &len) != 0) {
      throw py::value_error("decode_contacts_batch_parallel: expected bytes in list");
    }

    payload_infos[i] = {data, static_cast<std::size_t>(len)};
  }

  const unsigned int num_threads = std::min(std::max(1u, std::thread::hardware_concurrency()), static_cast<unsigned int>(num_payloads));

  py::gil_scoped_release release;

  std::vector<ContactsRecordColumns> decoded_batch(num_payloads);
  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  const std::size_t chunk_size = (num_payloads + num_threads - 1) / num_threads;

  for (unsigned int t = 0; t < num_threads; ++t) {
    const std::size_t start_idx = t * chunk_size;
    const std::size_t end_idx = std::min(start_idx + chunk_size, num_payloads);

    if (start_idx >= num_payloads) break;

    threads.emplace_back([&payload_infos, &decoded_batch, start_idx, end_idx]() {
      for (std::size_t i = start_idx; i < end_idx; ++i) {
        decoded_batch[i] = decode_contacts_binary_payload(
          payload_infos[i].data,
          payload_infos[i].size
        );
      }
    });
  }

  for (auto& thread : threads) {
    thread.join();
  }

  py::gil_scoped_acquire acquire;
  for (std::size_t i = 0; i < num_payloads; ++i) {
    if (columnar) {
      results[i] = decode_contacts_to_dict_columnar(std::move(decoded_batch[i]));
    } else {
      results[i] = decode_contacts_to_dict_direct(std::move(decoded_batch[i]));
    }
  }

  return results;
}

} // namespace

void bind_contacts(py::module_ &m) {

  m.def("decode_contacts_binary",          &decode_contacts_binary,          py::arg("payload")); // numpy output
  m.def("decode_contacts_binary_direct",   &decode_contacts_binary_direct,   py::arg("payload")); // list of dicts output, columnar=False
  m.def("decode_contacts_binary_columnar", &decode_contacts_binary_columnar, py::arg("payload")); // dict of lists output, columnar=True
  m.def("decode_contacts_batch_parallel",  &decode_contacts_batch_parallel,  py::arg("payloads"), py::arg("columnar") = true); // parallel batch decode

  py::enum_<analysis::contacts::ContactProvider>(m, "ContactProvider")
    .value("MolStar",     analysis::contacts::ContactProvider::MolStar)
    .value("Arpeggio",    analysis::contacts::ContactProvider::Arpeggio)
    .value("GetContacts", analysis::contacts::ContactProvider::GetContacts);

  py::enum_<Category>(m, "Category")
    .value("None_",                 Category::None)
    .value("Generic",               Category::Generic)
    .value("Hydrophobic",           Category::Hydrophobic)
    .value("Halogen",               Category::Halogen)
    .value("HydrogenBond",          Category::HydrogenBond)
    .value("WeakHydrogenBond",      Category::WeakHydrogenBond)
    .value("PolarHydrogenBond",     Category::PolarHydrogenBond)
    .value("WeakPolarHydrogenBond", Category::WeakPolarHydrogenBond)
    .value("Aromatic",              Category::Aromatic)
    .value("Ionic",                 Category::Ionic)
    .value("MetalCoordination",     Category::MetalCoordination)
    .value("CationPi",              Category::CationPi)
    .value("PiStacking",            Category::PiStacking)
    .value("Carbonyl",              Category::Carbonyl)
    .value("VanDerWaals",           Category::VanDerWaals)
    .value("DonorPi",               Category::DonorPi)
    .value("SulphurPi",             Category::SulphurPi)
    .value("CarbonPi",              Category::CarbonPi);

  py::enum_<Flavor>(m, "Flavor")
    .value("Default",  Flavor::Default)
    .value("Parallel", Flavor::Parallel)
    .value("TShape",   Flavor::TShape);

  py::class_<InteractionType>(m, "InteractionType")
    .def(py::init<Category, Flavor>(),
         py::arg_v("category", Category::None, "Category.None_"),
         py::arg_v("flavor",   Flavor::Default, "Flavor.Default"))
    .def_readonly_static("All",                   &InteractionType::All)
    .def_readonly_static("None_",                 &InteractionType::None)
    .def_readonly_static("Generic",               &InteractionType::Generic)
    .def_readonly_static("Hydrophobic",           &InteractionType::Hydrophobic)
    .def_readonly_static("Halogen",               &InteractionType::Halogen)
    .def_readonly_static("Ionic",                 &InteractionType::Ionic)
    .def_readonly_static("CationPi",              &InteractionType::CationPi)
    .def_readonly_static("HydrogenBond",          &InteractionType::HydrogenBond)
    .def_readonly_static("WeakHydrogenBond",      &InteractionType::WeakHydrogenBond)
    .def_readonly_static("PolarHydrogenBond",     &InteractionType::PolarHydrogenBond)
    .def_readonly_static("WeakPolarHydrogenBond", &InteractionType::WeakPolarHydrogenBond)
    .def_readonly_static("MetalCoordination",     &InteractionType::MetalCoordination)
    .def_readonly_static("Aromatic",              &InteractionType::Aromatic)
    .def_readonly_static("PiStacking",            &InteractionType::PiStacking)
    .def_readonly_static("PiStackingP",           &InteractionType::PiStackingP)
    .def_readonly_static("PiStackingT",           &InteractionType::PiStackingT)
    .def_readonly_static("Carbonyl",              &InteractionType::Carbonyl)
    .def_readonly_static("VanDerWaals",           &InteractionType::VanDerWaals)
    .def_readonly_static("DonorPi",               &InteractionType::DonorPi)
    .def_readonly_static("SulphurPi",             &InteractionType::SulphurPi)
    .def_readonly_static("CarbonPi",              &InteractionType::CarbonPi)

    .def_property_readonly("category", [](const InteractionType &self){ return self.category; })
    .def_property_readonly("flavor",   [](const InteractionType &self){ return self.flavor; })

    .def("__str__",  [](const InteractionType& self) { return interaction_type_to_string(self); })
    .def("__repr__", [](const InteractionType& self) { return "InteractionType(" + interaction_type_to_string(self) + ")"; })
    .def(
      "__int__",
      [](const InteractionType& self) { return static_cast<unsigned int>(static_cast<std::uint32_t>(self)); },
      R"doc(
Return the 32-bit packed code for this InteractionType.

Layout: [ flavor:16 | category:16 ]. Both Category and Flavor are 16-bit enums.
)doc")

    .def("__hash__", [](const InteractionType& self) {
        return static_cast<std::size_t>(static_cast<std::uint32_t>(self));
      })
    .def("__eq__", [](const InteractionType& self, py::object other) -> py::object {
        if (py::isinstance<InteractionType>(other)) {
          auto o = other.cast<InteractionType>();
          return py::bool_(self == o);
        }
        return py::reinterpret_borrow<py::object>(Py_NotImplemented);
      })
    .def("__ne__", [](const InteractionType& self, py::object other) -> py::object {
        if (py::isinstance<InteractionType>(other)) {
          auto o = other.cast<InteractionType>();
          return py::bool_(self != o);
        }
        return py::reinterpret_borrow<py::object>(Py_NotImplemented);
      })
    .def("__or__", [](const InteractionType& lhs, const InteractionType& rhs) {
        return InteractionTypeSet{lhs, rhs};
      }, py::is_operator())
    .def("__or__", [](const InteractionType& lhs, const InteractionTypeSet& rhs) {
        InteractionTypeSet out = rhs;
        out |= lhs;
        return out;
      }, py::is_operator())
    .def("__ror__", [](const InteractionType& rhs, const InteractionType& lhs) {
        return InteractionTypeSet{lhs, rhs};
      }, py::is_operator())
    .def("__ror__", [](const InteractionType& rhs, const InteractionTypeSet& lhs) {
        InteractionTypeSet out = lhs;
        out |= rhs;
        return out;
      }, py::is_operator());

  py::class_<InteractionTypeSet>(m, "InteractionTypeSet")
    .def(py::init<>())
    .def(py::init<InteractionType>())
    .def(py::init<std::initializer_list<InteractionType>>())
    .def_static("all", []() { return InteractionTypeSet::all(); })
    .def("is_all", &InteractionTypeSet::is_all)
    .def("empty", &InteractionTypeSet::empty)
    .def("members", [](const InteractionTypeSet& self) {
        auto members = self.members();
        py::list out;
        for (auto type : members) out.append(type);
        return out;
      })
    .def("__len__", [](const InteractionTypeSet& self) { return static_cast<std::size_t>(self.count()); })
    .def("__contains__", [](const InteractionTypeSet& self, const InteractionType& value) {
        return self.contains(value);
      })
    .def("__or__", [](const InteractionTypeSet& lhs, const InteractionTypeSet& rhs) {
        auto out = lhs;
        out |= rhs;
        return out;
      }, py::is_operator())
    .def("__or__", [](const InteractionTypeSet& lhs, const InteractionType& rhs) {
        auto out = lhs;
        out |= rhs;
        return out;
      }, py::is_operator())
    .def("__ror__", [](const InteractionTypeSet& rhs, const InteractionTypeSet& lhs) {
        auto out = lhs;
        out |= rhs;
        return out;
      }, py::is_operator())
    .def("__ror__", [](const InteractionTypeSet& rhs, const InteractionType& lhs) {
        auto out = rhs;
        out |= lhs;
        return out;
      }, py::is_operator())
    .def("__ior__", [](InteractionTypeSet& lhs, const InteractionTypeSet& rhs) -> InteractionTypeSet& {
        lhs |= rhs;
        return lhs;
      }, py::is_operator())
    .def("__ior__", [](InteractionTypeSet& lhs, const InteractionType& rhs) -> InteractionTypeSet& {
        lhs |= rhs;
        return lhs;
      }, py::is_operator())
    .def("__str__", [](const InteractionTypeSet& self) {
        return interaction_type_set_to_string(self, "|");
      })
    .def("__repr__", [](const InteractionTypeSet& self) {
        return std::string("InteractionTypeSet(") + interaction_type_set_to_string(self, "|") + ")";
      });

  // A contact between two entities with distance and interaction type
    py::class_<Contact>(m, "Contact")
      .def(py::init<EntityID, EntityID, float, InteractionType>(),
           py::arg("lhs"), py::arg("rhs"), py::arg("distance_sq"), py::arg("type"))
      .def_property(
        "lhs",
        [](const Contact &self) { return self.lhs; },
        [](Contact &self, const EntityID &value) { self.lhs = value; },
        "Left entity (EntityID)"
      )
      .def_property(
        "rhs",
        [](const Contact &self) { return self.rhs; },
        [](Contact &self, const EntityID &value) { self.rhs = value; },
        "Right entity (EntityID)"
      )
      .def_property("distance_sq",
        [](const Contact &self) { return self.distance; },
        [](Contact &self, float value) { self.distance = value; },
        "Distance squared between entities (A^2)" )
      .def_property("type", 
        [](const Contact &self) { return self.type; },
        [](Contact &self, const InteractionType &value) { self.type = value; },
        "Interaction type (InteractionType)")

      .def("__hash__", [](const Contact& c) {
          size_t h1 = std::hash<uint64_t>{}(c.lhs.raw);
          size_t h2 = std::hash<uint64_t>{}(c.rhs.raw);
          size_t h3 = std::hash<uint32_t>{}(static_cast<uint32_t>(c.type));
          return (h1 ^ (h2 << 1)) ^ (h3 << 2);
      })
      .def("__eq__", [](const Contact &self, py::object other) -> bool {
        if (!py::isinstance<Contact>(other)) return false;
        return self == other.cast<Contact>();
      }, py::is_operator())
      .def("__ne__", [](const Contact &self, py::object other) -> bool {
        if (!py::isinstance<Contact>(other)) return true;
        return self != other.cast<Contact>();
      }, py::is_operator())
      .def("__lt__", &Contact::operator<,  py::is_operator())
      .def("__repr__", [](const Contact &self) {
        return py::str("Contact(lhs={}, rhs={}, type={}, distance_sq={})")
          .format(self.lhs.to_string(), self.rhs.to_string(), interaction_type_to_string(self.type), self.distance);
      })
      .def("__str__", [](const Contact &self) {
        return py::str("Contact(lhs={}, rhs={}, type={}, distance_sq={})")
          .format(self.lhs.to_string(), self.rhs.to_string(), interaction_type_to_string(self.type), self.distance);
      })
      .def("to_dict", [](const Contact &self) {
        py::dict d;
        d[py::str("lhs_kind")]    = self.lhs.kind();
        d[py::str("lhs_index")]   = self.lhs.index();
        d[py::str("rhs_kind")]    = self.rhs.kind();
        d[py::str("rhs_index")]   = self.rhs.index();
        d[py::str("distance_sq")] = self.distance;
        d[py::str("type")]        = self.type;
        d[py::str("category")]    = self.type.category;
        d[py::str("flavor")]      = self.type.flavor;
        return d;
      }, "Return a dict representation suitable for testing/serialization")
      ;

  // Ordered container for Contact obj with set operations and no duplicates
  py::class_<ContactSet>(m, "ContactSet")
    .def(py::init<>())
    .def("data",         [](ContactSet &self) { return self.data(); })
    .def("insert",       py::overload_cast<const ContactSet &>(&ContactSet::insert), py::arg("contacts"))
    .def("insert",       py::overload_cast<const Contact    &>(&ContactSet::insert), py::arg("contact" ))

    .def("set_union",         &ContactSet::set_union, py::arg("other"))
    .def("set_intersection",  &ContactSet::set_intersection, py::arg("other"))
    .def("set_difference",    &ContactSet::set_difference, py::arg("other"))
    .def("set_symmetric_difference", &ContactSet::set_symmetric_difference, py::arg("other"))

    .def("size",         &ContactSet::size)
    .def("empty",        &ContactSet::empty)
    .def("make_generic", &ContactSet::make_generic)

    .def("__len__",      &ContactSet::size)
    .def("__and__",  [](const ContactSet &a, const ContactSet &b) { return a &  b; }, py::is_operator())
    .def("__or__",   [](const ContactSet &a, const ContactSet &b) { return a |  b; }, py::is_operator())
    .def("__sub__",  [](const ContactSet &a, const ContactSet &b) { return a -  b; }, py::is_operator())
    .def("__xor__",  [](const ContactSet &a, const ContactSet &b) { return a ^  b; }, py::is_operator())
    .def("__iand__", [](      ContactSet &a, const ContactSet &b) { return a &= b; }, py::is_operator())
    .def("__ior__",  [](      ContactSet &a, const ContactSet &b) { return a |= b; }, py::is_operator())
    .def("__isub__", [](      ContactSet &a, const ContactSet &b) { return a -= b; }, py::is_operator())
    .def("__ixor__", [](      ContactSet &a, const ContactSet &b) { return a ^= b; }, py::is_operator())

    .def("__getitem__", [](const ContactSet &self, int index) -> Contact {
      if (index < 0 || static_cast<size_t>(index) >= self.size()) throw py::index_error("Index out of range");
      return self.data()[static_cast<size_t>(index)];
    })
    .def("__iter__", [](const ContactSet &self) { return py::make_iterator(self.begin(), self.end()); }, py::keep_alive<0, 1>());

  using MsEngine = InteractionEngine<MolStarContactProvider>;
  using AgEngine = InteractionEngine<ArpeggioContactProvider>;
  using GcEngine = InteractionEngine<GetContactsProvider>;

  py::class_<MsEngine>(m, "MolStarContactsEngine")
    .def(py::init<>())
    .def("compute", [](const MsEngine& eng, const Topology& topo) {
        auto ts = compute::snapshot_of(topo, topo.conformer());
        return eng.compute(ts);
      }, py::arg("topology"))
    .def("compute", [](const MsEngine& eng, const Topology& topo, py::object only) {
        auto ts = compute::snapshot_of(topo, topo.conformer());
        auto filter = parse_optional_python_interactions(only);
        if (!filter || filter->is_all()) return eng.compute(ts);
        return eng.compute(ts, std::move(filter));
      }, py::arg("topology"), py::arg("only") = py::none());

  py::class_<AgEngine>(m, "ArpeggioContactsEngine")
    .def(py::init<>())
    .def("compute", [](const AgEngine& eng, const Topology& topo) {
        auto ts = compute::snapshot_of(topo, topo.conformer());
        return eng.compute(ts);
      }, py::arg("topology"))
    .def("compute", [](const AgEngine& eng, const Topology& topo, py::object only) {
        auto ts = compute::snapshot_of(topo, topo.conformer());
        auto filter = parse_optional_python_interactions(only);
        if (!filter || filter->is_all()) return eng.compute(ts);
        return eng.compute(ts, std::move(filter));
      }, py::arg("topology"), py::arg("only") = py::none());

  py::class_<GcEngine>(m, "GetContactsEngine")
    .def(py::init<>())
    .def("compute", [](const GcEngine& eng, const Topology& topo) {
        auto ts = compute::snapshot_of(topo, topo.conformer());
        return eng.compute(ts);
      }, py::arg("topology"))
    .def("compute", [](const GcEngine& eng, const Topology& topo, py::object only) {
        auto ts = compute::snapshot_of(topo, topo.conformer());
        auto filter = parse_optional_python_interactions(only);
        if (!filter || filter->is_all()) return eng.compute(ts);
        return eng.compute(ts, std::move(filter));
      }, py::arg("topology"), py::arg("only") = py::none());
}
} // namespace lahuta::bindings
