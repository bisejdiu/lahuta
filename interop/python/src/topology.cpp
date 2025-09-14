#include <memory>

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/typing.h>

#include "residues.hpp"
#include "topology.hpp"
#include "topology_flags.hpp"

namespace py = pybind11;

// clang-format off
namespace lahuta::bindings {

void bind_topology(py::module &m) {

  py::enum_<ContactComputerType>(m, "ContactComputerType", "Atom typing backends used when classifying atoms for contacts.")
    .value("None_",    ContactComputerType::None,     "Disable atom typing (no classification)")
    .value("Arpeggio", ContactComputerType::Arpeggio, "Use Arpeggio-style atom typing")
    .value("Molstar",  ContactComputerType::Molstar,  "Use Mol* atom typing");

  py::class_<TopologyBuildingOptions>(m, "TopologyBuildingOptions", "Options controlling topology construction.")
    .def(py::init<>(), "Create default options")
    .def_readwrite("cutoff",    &TopologyBuildingOptions::cutoff,    "Neighbor-search cutoff used in bond perception and neighbors (Å)")
    .def_readwrite("auto_heal", &TopologyBuildingOptions::auto_heal, "Resolve required dependencies automatically during execution")
    .def_readwrite("atom_typing_method",        &TopologyBuildingOptions::atom_typing_method,        "Backend used for atom typing")
    .def_readwrite("compute_nonstandard_bonds", &TopologyBuildingOptions::compute_nonstandard_bonds, "Include non-standard/metal coordination where applicable");

  py::enum_<TopologyComputation>(m, "TopologyComputers", "Bitmask flags representing topology stages. Combine with | and &.")
    .value("None_",            TopologyComputation::None,             "No stage selected")
    .value("Neighbors",        TopologyComputation::Neighbors,        "Compute neighbor lists (used by bond perception)")
    .value("Bonds",            TopologyComputation::Bonds,            "Perceive covalent bonds from neighbors/chemistry")
    .value("NonStandardBonds", TopologyComputation::NonStandardBonds, "Include non-standard bonds/coordination where supported")
    .value("Residues",         TopologyComputation::Residues,         "Assemble residue/group membership")
    .value("Rings",            TopologyComputation::Rings,            "Detect cycles and annotate aromatic rings")
    .value("AtomTyping",       TopologyComputation::AtomTyping,       "Assign atom types (per ContactComputerType)")
    .value("Basic",            TopologyComputation::Basic,            "Neighbors + Bonds")
    .value("Standard",         TopologyComputation::Standard,         "Basic + Residues")
    .value("Extended",         TopologyComputation::Extended,         "Standard + Rings")
    .value("Complete",         TopologyComputation::Complete,         "Extended + AtomTyping")
    .value("All",              TopologyComputation::All,              "Alias of Complete")
    .def("__or__",  [](TopologyComputation a, TopologyComputation b) { return a | b; }, "Bitwise OR: combine flags")
    .def("__and__", [](TopologyComputation a, TopologyComputation b) { return a & b; }, "Bitwise AND: intersect flags");

  py::class_<Residue>(m, "Residue", "Residue view with identity and atom membership. Comparable and hashable.")
    .def(py::init<>(), "Create empty residue")
    .def(py::init<const std::string &, int, const std::string &, const std::string &>(),
         py::arg("chain_id"), py::arg("number"), py::arg("name"), py::arg("alt_loc"),
         "Create residue with specified properties")
    .def_readwrite("idx",      &Residue::idx,      "0-based index of this residue in the parent molecule")
    .def_readwrite("chain_id", &Residue::chain_id, "Chain identifier")
    .def_readwrite("number",   &Residue::number,   "Residue number")
    .def_readwrite("name",     &Residue::name,     "Residue name")
    .def_readwrite("alt_loc",  &Residue::alt_loc,  "Alternative location identifier")
    .def_readwrite("atoms",    &Residue::atoms,    "Atom indices in this residue")

    .def("__hash__", [](const Residue& self) {
        return py::hash(py::make_tuple(self.chain_id, self.number, self.name, self.alt_loc));
    })
    .def("__eq__", [](const Residue& self, py::object other) -> py::object {
        if (py::isinstance<Residue>(other)) {
            const auto& o = other.cast<const Residue&>();
            return py::bool_(self == o);
        }
        return py::reinterpret_borrow<py::object>(Py_NotImplemented);
    })
    .def("__ne__", [](const Residue& self, py::object other) -> py::object {
        if (py::isinstance<Residue>(other)) {
            const auto& o = other.cast<const Residue&>();
            return py::bool_(self != o);
        }
        return py::reinterpret_borrow<py::object>(Py_NotImplemented);
    })
    .def("__repr__", [](const Residue& self) {
        return "Residue(chain_id='" + self.chain_id + "', number=" + std::to_string(self.number) +
               ", name='" + self.name + "', alt_loc='" + self.alt_loc + "')";
    })
    .def("__str__", [](const Residue& self) {
        return self.chain_id + ":" + std::to_string(self.number) + " " + self.name +
               (self.alt_loc.empty() ? "" : " (" + self.alt_loc + ")");
    });

  py::class_<Residues>(m, "Residues", "Random-access container over residues with functional helpers.")
    .def(py::init<const RDKit::RWMol &>(), py::arg("mol"), "Create residues from RDKit molecule")

    .def_property_readonly("residues", &Residues::get_residues, "Snapshot list of residues")
    .def("filter", [](const Residues &self, py::function func) {
        auto pred = [func](const Residue &r) {return func(r).cast<bool>();};
        return self.filter(pred);
     }, py::arg("func"),
      py::keep_alive<0, 1>(),
      "Return a new Residues filtered by predicate")

    .def("map", [](const Residues &self, py::function func) {
        std::vector<py::object> results;
        for (const auto &r : self.get_residues()) {
            results.push_back(func(r));
        }
        return results;
     }, py::arg("func"), "Apply function to each residue, returning a Python list")
    .def("get_atom_ids",      &Residues::get_atom_ids,      "Concatenated atom indices of all residues (order-preserving)")
    .def("get_residue_names", &Residues::get_residue_names, "Residue names aligned to container order")
    .def("__getitem__", [](Residues &self, size_t i) -> const Residue& {
        const auto &v = self.get_residues();  // NOTE: must be a reference-returning accessor
        if (i >= v.size()) throw py::index_error();
        return v[i];
      },
      py::return_value_policy::reference_internal,
      "Get residue by index (view kept alive by the container)")

    .def("__iter__", [](const Residues &self) {
        return py::make_iterator(self.begin(), self.end());
    }, py::keep_alive<0, 1>(), "Iterate over residues")
    .def("__len__", &Residues::size, "Get number of residues");

  py::class_<Topology, std::shared_ptr<Topology>>(m, "Topology", "Topology manager built over an RDKit molecule and conformer.")
    // Return Python lists of references tied to the parent Topology so that
    // individual records keep the owner alive even after the list is dropped.
    .def_property_readonly("atom_types", [](Topology &self) -> py::typing::List<AtomRec> {
        const auto &recs = self.records<AtomRec>();
        py::list out;
        for (const auto &rec : recs) {
          out.append(py::cast(&rec,
                              py::return_value_policy::reference_internal,
                              py::cast(&self, py::return_value_policy::reference)));
        }
        return out;
      },
      py::keep_alive<0, 1>(),
      "Per-atom typing records (size N_atoms)")
    .def_property_readonly("rings", [](Topology &self) -> py::typing::List<RingRec> {
        const auto &recs = self.records<RingRec>();
        py::list out;
        for (const auto &rec : recs) {
          out.append(py::cast(&rec,
                              py::return_value_policy::reference_internal,
                              py::cast(&self, py::return_value_policy::reference)));
        }
        return out;
      },
      py::keep_alive<0, 1>(),
      "Detected rings (may be empty)")
    .def_property_readonly("groups", [](Topology &self) -> py::typing::List<GroupRec> {
        const auto &recs = self.records<GroupRec>();
        py::list out;
        for (const auto &rec : recs) {
          out.append(py::cast(&rec,
                              py::return_value_policy::reference_internal,
                              py::cast(&self, py::return_value_policy::reference)));
        }
        return out;
      },
      py::keep_alive<0, 1>(),
      "Detected functional groups (may be empty)")

    .def_property_readonly("residues",
        [](Topology &top) -> const Residues& { return top.get_residues(); },
        py::return_value_policy::reference_internal,
        "Get residues container (kept alive by the Topology)")

    .def("atom_types_filter_by_fn", [](Topology &top, py::function func) {
        auto pred = [func](const AtomRec &r) {return func(r).cast<bool>();};
        auto recs = top.records<AtomRec>();
        auto filtered = std::vector<AtomRec>();
        filtered.reserve(recs.size());
        for (const auto &rec : recs) {
            if (pred(rec)) {
                filtered.push_back(rec);
            }
        }
        return filtered;
    }, py::arg("func"), "Return atom records for which the predicate returns True")

    .def("get_atom_ids", &Topology::get_atom_ids, "0-based atom indices present in topology")
    .def("build",     &Topology::build, py::arg("options"), "Build all enabled stages. Returns a boolean")

    .def("run_mask",  &Topology::run_mask, py::arg("mask"), "Run stages specified by a bitmask of TopologyComputers")

    .def("assign_molstar_typing",      &Topology::assign_molstar_typing,      "Assign per-atom types using Mol* rules; populates atom_types")
    .def("assign_arpeggio_atom_types", &Topology::assign_arpeggio_atom_types, "Assign per-atom types using Arpeggio rules; populates atom_types")

    .def("enable_computation",     &Topology::enable_computation,     py::arg("comp"), py::arg("enabled"), "Enable/disable a specific stage")
    .def("enable_only",            &Topology::enable_only,            py::arg("comps"),   "Enable only the provided bitmask; disables all others")
    .def("is_computation_enabled", &Topology::is_computation_enabled, py::arg("comp"),    "Whether a stage is enabled")
    .def("execute_computation",    &Topology::execute_computation,    py::arg("comp"),    "Run a single stage, resolving dependencies")
    .def("set_cutoff",             &Topology::set_cutoff,             py::arg("cutoff"),  "Neighbor cutoff used by bond perception (Å)")
    .def("set_atom_typing_method", &Topology::set_atom_typing_method, py::arg("method"),  "Set atom typing backend for future typing")
    .def("set_compute_nonstandard_bonds", &Topology::set_compute_nonstandard_bonds, py::arg("compute"), "Include metal/coordination bonds if True")

    .def("get_atom",  &Topology::atom,  py::arg("idx"), py::return_value_policy::reference_internal, "Atom record at atom index")
    .def("get_ring",  &Topology::ring,  py::arg("idx"), py::return_value_policy::reference_internal, "Ring record at index (0..n_rings-1)")
    .def("get_group", &Topology::group, py::arg("idx"), py::return_value_policy::reference_internal, "Group record at index (0..n_groups-1)")

    // Typed resolving helpers
    .def("resolve_atom", [](Topology &self, const EntityID &id) -> const AtomRec& {
        return self.resolve<Kind::Atom>(id);
      }, py::arg("id"), py::return_value_policy::reference_internal,
      "Resolve an Atom EntityID to AtomRec (reference, lifetime tied to Topology)")
    .def("resolve_ring", [](Topology &self, const EntityID &id) -> const RingRec& {
        return self.resolve<Kind::Ring>(id);
      }, py::arg("id"), py::return_value_policy::reference_internal,
      "Resolve a Ring EntityID to RingRec (reference, lifetime tied to Topology)")
    .def("resolve_group", [](Topology &self, const EntityID &id) -> const GroupRec& {
        return self.resolve<Kind::Group>(id);
      }, py::arg("id"), py::return_value_policy::reference_internal,
      "Resolve a Group EntityID to GroupRec (reference, lifetime tied to Topology)")

    .def("molecule",  [](const Topology &self) -> auto { return self.molecule();  }, py::return_value_policy::reference_internal, "Underlying RDKit molecule (borrowed reference)")
    .def("conformer", [](const Topology &self) -> auto { return self.conformer(); }, py::return_value_policy::reference_internal, "RDKit conformer used for positions (borrowed reference)");
}

} // namespace lahuta::bindings
