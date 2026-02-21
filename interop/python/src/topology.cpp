/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: [] {
 *   auto t1 = std::make_tuple("besian");
 *   auto t2 = std::make_tuple("sejdiu");
 *   auto t3 = std::make_tuple("@gmail.com");
 *   auto combined = std::tuple_cat(t1, t2, t3);
 *   return std::apply([](auto... args) {
 *     return (std::string{} + ... + std::string(args));
 *   }, combined);
 * }();
 *
 */

#include <memory>

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/typing.h>

#include "residues/residues.hpp"
#include "topology.hpp"
#include "topology_flags.hpp"
#include "typing/flags.hpp"

namespace py = pybind11;

// clang-format off
namespace lahuta::bindings {
namespace Flags = AtomTypeFlags;

// Dynamic views that compute geometry from the owning molecule
// and keep it alive independently of the Topology via shared_ptr.
struct RingView {
  std::shared_ptr<RDKit::RWMol> mol;
  std::vector<unsigned int> atom_idxs;
  bool aromatic_flag;

  std::size_t size() const { return atom_idxs.size(); }
  bool aromatic() const { return aromatic_flag; }

  RDGeom::Point3D center() const {
    const RDKit::Conformer &conf = mol->getConformer();
    RDGeom::Point3D c{0.0, 0.0, 0.0};
    if (atom_idxs.empty()) return c;
    for (auto idx : atom_idxs) c += conf.getAtomPos(idx);
    c /= static_cast<double>(atom_idxs.size());
    return c;
  }

  RDGeom::Point3D normal() const {
    RDGeom::Point3D nrm{0.0, 0.0, 0.0};
    if (atom_idxs.size() < 3) return nrm;
    const RDKit::Conformer &conf = mol->getConformer();
    auto p0 = conf.getAtomPos(atom_idxs[0]);
    auto p1 = conf.getAtomPos(atom_idxs[1]);
    auto p2 = conf.getAtomPos(atom_idxs[2]);
    auto v1 = p1 - p0;
    auto v2 = p2 - p0;
    nrm = v1.crossProduct(v2);
    const double len_sq = nrm.lengthSq();
    if (len_sq > 0.0) nrm /= std::sqrt(len_sq);
    return nrm;
  }

  py::list atoms() const {
    py::list out;
    for (auto idx : atom_idxs) {
      out.append(py::cast(mol->getAtomWithIdx(idx), py::return_value_policy::reference));
    }
    return out;
  }
};

struct GroupView {
  std::shared_ptr<RDKit::RWMol> mol;
  std::vector<unsigned int> atom_idxs;
  AtomType a_type_val;
  FeatureGroup type_val;

  AtomType a_type() const { return a_type_val; }
  FeatureGroup type() const { return type_val; }
  py::list atoms() const {
    py::list out;
    for (auto idx : atom_idxs) {
      out.append(py::cast(mol->getAtomWithIdx(idx), py::return_value_policy::reference));
    }
    return out;
  }
  RDGeom::Point3D center() const {
    const RDKit::Conformer &conf = mol->getConformer();
    RDGeom::Point3D c{0.0, 0.0, 0.0};
    if (atom_idxs.empty()) return c;
    for (auto idx : atom_idxs) c += conf.getAtomPos(idx);
    c /= static_cast<double>(atom_idxs.size());
    return c;
  }
};

void bind_topology(py::module &m) {
  // Bind dynamic views first
  py::class_<RingView>(m, "RingView", "Ring view with dynamic geometry accessors")
    .def_property_readonly("size",    &RingView::size,     "Number of atoms in the ring")
    .def_property_readonly("atoms",   &RingView::atoms,    "Atoms participating in the ring")
    .def_property_readonly("aromatic",&RingView::aromatic, "Whether the ring is aromatic")
    .def_property_readonly("center",  &RingView::center,   "Ring geometric center (A)")
    .def_property_readonly("normal",  &RingView::normal,   "Normal vector to ring plane")
    .def("__repr__", [](const RingView &self) {
      return py::str("RingView(atoms={}, aromatic={})").format(self.size(), self.aromatic());
    });

  py::class_<GroupView>(m, "GroupView", "Functional group view with dynamic geometry accessors")
    .def_property_readonly("a_type", &GroupView::a_type, "Atom type associated with this group")
    .def_property_readonly("type",   &GroupView::type,   "Functional group classification")
    .def_property_readonly("atoms",  &GroupView::atoms,  "Atoms participating in the group")
    .def_property_readonly("center", &GroupView::center, "Geometric center of the group (A)")
    .def("__repr__", [](const GroupView &self) {
      return py::str("GroupView(a_type={}, type={}, atoms={})").format(
        static_cast<uint32_t>(self.a_type()), static_cast<int>(self.type()), self.atoms().size());
    });

  py::enum_<AtomTypingMethod>(m, "AtomTypingMethod", "Atom typing backends used when classifying atoms for contacts.")
    .value("Arpeggio",    AtomTypingMethod::Arpeggio,    "Use Arpeggio-style atom typing")
    .value("MolStar",     AtomTypingMethod::Molstar,     "Use Mol* atom typing")
    .value("GetContacts", AtomTypingMethod::GetContacts, "Use getcontacts-style atom typing");

  py::class_<TopologyBuildingOptions>(m, "TopologyBuildingOptions", "Options controlling topology construction.")
    .def(py::init<>(), "Create default options")
    .def_readwrite("cutoff",    &TopologyBuildingOptions::cutoff,    "Neighbor-search cutoff used in bond perception and neighbors (A)")
    .def_readwrite("auto_heal", &TopologyBuildingOptions::auto_heal, "Resolve required dependencies automatically during execution")
    .def_readwrite("atom_typing_method",        &TopologyBuildingOptions::atom_typing_method,        "Backend used for atom typing")
    .def_readwrite("compute_nonstandard_bonds", &TopologyBuildingOptions::compute_nonstandard_bonds, "Include non-standard/metal coordination where applicable");

  py::enum_<TopologyComputation>(m, "TopologyComputers", "Bitmask flags representing topology stages. Combine with | and &.")
    .value("NoComp",           TopologyComputation::None,             "No stage selected")
    .value("Neighbors",        TopologyComputation::Neighbors,        "Compute neighbor lists (used by bond perception)")
    .value("Bonds",            TopologyComputation::Bonds,            "Perceive covalent bonds from neighbors/chemistry")
    .value("NonStandardBonds", TopologyComputation::NonStandardBonds, "Include non-standard bonds/coordination where supported")
    .value("Residues",         TopologyComputation::Residues,         "Assemble residue/group membership")
    .value("Rings",            TopologyComputation::Rings,            "Detect cycles and annotate aromatic rings")
    .value("AtomTyping",       TopologyComputation::AtomTyping,       "Assign atom types (per AtomTypingMethod)")
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
        auto pred = [func](const Residue &r) {
          py::object out = func(r);
          try {
            return out.cast<bool>();
          } catch (const py::cast_error &) {
            const std::string type_name = py::str(out.get_type().attr("__name__"));
            throw py::type_error(
              "Residues.filter predicate must return a bool (or bool-castable value), got '" +
              type_name + "'. Example: lambda r: r.name == 'ARG' or lambda r: True"
            );
          }
        };
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
    .def("residue_index_of_atom", &Residues::residue_index_of_atom, py::arg("atom_idx"), "Return residue index for an atom index")
    .def("residue_of_atom",       &Residues::residue_of_atom,       py::arg("atom_idx"),
      py::return_value_policy::reference_internal, "Return residue object for an atom index")
    .def_property_readonly("atom_to_residue_indices", &Residues::atom_to_residue_indices,
      "Atom-aligned residue index mapping (-1 for unmapped atoms)")
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
    .def_property_readonly("atom_records", [](Topology &self) -> py::typing::List<AtomRec> {
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
    .def_property_readonly("rings", [](Topology &self) -> py::typing::List<RingView> {
        const auto &recs = self.records<RingRec>();
        auto mol_ptr = self.molecule_ptr();
        py::list out;
        for (const auto &rec : recs) {
          std::vector<unsigned int> idxs;
          idxs.reserve(rec.atoms.size());
          for (const auto &a : rec.atoms) idxs.push_back(static_cast<unsigned int>(a.get().getIdx()));
          out.append(RingView{mol_ptr, std::move(idxs), rec.aromatic});
        }
        return out;
      }, py::keep_alive<0, 1>(),
      "Detected rings as dynamic views")
    .def_property_readonly("groups", [](Topology &self) -> py::typing::List<GroupView> {
        const auto &recs = self.records<GroupRec>();
        auto mol_ptr = self.molecule_ptr();
        py::list out;
        for (const auto &rec : recs) {
          std::vector<unsigned int> idxs;
          idxs.reserve(rec.atoms.size());
          for (const auto &a : rec.atoms) idxs.push_back(static_cast<unsigned int>(a.get().getIdx()));
          out.append(GroupView{mol_ptr, std::move(idxs), rec.a_type, rec.type});
        }
        return out;
      }, py::keep_alive<0, 1>(),
      "Detected functional groups as dynamic views")

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
    .def("atoms_with_type", [](Topology &self, AtomType type) -> py::typing::List<AtomRec> {
        const auto &recs = self.records<AtomRec>();
        py::list out;
        for (const auto &rec : recs) {
          if (Flags::has(rec.type, type)) {
            out.append(py::cast(&rec,
                                py::return_value_policy::reference_internal,
                                py::cast(&self, py::return_value_policy::reference)));
          }
        }
        return out;
      }, py::arg("type"),
      "Return atom records containing the specified type flag")

    .def("get_atom_ids", &Topology::get_atom_ids, "0-based atom indices present in topology")
    .def("build",     &Topology::build, py::arg("options"), "Build all enabled stages. Returns a boolean")

    .def("assign_typing", &Topology::assign_typing, py::arg("method"), "Assign per-atom types using specified method; populates atom_records")

    .def("has_computed",           &Topology::has_computed,           py::arg("comp"),    "Whether a stage completed successfully")
    .def("set_cutoff",             &Topology::set_cutoff,             py::arg("cutoff"),  "Neighbor cutoff used by bond perception (A)")
    .def("set_compute_nonstandard_bonds", &Topology::set_compute_nonstandard_bonds, py::arg("compute"), "Include metal/coordination bonds if True")

    .def("get_atom",  &Topology::atom,  py::arg("idx"), py::return_value_policy::reference_internal, "Atom record at atom index")
    .def("get_ring",  &Topology::ring,  py::arg("idx"), py::return_value_policy::reference_internal, "Ring record at index (0..n_rings-1)")
    .def("get_group", &Topology::group, py::arg("idx"), py::return_value_policy::reference_internal, "Group record at index (0..n_groups-1)")
    .def("residue_index_of_atom", &Topology::residue_index_of_atom, py::arg("atom_idx"), "Return residue index for an atom index")
    .def("residue_of_atom",       &Topology::residue_of_atom,       py::arg("atom_idx"),
      py::return_value_policy::reference_internal, "Return residue object for an atom index")
    .def_property_readonly("atom_to_residue_indices", &Topology::atom_to_residue_indices,
      "Atom-aligned residue index mapping (-1 for unmapped atoms)")

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

    .def("molecule",
         [](const Topology &self) -> const RDKit::RWMol & { return self.molecule(); },
         py::return_value_policy::reference_internal,
         "Underlying RDKit molecule (borrowed reference)")
    .def("conformer",
         [](const Topology &self, int idx) -> const RDKit::Conformer & { return self.molecule().getConformer(idx); },
         py::arg("idx") = 0,
         py::return_value_policy::reference_internal,
         "RDKit conformer by id (borrowed reference)");
}

} // namespace lahuta::bindings
