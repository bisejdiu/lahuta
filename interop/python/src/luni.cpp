#include <cstddef>
#include <memory>

#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/RWMol.h>

#include "convert.hpp"
#include "definitions.hpp"
#include "lahuta.hpp"
#include "logging/logging.hpp"
#include "numpy_utils.hpp"
#include "spatial/fastns.hpp"
#include "spatial/nsresults.hpp"
#include "topology.hpp"

namespace py = pybind11;

// clang-format off
namespace lahuta::bindings {

class LuniProperties {
private:
  Luni& luni_;
  py::object owner_;

public:
  LuniProperties(Luni& luni) : luni_(luni), owner_(py::cast(luni)) {}

  // returns a new numpy array a sa copy (float64, shape (n,3))
  auto positions() {
    const RDKit::Conformer &conf = luni_.get_conformer(/*id*/);
    const auto &coords = conf.getPositions();
    const size_t n_atoms = coords.size();
    auto result = py::array_t<double>(std::vector<py::ssize_t>{static_cast<py::ssize_t>(n_atoms), static_cast<py::ssize_t>(3)});
    auto buf = result.request();
    auto *dst = static_cast<double *>(buf.ptr);
    for (size_t i = 0; i < n_atoms; ++i) {
      dst[i * 3 + 0] = coords[i].x;
      dst[i * 3 + 1] = coords[i].y;
      dst[i * 3 + 2] = coords[i].z;
    }
    return result;
  }

  // zero-copy view of coordinates
  auto positions_view() {
    const RDKit::Conformer &conf = luni_.get_conformer(/*id*/);
    const auto &coords = conf.getPositions();
    // Tie lifetime to the parent LahutaSystem py object to avoid
    // re-casting the Conformer on every call (lower overhead, less jitter).
    return numpy::make_coordinates_view_f64(coords, owner_);
  }

  auto indices()     { return numpy::as_numpy_copy(  luni_.indices()); }
  auto atom_nums()   { return numpy::as_numpy_copy(  luni_.atomic_numbers()); }
  auto resids()      { return numpy::as_numpy_copy(  luni_.resids()); }
  auto resindices()  { return numpy::as_numpy_copy(  luni_.resindices()); }
  auto names()       { return numpy::string_array_1d(luni_.names()); }
  auto symbols()     { return numpy::string_array_1d(luni_.symbols()); }
  auto elements()    { return numpy::string_array_1d(luni_.elements()); }
  auto resnames()    { return numpy::string_array_1d(luni_.resnames()); }
  auto chainlabels() { return numpy::string_array_1d(luni_.chainlabels()); }
};

NSResults remove_adjascent_residueid_pairs(const Luni &luni, NSResults &results, int res_diff) {
  Pairs     f_pairs;
  Distances distances;

  for (const auto &[pair, dist] : results) {
    auto *fatom = luni.get_molecule().getAtomWithIdx(pair.first);
    auto *satom = luni.get_molecule().getAtomWithIdx(pair.second);

    auto *finfo = static_cast<const RDKit::AtomPDBResidueInfo *>(fatom->getMonomerInfo());
    auto *sinfo = static_cast<const RDKit::AtomPDBResidueInfo *>(satom->getMonomerInfo());

    if (fatom->getAtomicNum() == 1 || satom->getAtomicNum() == 1) continue;

    auto f_resid = finfo->getResidueNumber();
    auto s_resid = sinfo->getResidueNumber();

    auto is_either_nonprotein = !definitions::is_polymer(finfo->getResidueName()) || !definitions::is_polymer(sinfo->getResidueName());
    if (std::abs(f_resid - s_resid) > res_diff || is_either_nonprotein) {
      f_pairs.push_back(pair);
      distances.push_back(dist);
    }
  }
  return NSResults(f_pairs, distances);
}

void bind_luni(py::module &m) {

  py::class_<IR>(m, "IR", "Intermediate representation for molecular data. All arrays must be 0-based and aligned by atom index.")
    .def(py::init<>(), "Create empty IR")
    .def(py::init<
         std::vector<int>, std::vector<int>, std::vector<std::string>,
         std::vector<int>, std::vector<std::string>, std::vector<std::string>, std::vector<std::vector<float>>>(),
         py::arg("atom_indices"), py::arg("atomic_numbers"), py::arg("atom_names"),
         py::arg("resids"), py::arg("resnames"), py::arg("chainlabels"), py::arg("positions"),
         "Create IR from molecular data arrays.")
    .def_readwrite("atom_indices",   &IR::atom_indices,   "0-based atom indices (0..N-1)")
    .def_readwrite("atomic_numbers", &IR::atomic_numbers, "Atomic numbers (Z)")
    .def_readwrite("atom_names",     &IR::atom_names,     "PDB atom names (e.g., 'CA')")
    .def_readwrite("resids",         &IR::resids,         "Residue sequence identifiers")
    .def_readwrite("resnames",       &IR::resnames,       "Residue names (e.g., 'ALA')")
    .def_readwrite("chainlabels",    &IR::chainlabels,    "Chain labels (e.g., 'A')")
    .def_readwrite("positions",      &IR::positions,      "Cartesian coordinates in A, shape (N, 3)");

  py::class_<LuniProperties>(m, "LahutaSystemProperties", "Wrapper for accessing molecular properties")
    .def_property_readonly("positions",       &LuniProperties::positions,      "Atom coordinates (copy, float64, shape (n,3))")
    .def_property_readonly("positions_view",  &LuniProperties::positions_view, "Atom coordinates view (zero-copy, float64, shape (n,3))")
    .def_property_readonly("indices",     &LuniProperties::indices,     "Atom indices")
    .def_property_readonly("atom_nums",   &LuniProperties::atom_nums,   "Atomic numbers")
    .def_property_readonly("resids",      &LuniProperties::resids,      "Residue IDs")
    .def_property_readonly("resindices",  &LuniProperties::resindices,  "Residue indices")
    .def_property_readonly("names",       &LuniProperties::names,       "Atom names")
    .def_property_readonly("symbols",     &LuniProperties::symbols,     "Element symbols")
    .def_property_readonly("elements",    &LuniProperties::elements,    "Element names")
    .def_property_readonly("resnames",    &LuniProperties::resnames,    "Residue names")
    .def_property_readonly("chainlabels", &LuniProperties::chainlabels, "Chain labels")
    ;

  py::class_<Luni, std::shared_ptr<Luni>> (m, "LahutaSystem", "Main class storing the parsed molecular structure")
    .def(py::init<std::string>(), py::arg("file_name"), "Create a LahutaSystem object from a molecular structure file")
    .def_static("from_model_file", &Luni::from_model_file, py::arg("file_name"), "Create a LahutaSystem object from a model file (fast model pathway)")
    .def_property_readonly("is_model", &Luni::is_model_origin, "Whether the system originated from a model input")

    .def_static("create", py::overload_cast<const IR &>(&Luni::create),                    py::arg("ir"),  "Create from an intermediate representation")
    .def_static("create", py::overload_cast<std::shared_ptr<RDKit::RWMol>>(&Luni::create), py::arg("mol"), "Create from an RDKit molecule")
    // .def_static("create", py::overload_cast<const gemmi::Structure &>(&Luni::create),      py::arg("st"),  "Create from a Gemmi structure")

    .def("build_topology", &Luni::build_topology, py::arg("t_opts") = std::nullopt, "Build the topology with optional configuration")
    .def("has_topology_built", &Luni::has_topology_built, "Check if the topology has been successfully built")

    .def("enable_computation",          &Luni::enable_computation,          py::arg("comp"), py::arg("enabled"), "Enable or disable a specific topology computation")
    .def("enable_only",                 &Luni::enable_only,                 py::arg("comps"),  "Enable only the specified topology computations (disabling all others)")
    .def("is_computation_enabled",      &Luni::is_computation_enabled,      py::arg("comp"),   "Check if a specific topology computation is enabled")
    .def("execute_computation",         &Luni::execute_computation,         py::arg("comp"),   "Execute a specific topology computation with its dependencies")
    .def("set_search_cutoff_for_bonds", &Luni::set_search_cutoff_for_bonds, py::arg("cutoff"), "Set the cutoff distance for neighbor search in bond perception")
    .def("set_atom_typing_method",      &Luni::set_atom_typing_method,      py::arg("method"), "Set the atom typing method (Arpeggio or Molstar)")

    // filtering
    .def("filter", &Luni::filter, py::arg("atom_indices"), "Create a filtered copy of the molecule with specified atom indices")

    .def("get_topology",
       [](Luni &self) -> const Topology& { return *self.get_topology(); },
       py::return_value_policy::reference_internal,
       "Get the topology object (keeps the parent LahutaSystem alive)")

    .def("get_molecule",  &Luni::get_molecule,                      py::return_value_policy::reference_internal, "Get RDKit molecule object")
    .def("get_conformer", &Luni::get_conformer, py::arg("id") = -1, py::return_value_policy::reference_internal, "Get conformer by ID")
    .def("get_atom",      &Luni::get_atom,      py::arg("idx"),     py::return_value_policy::reference_internal, "Get atom by index")

    // TODO: need to bind pdb info object
    // .def("get_info", &Luni::get_info, py::arg("idx"), py::return_value_policy::reference_internal, "Get PDB residue info for atom by index")

    .def_property_readonly("props", [](Luni &luni) {
        return LuniProperties(luni);
    }, py::return_value_policy::reference_internal, "Molecular properties")

    .def_property_readonly("n_atoms",   [](class Luni &luni) { return luni.n_atoms(); },       "Number of atoms")
    .def_property_readonly("file_name", [](class Luni &luni) { return luni.get_file_name(); }, "Source file name")
    // very simple neighbor search
    .def("find_neighbors", [](class Luni &luni, double cutoff, int res_dif){
        auto grid = FastNS(luni.get_conformer().getPositions());
        if (!grid.build(cutoff)) {
          Logger::get_logger()->error("Failed to update the grid with the given cutoff");
          return NSResults();
        }
        auto ns = grid.self_search();
        if (res_dif > 0) {
          ns = remove_adjascent_residueid_pairs(luni, ns, res_dif);
        }
        return ns;
    }, py::arg("cutoff"), py::arg("res_dif") = 0,
    "Find neighboring atoms within cutoff distance. If res_dif > 0, exclude neighbors from the same residue or within res_dif residues.")
    ;
}
} // namespace lahuta::bindings
