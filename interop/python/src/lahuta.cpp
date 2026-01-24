#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "bindings.hpp"
#include "properties/luni_props.hpp"
#include "runtime.hpp"
#include "version.hpp"

// clang-format off
PYBIND11_MODULE(lahuta, m) {
  m.doc() = "lahuta: A Python binding for the Lahuta library";

  // Surface the Lahuta core version in the extension module.
  m.attr("__version__") = std::string(lahuta::version);
  m.attr("__version_info__") = py::make_tuple(
      lahuta::version_major,
      lahuta::version_minor,
      lahuta::version_patch,
      std::string(lahuta::version_suffix)
  );

  // Pre-import numpy so first array creation in property getters doesn't pay
  // the import/initialization cost on the first call (**should** help move cost to module import).
  try { (void)py::module_::import("numpy"); } catch (...) { /**/ }

  // Initialize thread-dependent resources for single thread default
  try { lahuta::LahutaRuntime::ensure_initialized(1); } catch (...) {}

  lahuta::LuniProperties::initialize();

  py::options options;
  options.disable_enum_members_docstring();

  // Binding of modules in logical order is important due to dependencies
  lahuta::bindings::bind_rdkit(m);          // RDKit fundamental classes (Atom, Point3D, etc.)

  lahuta::bindings::bind_atom_types(m);     // AtomType enum and related functions
  lahuta::bindings::bind_elements(m);       // Element helpers
  lahuta::bindings::bind_records(m);        // Entity records (AtomRec, RingRec, GroupRec)
  lahuta::bindings::bind_topology(m);       // Topology and related structures
  lahuta::bindings::bind_entity_id(m);      // EntityID class
  lahuta::bindings::bind_contacts(m);       // Interaction types, Contact/ContactSet and engines
  lahuta::bindings::bind_entities(m);       // FeatureGroup enum, EntityID, SearchOptions, find_contacts
  lahuta::bindings::bind_grid_search(m);    // Grid search class for spatial queries
  lahuta::bindings::bind_distance(m);       // Distance computation and neighbor search
  lahuta::bindings::bind_luni(m);           // Main Luni class

  lahuta::bindings::bind_utilities(m);      // Utility functions (compute_angles, etc.)

  lahuta::bindings::bind_logger(m);
  lahuta::bindings::bind_properties(m);

  lahuta::bindings::bind_pipeline_dynamic(m);
  lahuta::bindings::bind_db(m);
}
