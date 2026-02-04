/**
 * Lahuta - a performant and scalable library for structural biology and bioinformatics
 *
 * Copyright (c) Besian I. Sejdiu (@bisejdiu)
 * License: TBD (see LICENSE file for more info).
 *
 * Contact: std::string{"besian"} + "sejdiu" + "@gmail.com";
 *
 */

#include <sstream>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/Bond.h>
#include <rdkit/GraphMol/Conformer.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/RWMol.h>

#include "numpy_utils.hpp"

namespace py = pybind11;

// clang-format off
namespace lahuta::bindings {
void bind_rdkit(py::module &m) {
  py::module_ rdkit = m.def_submodule("rdkit", "Bindings for RDKit core types");

  py::class_<RDKit::Atom>      Atom(rdkit, "Atom");
  Atom
    .def("getAtomicNum",        &RDKit::Atom::getAtomicNum)
    .def("setAtomicNum",        &RDKit::Atom::setAtomicNum)
    .def("getSymbol",           &RDKit::Atom::getSymbol)
    .def("getIdx",              &RDKit::Atom::getIdx)
    .def("getDegree",           &RDKit::Atom::getDegree)
    .def("calcExplicitValence", &RDKit::Atom::calcExplicitValence, py::arg("strict") = true)
    .def("calcImplicitValence", &RDKit::Atom::calcImplicitValence, py::arg("strict") = true)
    .def("getFormalCharge",     &RDKit::Atom::getFormalCharge)
    .def("getHybridization",    &RDKit::Atom::getHybridization)
    .def("getIsAromatic",       &RDKit::Atom::getIsAromatic)
    .def("getMass",             &RDKit::Atom::getMass)
    .def("getNumExplicitHs",    &RDKit::Atom::getNumExplicitHs)
    .def("getNumImplicitHs",    &RDKit::Atom::getNumImplicitHs)
    .def("getTotalNumHs",       &RDKit::Atom::getTotalNumHs, py::arg("include_neighbors") = false)
    .def("getTotalValence",     &RDKit::Atom::getTotalValence)
    .def("getExplicitValence",  &RDKit::Atom::getExplicitValence)
    .def("getImplicitValence",  &RDKit::Atom::getImplicitValence);

  // Monomer info. Instances are owned by the Atom.
  py::class_<RDKit::AtomMonomerInfo> MonomerInfo(rdkit, "AtomMonomerInfo", "Per-atom monomer metadata; owned by the Atom.");

  py::enum_<RDKit::AtomMonomerInfo::AtomMonomerType>(MonomerInfo, "AtomMonomerType")
    .value("UNKNOWN",    RDKit::AtomMonomerInfo::AtomMonomerType::UNKNOWN)
    .value("PDBRESIDUE", RDKit::AtomMonomerInfo::AtomMonomerType::PDBRESIDUE)
    .value("OTHER",      RDKit::AtomMonomerInfo::AtomMonomerType::OTHER)
    .export_values();

  MonomerInfo
    .def("getName",        &RDKit::AtomMonomerInfo::getName, "Monomer name (e.g., atom label)")
    .def("setName",        &RDKit::AtomMonomerInfo::setName, py::arg("name"), "Set the monomer name")
    .def("getMonomerType", &RDKit::AtomMonomerInfo::getMonomerType, "Kind of monomer annotation")
    .def("setMonomerType", &RDKit::AtomMonomerInfo::setMonomerType)
    .def("copy",           &RDKit::AtomMonomerInfo::copy, py::return_value_policy::take_ownership, "Deep-copy the monomer info");

  py::class_<RDKit::AtomPDBResidueInfo, RDKit::AtomMonomerInfo>(rdkit, "AtomPDBResidueInfo", "PDB residue-level annotation for an atom (serial, resname, chain, etc.).")
    .def(
         py::init([](const std::string &atomName, int serialNumber, const std::string &residueName, int residueNumber, const std::string &chainId,
                    const std::string &altLoc, const std::string &insertionCode, double occupancy, double tempFactor,
                    bool isHet, unsigned int secondaryStructure, unsigned int segmentNumber) {
           return RDKit::AtomPDBResidueInfo(
               atomName, serialNumber, altLoc, residueName, residueNumber, chainId,
               insertionCode, occupancy, tempFactor, isHet, secondaryStructure, segmentNumber);
         }),
         py::arg("atomName"),
         py::arg("serialNumber"),
         py::arg("residueName"),
         py::arg("residueNumber"),
         py::arg("chainId"),
         py::arg("altLoc") = std::string{},
         py::arg("insertionCode") = std::string{},
         py::arg("occupancy") = 1.0,
         py::arg("tempFactor") = 0.0,
         py::arg("isHet") = false,
         py::arg("secondaryStructure") = 0u,
         py::arg("segmentNumber") = 0u,
         "Create a PDB residue annotation (requires 5 fields: atomName, serialNumber, residueName, residueNumber, chainId).")
    .def("getSerialNumber", &RDKit::AtomPDBResidueInfo::getSerialNumber,  "Get the PDB serial number of the atom")
    .def("getAltLoc",       &RDKit::AtomPDBResidueInfo::getAltLoc,        "Get the alternate location identifier (altLoc) for the atom")
    .def("getResidueName",  &RDKit::AtomPDBResidueInfo::getResidueName,   "Get the residue name")
    .def("getResidueIndex", &RDKit::AtomPDBResidueInfo::getResidueIndex,  "Index of the residue in source (non-PDB standard)")
    .def("getResidueNumber",&RDKit::AtomPDBResidueInfo::getResidueNumber, "Get the residue number for the atom")
    .def("getChainId",      &RDKit::AtomPDBResidueInfo::getChainId,       "Get the chain identifier for the atom")
    .def("getInsertionCode",&RDKit::AtomPDBResidueInfo::getInsertionCode, "Get the insertion code for the atom")
    .def("getOccupancy",    &RDKit::AtomPDBResidueInfo::getOccupancy,     "Get the occupancy value for the atom")
    .def("getTempFactor",   &RDKit::AtomPDBResidueInfo::getTempFactor,    "Get the temperature factor (B-factor) for the atom")
    .def("getIsHeteroAtom", &RDKit::AtomPDBResidueInfo::getIsHeteroAtom,  "Check if the atom is a hetero atom (non-standard residue)")
    .def("getSegmentNumber",&RDKit::AtomPDBResidueInfo::getSegmentNumber, "Get the segment number for the atom")
    .def("getSecondaryStructure", &RDKit::AtomPDBResidueInfo::getSecondaryStructure, "Get the secondary structure identifier for the atom")

    .def("setSerialNumber", &RDKit::AtomPDBResidueInfo::setSerialNumber,  py::arg("serialNumber"),  "Set the PDB serial number of the atom")
    .def("setAltLoc",       &RDKit::AtomPDBResidueInfo::setAltLoc,        py::arg("altLoc"),        "Set the alternate location identifier (altLoc) for the atom")
    .def("setResidueName",  &RDKit::AtomPDBResidueInfo::setResidueName,   py::arg("residueName"),   "Set the residue name for the atom")
    .def("setResidueIndex", &RDKit::AtomPDBResidueInfo::setResidueIndex,  py::arg("residueIndex"),  "Set the residue index for the atom (non-PDB standard)")
    .def("setResidueNumber",&RDKit::AtomPDBResidueInfo::setResidueNumber, py::arg("residueNumber"), "Set the residue number for the atom")
    .def("setChainId",      &RDKit::AtomPDBResidueInfo::setChainId,       py::arg("chainId"),       "Set the chain identifier for the atom")
    .def("setInsertionCode",&RDKit::AtomPDBResidueInfo::setInsertionCode, py::arg("insertionCode"), "Set the insertion code for the atom")
    .def("setOccupancy",    &RDKit::AtomPDBResidueInfo::setOccupancy,     py::arg("occupancy"),     "Set the occupancy value for the atom")
    .def("setTempFactor",   &RDKit::AtomPDBResidueInfo::setTempFactor,    py::arg("tempFactor"),    "Set the temperature factor (B-factor) for the atom")
    .def("setIsHeteroAtom", &RDKit::AtomPDBResidueInfo::setIsHeteroAtom,  py::arg("isHeteroAtom"),  "Set whether the atom is a hetero atom (non-standard residue)")
    .def("setSegmentNumber",&RDKit::AtomPDBResidueInfo::setSegmentNumber, py::arg("segmentNumber"), "Set the segment number for the atom")
    .def("setSecondaryStructure", &RDKit::AtomPDBResidueInfo::setSecondaryStructure, py::arg("secondaryStructure"), "Set the secondary structure identifier for the atom");

  // Probably confusing for users. No idea if we should expose this at all.
  // We may wrap this in an RDKit::AtomPDBResidueInfo instead when we return values to the user, to hide the underlying pool mechanics.
  py::class_<RDKit::pAtomPDBResidueInfo, RDKit::AtomPDBResidueInfo>(rdkit, "pAtomPDBResidueInfo", "Pooled AtomPDBResidueInfo optimized for Lahuta's object pool.")
    .def(py::init<>())
    .def(py::init<const char*, int, const char*, int>(), py::arg("atomName"), py::arg("serialNumber"), py::arg("residueName"), py::arg("residueNumber"),
         "Construct with common PDB fields. Other fields default.")
    .def("initialize", &RDKit::pAtomPDBResidueInfo::initialize, py::arg("atomName"), py::arg("serialNumber"), py::arg("residueName"), py::arg("residueNumber"),
         "Reinitialize fields without allocating a new object.")
    .def("resetState", &RDKit::pAtomPDBResidueInfo::resetState, "Reset fields to defaults (pool-friendly).");

  // Expose a view from Atom to its (optional) monomer info.
  // For Python ergonomics, return AtomPDBResidueInfo when present, else None.
  // FIX: we also need to consider the pooled variant here!
  Atom.def("getMonomerInfo",
    [](RDKit::Atom &self) -> RDKit::AtomPDBResidueInfo * {
      auto *info = dynamic_cast<RDKit::AtomPDBResidueInfo *>(self.getMonomerInfo());
      if (!info) return nullptr;
      return info;
    },
    py::return_value_policy::reference,
    "PDB residue annotation for this atom, or None if absent or non-PDB type.");

  py::class_<RDKit::Conformer>(rdkit, "Conformer")
    .def(py::init<>(), "Create an empty conformer")
    .def(py::init<unsigned int>(), py::arg("numAtoms"), "Create with N atoms initialized to (0,0,0)")
    .def("getId",       &RDKit::Conformer::getId, "Conformer identifier")
    .def("setId",       &RDKit::Conformer::setId)
    .def("is3D",        &RDKit::Conformer::is3D, "Whether positions are 3D")
    .def("set3D",       &RDKit::Conformer::set3D)
    .def("getNumAtoms", &RDKit::Conformer::getNumAtoms, "Number of atom positions")
    .def("resize",      &RDKit::Conformer::resize, py::arg("size"), "Resize internal position vector")
    .def("reserve",     &RDKit::Conformer::reserve, py::arg("size"), "Reserve capacity for positions")
    .def("getPositions", [](RDKit::Conformer &self) {
      const auto &coords = self.getPositions();
      if (coords.empty()) {
        return py::array(py::dtype::of<double>(), std::vector<py::ssize_t>{0, 3});
      }
      return lahuta::numpy::make_coordinates_view_f64(coords, py::cast(self));
    }, "Atomic positions as a zero-copy NumPy view")
    .def("getAtomPos",
         [](const RDKit::Conformer &self, unsigned int idx) -> const RDGeom::Point3D & {
           return self.getAtomPos(idx);
         },
         py::arg("idx"), py::return_value_policy::reference,
         "Get atom position by index")
    .def("setAtomPos",
         [](RDKit::Conformer &self, unsigned int idx,
            lahuta::numpy::np_f64 arr) {
           lahuta::numpy::require_shape_1d_len(arr, 3, "Input array must have shape (3,)");
           auto v = arr.unchecked<1>();
           self.setAtomPos(idx, RDGeom::Point3D(v(0), v(1), v(2)));
         }, py::arg("idx"), py::arg("pos"), "Set atom position from numpy array of shape (3,)")
    .def("setAtomPos",
         [](RDKit::Conformer &self, unsigned int idx, const RDGeom::Point3D &pos) {
           self.setAtomPos(idx, pos);
         },
         py::arg("idx"), py::arg("pos"),
         "Set atom position from Point3D")
    .def("setAllAtomPositions",
         [](RDKit::Conformer &self,
            lahuta::numpy::np_f64 arr) {
           lahuta::numpy::require_shape_2d_cols(arr, 3, "Input array must have shape (N,3)");
           auto points = lahuta::numpy::to_point3d_vect(arr);
           RDGeom::POINT3D_VECT &positions = self.getPositions();
           positions = std::move(points);
         }, py::arg("positions"), "Replace all positions from numpy array (N,3)")
    .def("hasOwningMol", &RDKit::Conformer::hasOwningMol, "Whether this conformer is attached to a molecule");

  py::enum_<RDKit::Atom::HybridizationType>(Atom, "HybridizationType")
    .value("UNSPECIFIED", RDKit::Atom::HybridizationType::UNSPECIFIED)
    .value("S",           RDKit::Atom::HybridizationType::S)
    .value("SP",          RDKit::Atom::HybridizationType::SP)
    .value("SP2",         RDKit::Atom::HybridizationType::SP2)
    .value("SP3",         RDKit::Atom::HybridizationType::SP3)
    .value("SP2D",        RDKit::Atom::HybridizationType::SP2D)
    .value("SP3D",        RDKit::Atom::HybridizationType::SP3D)
    .value("SP3D2",       RDKit::Atom::HybridizationType::SP3D2)
    .value("OTHER",       RDKit::Atom::HybridizationType::OTHER)
    .export_values();

  py::enum_<RDKit::Bond::BondType>(rdkit, "BondType")
    .value("UNSPECIFIED", RDKit::Bond::BondType::UNSPECIFIED)
    .value("SINGLE",      RDKit::Bond::BondType::SINGLE)
    .value("DOUBLE",      RDKit::Bond::BondType::DOUBLE)
    .value("TRIPLE",      RDKit::Bond::BondType::TRIPLE)
    .value("AROMATIC",    RDKit::Bond::BondType::AROMATIC)
    .export_values();

  py::enum_<RDKit::Bond::BondDir>(rdkit, "BondDir")
    .value("NONE",         RDKit::Bond::BondDir::NONE)
    .value("BEGINWEDGE",   RDKit::Bond::BondDir::BEGINWEDGE)
    .value("BEGINDASH",    RDKit::Bond::BondDir::BEGINDASH)
    .value("ENDDOWNRIGHT", RDKit::Bond::BondDir::ENDDOWNRIGHT)
    .value("ENDUPRIGHT",   RDKit::Bond::BondDir::ENDUPRIGHT)
    .value("EITHERDOUBLE", RDKit::Bond::BondDir::EITHERDOUBLE)
    .value("UNKNOWN",      RDKit::Bond::BondDir::UNKNOWN)
    .export_values();

  py::enum_<RDKit::Bond::BondStereo>(rdkit, "BondStereo")
    .value("STEREONONE",     RDKit::Bond::BondStereo::STEREONONE)
    .value("STEREOANY",      RDKit::Bond::BondStereo::STEREOANY)
    .value("STEREOZ",        RDKit::Bond::BondStereo::STEREOZ)
    .value("STEREOE",        RDKit::Bond::BondStereo::STEREOE)
    .value("STEREOCIS",      RDKit::Bond::BondStereo::STEREOCIS)
    .value("STEREOTRANS",    RDKit::Bond::BondStereo::STEREOTRANS)
    .value("STEREOATROPCW",  RDKit::Bond::BondStereo::STEREOATROPCW)
    .value("STEREOATROPCCW", RDKit::Bond::BondStereo::STEREOATROPCCW)
    .export_values();

  py::class_<RDKit::Bond> Bond(rdkit, "Bond", "RDKit bond object (owned by a molecule)");
  Bond
    .def("getIdx",          &RDKit::Bond::getIdx, "Index within owning molecule")
    .def("getBondType",     &RDKit::Bond::getBondType, "Bond type (order)")
    .def("setBondType",     &RDKit::Bond::setBondType)
    .def("getBondTypeAsDouble", &RDKit::Bond::getBondTypeAsDouble, "Bond order as float (e.g., 1.0, 1.5, 2.0)")
    .def("getIsAromatic",   &RDKit::Bond::getIsAromatic)
    .def("setIsAromatic",   &RDKit::Bond::setIsAromatic)
    .def("getIsConjugated", &RDKit::Bond::getIsConjugated)
    .def("setIsConjugated", &RDKit::Bond::setIsConjugated)
    .def("getBondDir",      &RDKit::Bond::getBondDir, "Directional annotation")
    .def("setBondDir",      &RDKit::Bond::setBondDir)
    .def("getStereo",       &RDKit::Bond::getStereo, "Stereo annotation")
    .def("setStereo",       &RDKit::Bond::setStereo)
    .def("getBeginAtomIdx", &RDKit::Bond::getBeginAtomIdx)
    .def("getEndAtomIdx",   &RDKit::Bond::getEndAtomIdx)
    .def("getOtherAtomIdx", &RDKit::Bond::getOtherAtomIdx, py::arg("this_idx"))
    .def("getBeginAtom",    &RDKit::Bond::getBeginAtom,    py::return_value_policy::reference, "Begin atom (by reference)")
    .def("getEndAtom",      &RDKit::Bond::getEndAtom,      py::return_value_policy::reference, "End atom (by reference)")
    .def("setStereoAtoms",  &RDKit::Bond::setStereoAtoms,   py::arg("begin_neighbor_idx"), py::arg("end_neighbor_idx"), "Set neighboring atoms used as reference for cis/trans")
    .def("getStereoAtoms", [](RDKit::Bond &self) {
          py::list res;
          for (auto v : self.getStereoAtoms()) res.append(v);
          return res;
        }, "Indices of stereo reference atoms (size 0 or 2)")
    .def("hasOwningMol", &RDKit::Bond::hasOwningMol)
    .def("__str__", [](const RDKit::Bond &b) { std::ostringstream oss; oss << b; return oss.str(); });

  // RWMol useful Python facing API
  py::class_<RDKit::RWMol>(rdkit, "RWMol")
    .def(py::init<>())

    .def("getNumAtoms", [](const RDKit::RWMol &self) { return self.getNumAtoms(); }, "Total number of atoms (explicit only)")
    .def("getNumAtoms", [](const RDKit::RWMol &self, bool only_explicit) { return self.getNumAtoms(only_explicit); }, py::arg("only_explicit"), "Number of atoms; include implicit Hs when only_explicit is False")
    .def("getNumBonds", [](const RDKit::RWMol &self) { return self.getNumBonds(); }, "Number of bonds (heavy atoms only)")
    .def("getNumBonds", [](const RDKit::RWMol &self, bool only_heavy) { return self.getNumBonds(only_heavy); }, py::arg("only_heavy") = true, "Number of bonds; include H-bonds when only_heavy is False")
    .def("getNumConformers", &RDKit::ROMol::getNumConformers, "Number of conformers in the molecule")

    // indexing/iteration
    .def("__getitem__", [](RDKit::RWMol &self, unsigned int idx) -> RDKit::Atom * { return self.getAtomWithIdx(idx); }, py::arg("idx"), py::return_value_policy::reference)
    .def("__iter__", [](RDKit::RWMol &self) {
        py::list res;
        for (auto a : self.atoms()) {
          res.append(py::cast(a, py::return_value_policy::reference));
        }
        return res.attr("__iter__")();
      })
    .def("atoms", [](RDKit::RWMol &self) {
        py::list res;
        for (auto a : self.atoms()) {
          res.append(py::cast(a, py::return_value_policy::reference));
        }
        return res;
      }, "List of Atom objects (by reference)")
    .def("bondPairs", [](RDKit::RWMol &self) {
        py::list res;
        for (auto b : self.bonds()) {
          res.append(py::make_tuple(b->getBeginAtomIdx(), b->getEndAtomIdx()));
        }
        return res;
      }, "List of (begin_idx, end_idx) for all bonds")
    .def("bonds", [](RDKit::RWMol &self) {
        py::list res;
        for (auto b : self.bonds()) {
          res.append(py::cast(b, py::return_value_policy::reference));
        }
        return res;
      }, "List of Bond objects (by reference)")

    .def("addAtom", [](RDKit::RWMol &self, bool update_label) {
          return self.addAtom(update_label);
        },
        py::arg("update_label") = true,
        "Add a new (empty) atom; returns new atom index")
    .def("addAtom", [](RDKit::RWMol &self, int atomic_number, bool update_label) {
          const auto idx = self.addAtom(update_label);
          auto *a = self.getAtomWithIdx(idx);
          a->setAtomicNum(atomic_number);
          return idx;
        },
        py::arg("atomic_number"), py::arg("update_label") = true,
        "Add a new atom with a given atomic number; returns new atom index")
    // Overloads to avoid default of unregistered enum types
    .def("addBond", [](RDKit::RWMol &self, unsigned int begin_idx, unsigned int end_idx) {
        return self.addBond(begin_idx, end_idx, RDKit::Bond::BondType::UNSPECIFIED);
      }, py::arg("begin_idx"), py::arg("end_idx"),
      "Add a bond with unspecified order. Returns new number of bonds")
    .def("addBond", [](RDKit::RWMol &self, unsigned int begin_idx, unsigned int end_idx, RDKit::Bond::BondType order) {
        return self.addBond(begin_idx, end_idx, order);
      }, py::arg("begin_idx"), py::arg("end_idx"), py::arg("order"),
      "Add a bond; returns new number of bonds")
    .def("removeAtom", py::overload_cast<unsigned int>(&RDKit::RWMol::removeAtom), py::arg("idx"), "Remove an atom by index")
    .def("removeBond", &RDKit::RWMol::removeBond, py::arg("begin_idx"), py::arg("end_idx"), "Remove a bond between two atom indices")

    .def("atomNeighbors", [](RDKit::RWMol &self, unsigned int idx) {
        auto *at = self.getAtomWithIdx(idx);
        py::list res;
        for (auto nbr : self.atomNeighbors(at)) {
          res.append(nbr->getIdx());
        }
        return res;
      }, py::arg("idx"),
      "Neighbor atom indices of the given atom")
    .def("atomBonds", [](RDKit::RWMol &self, unsigned int idx) {
        auto *at = self.getAtomWithIdx(idx);
        py::list res;
        for (auto b : self.atomBonds(at)) {
          res.append(py::make_tuple(b->getBeginAtomIdx(), b->getEndAtomIdx()));
        }
        return res;
      }, py::arg("idx"),
      "Bond endpoint index pairs for bonds incident to the given atom")
    .def("atomBondObjects", [](RDKit::RWMol &self, unsigned int idx) {
        auto *at = self.getAtomWithIdx(idx);
        py::list res;
        for (auto b : self.atomBonds(at)) {
          res.append(py::cast(b, py::return_value_policy::reference));
        }
        return res;
      }, py::arg("idx"), "Bond objects incident to the given atom (by reference)")

    .def("getConformer",
         [](RDKit::RWMol &self, int id) -> RDKit::Conformer & { return self.getConformer(id); },
         py::arg("id") = -1, py::return_value_policy::reference,
         "Get a conformer by id (or the first if id < 0)")
    .def("getPositions",
         [](RDKit::RWMol &self, int conf_id) {
           const RDKit::Conformer &conf = self.getConformer(conf_id);
           const auto &coords = conf.getPositions();
           if (coords.empty()) {
             return py::array(py::dtype::of<double>(), std::vector<py::ssize_t>{0, 3});
           }
           return lahuta::numpy::make_coordinates_view_f64(coords, py::cast(conf));
         },
         py::arg("conf_id") = -1,
         "Get atomic positions as a zero-copy NumPy view for the specified conformer")
    .def("addConformer",
         [](RDKit::RWMol &self, const RDKit::Conformer &conf, bool assign_id) {
           auto *dup = new RDKit::Conformer(conf);
           return self.addConformer(dup, assign_id);
         },
         py::arg("conformer"), py::arg("assign_id") = false,
         "Add a conformer (copied) to the molecule; returns its id")
    .def("removeConformer", &RDKit::ROMol::removeConformer, py::arg("id"), "Remove conformer with a given id")
    .def("clearConformers", &RDKit::ROMol::clearConformers, "Remove all conformers")

    .def("getAtomWithIdx",
         [](RDKit::RWMol &self, unsigned int idx) -> RDKit::Atom* {
           return self.getAtomWithIdx(idx);
         }, py::arg("idx"), py::return_value_policy::reference,
         "Get atom by index (by reference)")

    .def("getBondWithIdx",
         [](RDKit::RWMol &self, unsigned int idx) -> RDKit::Bond* {
           return self.getBondWithIdx(idx);
         },
         py::arg("idx"), py::return_value_policy::reference,
         "Get bond by index (by reference)")
    .def("getBondBetweenAtoms",
         [](RDKit::RWMol &self, unsigned int i, unsigned int j) -> RDKit::Bond* {
           return self.getBondBetweenAtoms(i, j);
         }, py::arg("i"), py::arg("j"), py::return_value_policy::reference,
         "Get bond between atom indices i and j, or None")

    .def("updatePropertyCache", &RDKit::ROMol::updatePropertyCache, py::arg("strict")        = true, "Recompute atom/bond property caches")
    .def("clearComputedProps",  &RDKit::ROMol::clearComputedProps,  py::arg("include_rings") = true, "Clear computed properties and optionally ring info")
    .def("getConnectedComponents", [](const RDKit::RWMol &self) {
        std::vector<int> mapping(self.getNumAtoms());
        int n = self.getConnectedComponents(mapping);
        return py::make_tuple(n, mapping);
      }, "Return (num_components, component_ids)")
    .def("debugMol", [](const RDKit::RWMol &self) {
        std::ostringstream oss;
        self.debugMol(oss);
        return oss.str();
      }, "Return a debug string for the molecule");

  py::class_<RDGeom::Point3D>(rdkit, "Point3D")
      .def(py::init<>())
      .def(py::init<double, double, double>())
      .def(py::init([](py::array_t<double, py::array::c_style | py::array::forcecast> coords) {
        if (coords.ndim() != 1 || coords.shape(0) != 3)
          throw std::invalid_argument("Input numpy array must have shape (3)");
        auto v = coords.unchecked<1>();
        return RDGeom::Point3D(v(0), v(1), v(2));
      }))
      .def_readwrite("x", &RDGeom::Point3D::x)
      .def_readwrite("y", &RDGeom::Point3D::y)
      .def_readwrite("z", &RDGeom::Point3D::z)
      .def("__str__",  [](const RDGeom::Point3D &p) { return "(" + std::to_string(p.x) + ", " + std::to_string(p.y) + ", " + std::to_string(p.z) + ")"; });

  rdkit.def("hasNonZeroZCoords", &RDKit::hasNonZeroZCoords, py::arg("conformer"), "Return True if any z coordinate has non-zero magnitude.");
}
} // namespace lahuta::bindings
