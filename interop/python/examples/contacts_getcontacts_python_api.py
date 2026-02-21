# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: Apache License 2.0 (see LICENSE file for more info).
#
# Contact:
#     class Email:
#         def __init__(self, v): self.v = v
#         def bind(self, f): return Email(f(self.v))
#         def get(self): return self.v
#     print(
#         Email("")
#         .bind(lambda s: s + "besian")
#         .bind(lambda s: s + "sejdiu")
#         .bind(lambda s: s + "@gmail.com")
#         .get()
#     )
#
"""Demonstrate getcontacts-style contact detection via Python APIs only."""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, TypeAlias

from lahuta import (
    AtomType,
    AtomTypingMethod,
    ContactSet,
    InteractionType,
    LahutaSystem,
    Topology,
    vdw_radius,
)
from lahuta.entities import atoms, find_contacts, rings
from lahuta.rdkit import Atom, Conformer, Point3D, RWMol

# fmt: off
Vec3D: TypeAlias = tuple[float, float, float]
AtomIndex:  TypeAlias = int
ResidueKey: TypeAlias = tuple[str, int, str]
ChainResidueKey: TypeAlias = tuple[str, int]
DisulfidePair:   TypeAlias = tuple[AtomIndex, AtomIndex]

# distance_sq is the squared distance (A^2).
ContactClassifier:      TypeAlias = Callable[[AtomIndex, AtomIndex, float],               InteractionType | None]
HydrogenBondClassifier: TypeAlias = Callable[[AtomIndex, AtomIndex, float, float | None], InteractionType | None]


def _vector_sub(a: Vec3D, b: Vec3D) -> Vec3D:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _dot(a: Vec3D, b: Vec3D) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _length_sq(vec: Vec3D) -> float:
    return _dot(vec, vec)


def _normalize(vec: Vec3D) -> Vec3D:
    len_sq = _length_sq(vec)
    if len_sq <= 1e-16:  # squared epsilon
        return (0.0, 0.0, 0.0)
    length = math.sqrt(len_sq)
    inv_length = 1.0 / length
    return (vec[0] * inv_length, vec[1] * inv_length, vec[2] * inv_length)


def _clamp(value: float, low: float = -1.0, high: float = 1.0) -> float:
    return max(low, min(high, value))


def _vector_angle_deg(a: Vec3D, b: Vec3D) -> float:
    len_sq_a = _length_sq(a)
    len_sq_b = _length_sq(b)
    if len_sq_a <= 1e-16 or len_sq_b <= 1e-16:
        return 180.0
    inv_len_product = 1.0 / math.sqrt(len_sq_a * len_sq_b)
    cos_val = _clamp(_dot(a, b) * inv_len_product)
    return math.degrees(math.acos(cos_val))


def _psi_angle_deg(center_a: Vec3D, center_b: Vec3D, normal_a: Vec3D) -> float:
    vec   = _vector_sub(center_b, center_a)
    angle = _vector_angle_deg(normal_a, vec)
    return angle if angle <= 90.0 else 180.0 - angle


def _point3d_to_tuple(pt: Point3D) -> Vec3D:
    return (float(pt.x), float(pt.y), float(pt.z))


def _ring_geometries(topology: Topology) -> list[RingGeometry]:
    geoms: list[RingGeometry] = []
    for ring_view in topology.rings:
        center = _point3d_to_tuple(ring_view.center)
        normal = _normalize(_point3d_to_tuple(ring_view.normal))
        atoms = list(ring_view.atoms)
        first_idx = atoms[0].getIdx() if atoms else None
        geoms.append(RingGeometry(center=center, normal=normal, first_atom_idx=first_idx))
    return geoms


class Constants:
    DATA_ROOT: Path = Path(__file__).resolve().parents[3] / "data"
    DEFAULT_STRUCTURE: Path = DATA_ROOT / "ubi.cif"

    POSITIVE_RESIDUE_ATOMS: dict[str, set[str]] = {
        "LYS": {"NZ"},
        "ARG": {"NH1", "NH2"},
    }

    HISTIDINE_RESIDUE_NAMES:  set[str] = {"HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP"}
    HISTIDINE_POSITIVE_ATOMS: set[str] = {"ND1", "NE2"}

    NEGATIVE_RESIDUE_ATOMS: dict[str, set[str]] = {
        "ASP": {"OD1", "OD2"},
        "GLU": {"OE1", "OE2"},
    }

    NUCLEIC_RESIDUES: set[str] = {"A", "C", "G", "U", "DA", "DC", "DG", "DT"}
    PHOSPHATE_ATOMS:  set[str] = {"OP1", "OP2", "O1P", "O2P"}

    STANDARD_PROTEIN_RESIDUES: set[str] = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
        "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
        "TYR", "VAL", "SEC", "PYL",
    }

    CHARGE_SKIP_RESIDUES: set[str] = (
           set(POSITIVE_RESIDUE_ATOMS.keys())
        .union(NEGATIVE_RESIDUE_ATOMS.keys())
        .union(HISTIDINE_RESIDUE_NAMES)
        .union(NUCLEIC_RESIDUES)
        .union(STANDARD_PROTEIN_RESIDUES)
    )


@dataclass(frozen=True, slots=True)
class ResidueInfo:
    """Residue identification information."""
    chain_id:       str
    residue_number: int
    residue_name:   str


@dataclass(slots=True)
class AtomMetadata:
    """Precomputed atom metadata for fast contact detection."""
    atomic_numbers:  list[int]
    vdw_radii:       list[float]
    residue_infos:   list[ResidueInfo | None]
    cys_sg_lookup:   dict[ChainResidueKey, AtomIndex]
    disulfide_pairs: set[DisulfidePair]
    negative_atoms:  set[AtomIndex]
    positive_atoms:  set[AtomIndex]

@dataclass(frozen=True, slots=True)
class RingGeometry:
    """Precomputed ring geometry for pi-stacking calculations."""
    center: Vec3D
    normal: Vec3D
    first_atom_idx: AtomIndex | None


def load_topology(structure: Path) -> tuple[LahutaSystem, Topology]:
    if not structure.exists():
        raise SystemExit(f"Structure file {structure} not found")

    system = LahutaSystem(str(structure))
    if not system.build_topology():
        raise SystemExit("Failed to build topology for the input structure")
    topology = system.get_topology()
    topology.assign_typing(AtomTypingMethod.GetContacts)
    return system, topology


def prepare_atom_metadata(topology: Topology) -> AtomMetadata:
    mol = topology.molecule()
    atomic_numbers: list[int] = []
    vdw_radii:      list[float] = []
    residue_infos:  list[ResidueInfo | None] = []
    cys_sg_lookup:  dict[ChainResidueKey, AtomIndex] = {}
    negative_atoms: set[AtomIndex] = set()
    positive_atoms: set[AtomIndex] = set()

    for atom in mol.atoms():
        idx = atom.getIdx()
        z   = atom.getAtomicNum()
        atomic_numbers.append(z)

        radius = vdw_radius(z) if z else 0.0
        vdw_radii.append(float(radius) if radius else 0.0)

        info = atom.getMonomerInfo()
        if info is not None:
            chain_id = info.getChainId()
            residue_number = info.getResidueNumber()
            resname        = info.getResidueName()
            atom_name      = info.getName()

            residue_infos.append(ResidueInfo(chain_id=chain_id, residue_number=residue_number, residue_name=resname))

            if resname in Constants.POSITIVE_RESIDUE_ATOMS and atom_name in Constants.POSITIVE_RESIDUE_ATOMS[resname]:
                positive_atoms.add(idx)
            elif resname in Constants.HISTIDINE_RESIDUE_NAMES and atom_name in Constants.HISTIDINE_POSITIVE_ATOMS:
                positive_atoms.add(idx)

            if resname in Constants.NEGATIVE_RESIDUE_ATOMS and atom_name in Constants.NEGATIVE_RESIDUE_ATOMS[resname]:
                negative_atoms.add(idx)
            elif resname in Constants.NUCLEIC_RESIDUES and atom_name in Constants.PHOSPHATE_ATOMS:
                negative_atoms.add(idx)

            if atom_name == "SG" and resname == "CYS":
                cys_sg_lookup[(chain_id, residue_number)] = idx

            if resname not in Constants.CHARGE_SKIP_RESIDUES:
                charge = atom.getFormalCharge()
                if charge > 0:
                    positive_atoms.add(idx)
                elif charge < 0:
                    negative_atoms.add(idx)
        else:
            residue_infos.append(None)
            charge = atom.getFormalCharge()
            if charge > 0:
                positive_atoms.add(idx)
            elif charge < 0:
                negative_atoms.add(idx)

    disulfide_pairs: set[DisulfidePair] = set()
    for bond in mol.bonds():
        begin = bond.getBeginAtom()
        end   = bond.getEndAtom()
        if begin.getAtomicNum() == 16 and end.getAtomicNum() == 16:
            disulfide_pairs.add(tuple(sorted((begin.getIdx(), end.getIdx()))))

    return AtomMetadata(
        atomic_numbers=atomic_numbers,
        vdw_radii=vdw_radii,
        residue_infos=residue_infos,
        cys_sg_lookup=cys_sg_lookup,
        disulfide_pairs=disulfide_pairs,
        negative_atoms=negative_atoms,
        positive_atoms=positive_atoms,
    )


def _residue_key(atom: Atom | None) -> ResidueKey | None:
    if atom is None:
        return None
    info = atom.getMonomerInfo()
    if info is None:
        return None
    return (info.getChainId(), info.getResidueNumber(), info.getResidueName().strip())


def _rings_same_residue(mol: RWMol, first_idx_a: AtomIndex | None, first_idx_b: AtomIndex | None) -> bool:
    if first_idx_a is None or first_idx_b is None:
        return False
    atom_a = mol.getAtomWithIdx(first_idx_a)
    atom_b = mol.getAtomWithIdx(first_idx_b)
    key_a = _residue_key(atom_a)
    key_b = _residue_key(atom_b)
    return key_a is not None and key_a == key_b


def _atoms_same_residue(mol: RWMol, idx_a: AtomIndex, idx_b: AtomIndex | None) -> bool:
    if idx_b is None:
        return False
    atom_a = mol.getAtomWithIdx(idx_a)
    atom_b = mol.getAtomWithIdx(idx_b)
    key_a = _residue_key(atom_a)
    key_b = _residue_key(atom_b)
    return key_a is not None and key_a == key_b


def compute_pi_stacking(topology: Topology) -> ContactSet:
    """Compute parallel pi-pi stacking interactions."""
    params_angle_cutoff    = 30.0
    params_psi_cutoff      = 45.0
    params_distance_max    = 10.0
    params_centroid_cutoff = 7.0
    params_centroid_cutoff_sq = params_centroid_cutoff * params_centroid_cutoff

    ring_geoms = _ring_geometries(topology)
    mol = topology.molecule()

    aromatic_ring_selector = rings(lambda r: r.aromatic)

    def tester(idx_a: int, idx_b: int, _d2: float) -> InteractionType:
        geom_a = ring_geoms[idx_a]
        geom_b = ring_geoms[idx_b]

        if _rings_same_residue(mol, geom_a.first_atom_idx, geom_b.first_atom_idx):
            return InteractionType.NoInteraction

        center_a = geom_a.center
        center_b = geom_b.center
        normal_a = geom_a.normal
        normal_b = geom_b.normal

        vec_ab = _vector_sub(center_a, center_b)
        center_distance_sq = _length_sq(vec_ab)
        if center_distance_sq <= 1e-16:
            return InteractionType.NoInteraction
        if center_distance_sq > params_centroid_cutoff_sq:
            return InteractionType.NoInteraction

        normal_angle = _vector_angle_deg(normal_a, normal_b)
        aligned = normal_angle if normal_angle <= 90.0 else 180.0 - normal_angle
        if aligned > params_angle_cutoff:
            return InteractionType.NoInteraction

        psi_a = _psi_angle_deg(center_a, center_b, normal_a)
        psi_b = _psi_angle_deg(center_b, center_a, normal_b)
        if min(psi_a, psi_b) > params_psi_cutoff:
            return InteractionType.NoInteraction

        return InteractionType.PiStackingP

    return find_contacts(topology, aromatic_ring_selector, tester=tester, distance_max=params_distance_max)


def compute_t_stacking(topology: Topology) -> ContactSet:
    """Compute T-shaped pi-pi stacking interactions."""
    params_angle_cutoff  = 30.0
    params_psi_cutoff    = 45.0
    params_distance_max  = 10.0
    params_centroid_cutoff = 5.0
    params_centroid_cutoff_sq = params_centroid_cutoff * params_centroid_cutoff

    ring_geoms = _ring_geometries(topology)
    mol = topology.molecule()
    aromatic_ring_selector = rings(lambda r: r.aromatic)

    def tester(idx_a: int, idx_b: int, _d2: float) -> InteractionType:
        geom_a = ring_geoms[idx_a]
        geom_b = ring_geoms[idx_b]

        if _rings_same_residue(mol, geom_a.first_atom_idx, geom_b.first_atom_idx):
            return InteractionType.NoInteraction

        center_a = geom_a.center
        center_b = geom_b.center
        normal_a = geom_a.normal
        normal_b = geom_b.normal

        vec_ab = _vector_sub(center_a, center_b)
        center_distance_sq = _length_sq(vec_ab)
        if center_distance_sq <= 1e-16:
            return InteractionType.NoInteraction
        if center_distance_sq > params_centroid_cutoff_sq:
            return InteractionType.NoInteraction

        normal_angle = _vector_angle_deg(normal_a, normal_b)
        if abs(normal_angle - 90.0) > params_angle_cutoff:
            return InteractionType.NoInteraction

        psi_a = _psi_angle_deg(center_a, center_b, normal_a)
        psi_b = _psi_angle_deg(center_b, center_a, normal_b)
        if min(psi_a, psi_b) > params_psi_cutoff:
            return InteractionType.NoInteraction

        return InteractionType.PiStackingT

    return find_contacts(topology, aromatic_ring_selector, tester=tester, distance_max=params_distance_max)


def compute_pi_cation(topology: Topology, metadata: AtomMetadata) -> ContactSet:
    """Compute cation-pi interactions."""
    params_centroid_cutoff = 6.0
    params_centroid_cutoff_sq = params_centroid_cutoff * params_centroid_cutoff
    params_angle_cutoff    = 60.0
    params_distance_max    = 10.0

    ring_geoms = _ring_geometries(topology)
    mol  = topology.molecule()
    conf = topology.conformer(0)

    positive_indices = metadata.positive_atoms
    aromatic_ring_selector = rings(lambda r: r.aromatic)
    cation_atom_selector   = atoms(lambda rec: rec.idx in positive_indices)

    def tester(idx_cation: int, idx_ring: int, _d2: float) -> InteractionType:
        geom_ring = ring_geoms[idx_ring]

        if _atoms_same_residue(mol, idx_cation, geom_ring.first_atom_idx):
            return InteractionType.NoInteraction

        center = geom_ring.center
        normal = geom_ring.normal

        cation_pos = _point3d_to_tuple(conf.getAtomPos(idx_cation))
        vec_center_to_cation = _vector_sub(cation_pos, center)

        center_distance_sq = _length_sq(vec_center_to_cation)
        if center_distance_sq <= 1e-16:
            return InteractionType.NoInteraction
        if center_distance_sq > params_centroid_cutoff_sq:
            return InteractionType.NoInteraction

        raw_angle = _vector_angle_deg(normal, vec_center_to_cation)
        angle = raw_angle if raw_angle <= 90.0 else 180.0 - raw_angle
        if angle > params_angle_cutoff:
            return InteractionType.NoInteraction

        return InteractionType.CationPi

    return find_contacts(
        topology,
        cation_atom_selector,
        aromatic_ring_selector,
        tester=tester,
        distance_max=params_distance_max,
    )


def _is_hydrophobic_carbon(mol: RWMol, idx: AtomIndex) -> bool:
    atom = mol.getAtomWithIdx(idx)
    if atom.getAtomicNum() != 6:
        return False
    for nbr_idx in mol.atomNeighbors(idx):
        neighbor = mol.getAtomWithIdx(nbr_idx)
        z = neighbor.getAtomicNum()
        if z not in (6, 1): # carbon or hydrogen
            return False
    return True


def _best_dha_angle_deg(mol: RWMol, conf: Conformer, donor_idx: AtomIndex, acceptor_idx: AtomIndex) -> float | None:
    donor_pos    = _point3d_to_tuple(conf.getAtomPos(donor_idx))
    acceptor_pos = _point3d_to_tuple(conf.getAtomPos(acceptor_idx))

    h_pos_cache = {}
    best_angle: float | None = None
    for neighbor_idx in mol.atomNeighbors(donor_idx):
        if neighbor_idx == acceptor_idx:
            continue
        neighbor_atom = mol.getAtomWithIdx(neighbor_idx)
        if neighbor_atom.getAtomicNum() != 1:
            continue

        if neighbor_idx not in h_pos_cache:
            h_pos_cache[neighbor_idx] = _point3d_to_tuple(conf.getAtomPos(neighbor_idx))
        h_pos = h_pos_cache[neighbor_idx]

        vec_d = _vector_sub(donor_pos, h_pos)
        vec_a = _vector_sub(acceptor_pos, h_pos)

        angle = _vector_angle_deg(vec_d, vec_a)
        if best_angle is None or angle > best_angle:
            best_angle = angle

    return best_angle


def compute_hydrophobic(topology: Topology, metadata: AtomMetadata, classifier: ContactClassifier | None = None) -> ContactSet:
    """Compute hydrophobic contacts between carbon atoms."""
    params_epsilon        = 0.5
    params_min_res_offset = 2
    params_distance_max   = 3.9
    params_distance_max_sq = params_distance_max * params_distance_max

    mol = topology.molecule()
    hydrophobic_selector = atoms(lambda rec: rec.type.has(AtomType.Hydrophobic))

    def tester(idx_a: int, idx_b: int, dist_sq: float) -> InteractionType:
        if residues_too_close(idx_a, idx_b, metadata.residue_infos, params_min_res_offset):
            return InteractionType.NoInteraction

        if not _is_hydrophobic_carbon(mol, idx_a):
            return InteractionType.NoInteraction
        if not _is_hydrophobic_carbon(mol, idx_b):
            return InteractionType.NoInteraction

        if is_cys_disulfide_contact(idx_a, idx_b, metadata.residue_infos, metadata.cys_sg_lookup, mol):
            return InteractionType.NoInteraction

        if dist_sq > params_distance_max_sq:
            return InteractionType.NoInteraction

        cutoff = metadata.vdw_radii[idx_a] + metadata.vdw_radii[idx_b] + params_epsilon
        cutoff_sq = cutoff * cutoff
        if dist_sq >= cutoff_sq:
            return InteractionType.NoInteraction

        if classifier is not None:
            result = classifier(idx_a, idx_b, dist_sq)
            if result is None:
                return InteractionType.NoInteraction
            return result

        return InteractionType.Hydrophobic

    return find_contacts(topology, hydrophobic_selector, tester=tester, distance_max=params_distance_max)


def compute_hbonds(topology: Topology, metadata: AtomMetadata, classifier: HydrogenBondClassifier | None = None) -> ContactSet:
    """Compute hydrogen bonds between donor and acceptor atoms."""
    params_distance_max    = 4.5
    params_distance_cutoff = 3.5
    params_distance_cutoff_sq = params_distance_cutoff * params_distance_cutoff
    params_angle_tolerance    = 180.0
    params_min_residue_offset = 1
    exclude_sulfur = True

    mol  = topology.molecule()
    conf = topology.conformer(0)

    donor_selector    = atoms(lambda rec: rec.type.has(AtomType.HbondDonor))
    acceptor_selector = atoms(lambda rec: rec.type.has(AtomType.HbondAcceptor))

    def tester(idx_d: int, idx_a: int, dist_sq: float) -> InteractionType:
        if residues_too_close(idx_d, idx_a, metadata.residue_infos, params_min_residue_offset):
            return InteractionType.NoInteraction

        atomic_d = metadata.atomic_numbers[idx_d]
        atomic_a = metadata.atomic_numbers[idx_a]
        if exclude_sulfur and (atomic_d == 16 or atomic_a == 16):
            return InteractionType.NoInteraction

        if dist_sq > params_distance_cutoff_sq:
            return InteractionType.NoInteraction

        angle = _best_dha_angle_deg(mol, conf, idx_d, idx_a)
        if angle is not None:
            min_angle = 180.0 - params_angle_tolerance
            if angle < min_angle:
                return InteractionType.NoInteraction

        if classifier is not None:
            result = classifier(idx_d, idx_a, dist_sq, angle)
            if result is None:
                return InteractionType.NoInteraction
            return result

        return InteractionType.HydrogenBond

    return find_contacts(
        topology,
        donor_selector,
        acceptor_selector,
        tester=tester,
        distance_max=params_distance_max,
    )


def residues_too_close(
    idx_a: AtomIndex,
    idx_b: AtomIndex,
    residue_infos: list[ResidueInfo | None],
    min_offset: int,
) -> bool:
    if min_offset <= 0:
        return False

    info_a = residue_infos[idx_a]
    info_b = residue_infos[idx_b]
    if info_a is None or info_b is None:
        return False
    if info_a.chain_id != info_b.chain_id:
        return False

    threshold = min_offset - 1
    return abs(info_a.residue_number - info_b.residue_number) <= threshold


def is_cys_disulfide_contact(
    idx_a: AtomIndex,
    idx_b: AtomIndex,
    residue_infos: list[ResidueInfo | None],
    cys_sg_lookup: dict[ChainResidueKey, AtomIndex],
    mol: RWMol,
) -> bool:
    """Check if two atoms are part of a disulfide bond between CYS residues."""
    info_a = residue_infos[idx_a]
    info_b = residue_infos[idx_b]

    if info_a is None or info_b is None:
        return False

    if info_a.residue_name != "CYS" or info_b.residue_name != "CYS":
        return False

    if info_a.chain_id == info_b.chain_id and info_a.residue_number == info_b.residue_number:
        return False

    key_a = (info_a.chain_id, info_a.residue_number)
    key_b = (info_b.chain_id, info_b.residue_number)
    sg_a = cys_sg_lookup.get(key_a)
    sg_b = cys_sg_lookup.get(key_b)

    if sg_a is None or sg_b is None:
        return False

    return mol.getBondBetweenAtoms(sg_a, sg_b) is not None


def compute_salt_bridges(topology: Topology, meta: AtomMetadata) -> ContactSet:
    """Compute salt bridges between positive and negative atoms."""
    distance_cutoff = 4.0
    distance_cutoff_sq = distance_cutoff * distance_cutoff

    negative_indices = meta.negative_atoms
    positive_indices = meta.positive_atoms

    negative_atom_selector = atoms(lambda rec: rec.idx in negative_indices)
    positive_atom_selector = atoms(lambda rec: rec.idx in positive_indices)

    def tester(idx_anion: int, idx_cation: int, dist_sq: float) -> InteractionType:
        if dist_sq > distance_cutoff_sq:
            return InteractionType.NoInteraction
        if idx_anion == idx_cation:
            return InteractionType.NoInteraction

        pair = tuple(sorted((idx_anion, idx_cation)))
        if pair in meta.disulfide_pairs:
            return InteractionType.NoInteraction

        return InteractionType.Ionic

    return find_contacts(topology, negative_atom_selector, positive_atom_selector, tester=tester, distance_max=5.0)


def compute_vdw_contacts(topology: Topology, meta: AtomMetadata) -> ContactSet:
    """Compute van der Waals contacts between heavy atoms."""
    mol = topology.molecule()
    epsilon = 0.5
    distance_max = 6.0
    min_residue_offset = 2

    heavy_atom_selector = atoms(lambda rec: meta.atomic_numbers[rec.idx] != 1)

    def tester(idx_a: int, idx_b: int, dist_sq: float) -> InteractionType:
        if idx_a == idx_b:
            return InteractionType.NoInteraction

        if meta.atomic_numbers[idx_a] == 1 or meta.atomic_numbers[idx_b] == 1:
            return InteractionType.NoInteraction

        if residues_too_close(idx_a, idx_b, meta.residue_infos, min_residue_offset):
            return InteractionType.NoInteraction

        if is_cys_disulfide_contact(idx_a, idx_b, meta.residue_infos, meta.cys_sg_lookup, mol):
            return InteractionType.NoInteraction

        radius_sum = meta.vdw_radii[idx_a] + meta.vdw_radii[idx_b] + epsilon
        if radius_sum <= 0.0:
            return InteractionType.NoInteraction

        if dist_sq > radius_sum * radius_sum:
            return InteractionType.NoInteraction

        return InteractionType.VanDerWaals

    return find_contacts(topology, heavy_atom_selector, tester=tester, distance_max=distance_max)


def summarize_contacts(
    label: str,
    contacts: ContactSet,
    topology: Topology,
    limit: int,
) -> None:
    print(f"{label}: {contacts.size()} contacts")
    if contacts.empty():
        return

    count = min(limit, contacts.size())
    for idx in range(count):
        contact = contacts[idx]
        print(f"  {contact.describe(topology)}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "structure",
        nargs="?",
        type=Path,
        default=Constants.DEFAULT_STRUCTURE,
        help=f"Input structure readable by Lahuta (default: {Constants.DEFAULT_STRUCTURE.name})",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=10,
        help="Number of exemplar contacts to print for each interaction",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    _, topology = load_topology(args.structure)

    metadata = prepare_atom_metadata(topology)

    salt_bridges = compute_salt_bridges(topology, metadata)
    vdw_contacts = compute_vdw_contacts(topology, metadata)
    hydrophobic  = compute_hydrophobic(topology, metadata)
    pi_cation    = compute_pi_cation(topology, metadata)
    hbonds       = compute_hbonds(topology, metadata)
    pi_parallel  = compute_pi_stacking(topology)
    pi_t         = compute_t_stacking(topology)

    summarize_contacts("Salt bridges",   salt_bridges, topology, args.limit)
    summarize_contacts("Van der Waals",  vdw_contacts, topology, args.limit)
    summarize_contacts("Hydrophobic",    hydrophobic,  topology, args.limit)
    summarize_contacts("Hydrogen bonds", hbonds,       topology, args.limit)
    summarize_contacts("Pi stacking",    pi_parallel,  topology, args.limit)
    summarize_contacts("T stacking",     pi_t,         topology, args.limit)
    summarize_contacts("Pi-cation",      pi_cation,    topology, args.limit)


if __name__ == "__main__":
    main()
