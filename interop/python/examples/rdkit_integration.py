# Lahuta - a performant and scalable library for structural biology and bioinformatics
#
# Copyright (c) Besian I. Sejdiu (@bisejdiu)
# License: TBD (see LICENSE file for more info).
#
# Contact:
#     T = type("", (), {"__getattr__": lambda self, name: {
#         "f": "besian", "l": "sejdiu", "d": "@gmail.com"
#     }[name]})
#     print(T().f + T().l + T().d)
#
"""Example of using (modified) RDKit objects via Lahuta's RDKit integration."""

from pathlib import Path

import numpy as np

from lahuta import LahutaSystem, logging

DATA = Path(__file__).resolve().parents[3] / "core" / "data" / "ubi.cif"


# fmt: off
def rdkit_integration() -> tuple[int, int]:
    sys = LahutaSystem(str(DATA))

    # Accessing RDKit requires building the topology first.
    topo = sys.get_or_build_topology()
    mol  = topo.molecule()
    conf = topo.conformer(0)  # first conformer

    logging.info(f"rdkit: atoms={mol.getNumAtoms()}, conformers={mol.getNumConformers()}, bonds={mol.getNumBonds()}")
    logging.info(f"rdkit: conformer id={conf.getId()}, atoms={conf.getNumAtoms()}, position shape={conf.getPositions().shape}")

    # get monomerinfo for first atom
    at0 = mol.getAtomWithIdx(0)
    info = at0.getMonomerInfo()
    if info:
        logging.info(f"MonomerInfo from Atom: name={info.getName()} res={info.getResidueName()}{info.getResidueNumber()} serial={info.getSerialNumber()}")
    else:
        logging.info("No MonomerInfo attached to the Atom")
    return mol.getNumAtoms(), mol.getNumBonds()


def point3d_and_conformer() -> None:
    from lahuta.rdkit import Conformer, Point3D

    conf = Conformer(numAtoms=3)
    conf.setId(0)
    conf.set3D(True)

    # either Point3D or numpy arrays
    conf.setAtomPos(0, Point3D(0.0, 0.0, 0.0))
    conf.setAtomPos(1, np.array([1.2, 0.0, 0.0], dtype=np.float64))
    conf.setAtomPos(2, np.array([0.0, 1.1, 0.0], dtype=np.float64))

    logging.info(f"Conformer: id={conf.getId()}, num_atoms={conf.getNumAtoms()}, is3D={conf.is3D()}")

    pos = conf.getPositions()
    logging.info(f"Positions view shape={pos.shape}, dtype={pos.dtype} (zero-copy)")


def rwmol_atoms_bonds_and_conformers() -> None:
    from lahuta.rdkit import BondType, Conformer, RWMol

    mol = RWMol()

    a0 = mol.addAtom(6)  # returns atom index
    a1 = mol.addAtom(6)
    a2 = mol.addAtom(8)
    mol.addBond(a0, a1, BondType.SINGLE)
    mol.addBond(a1, a2, BondType.SINGLE)

    logging.info(f"RWMol: num_atoms={mol.getNumAtoms()} num_bonds={mol.getNumBonds()}")

    at0 = mol.getAtomWithIdx(0)
    logging.info(f"Atom[0]: Z={at0.getAtomicNum()} symbol={at0.getSymbol()} degree={at0.getDegree()} aromatic={at0.getIsAromatic()}",)

    if mol.getNumBonds() > 0:
        b0 = mol.getBondWithIdx(0)
        logging.info(f"Bond[0]: ({b0.getBeginAtomIdx()}-{b0.getEndAtomIdx()}) type={b0.getBondType().name}")

    conf = Conformer(mol.getNumAtoms())
    conf.setAllAtomPositions(np.array([[0.1, 0.2, 0.3], [1.4, 0.0, 0.0], [2.8, 0.0, 0.0]], dtype=np.float64))
    cid = mol.addConformer(conf, assign_id=True)
    conf_ref = mol.getConformer(cid)
    logging.info(f"Attached conformer: id={conf_ref.getId()}, num_atoms={conf_ref.getNumAtoms()}")
    logging.info(f"Conformer first coord: {conf_ref.getPositions()[0].tolist()}")

def monomer_info() -> None:
    from lahuta.rdkit import AtomPDBResidueInfo

    info = AtomPDBResidueInfo(atomName="CA", serialNumber=1, residueName="ALA", residueNumber=42, chainId="A", )
    logging.info(f"AtomPDBResidueInfo: name={info.getName()} res={info.getResidueName()}{info.getResidueNumber()} serial={info.getSerialNumber()}",)


if __name__ == "__main__":
    rdkit_integration()
    point3d_and_conformer()
    rwmol_atoms_bonds_and_conformers()
    monomer_info()
