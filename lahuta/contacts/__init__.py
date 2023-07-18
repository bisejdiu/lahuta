from .aromatic import AromaticContacts
from .carbonyl import CarbonylContacts
from .contacts import (
    aromatic_neighbors,
    carbonyl_neighbors,
    covalent_neighbors,
    hbond_neighbors,
    hydrophobic_neighbors,
    ionic_neighbors,
    metalic_neighbors,
    polar_hbond_neighbors,
    vdw_neighbors,
    weak_hbond_neighbors,
    weak_polar_hbond_neighbors,
)
from .covalent import CovalentContacts
from .hbonds import (
    HBondContacts,
    PolarHBondContacts,
    WeakHBondContacts,
    WeakPolarHBondContacts,
)
from .hydrophobic import HydrophobicContacts
from .ionic import IonicContacts
from .metal import MetalicContacts
from .vdw import VanDerWaalsContacts

funcs = [
    aromatic_neighbors,
    carbonyl_neighbors,
    covalent_neighbors,
    hbond_neighbors,
    hydrophobic_neighbors,
    ionic_neighbors,
    metalic_neighbors,
    polar_hbond_neighbors,
    vdw_neighbors,
    weak_hbond_neighbors,
    weak_polar_hbond_neighbors,
]


F = type("F", (), {})
for func in funcs:
    setattr(F, func.__name__, func)

__all__ = [
    "aromatic_neighbors",
    "carbonyl_neighbors",
    "covalent_neighbors",
    "hbond_neighbors",
    "hydrophobic_neighbors",
    "ionic_neighbors",
    "metalic_neighbors",
    "polar_hbond_neighbors",
    "vdw_neighbors",
    "weak_hbond_neighbors",
    "weak_polar_hbond_neighbors",
    "AromaticContacts",
    "CarbonylContacts",
    "CovalentContacts",
    "HBondContacts",
    "HydrophobicContacts",
    "IonicContacts",
    "MetalContacts",
    "PolarHBondContacts",
    "VanDerWaalsContacts",
    "WeakHBondContacts",
    "WeakPolarHBondContacts",
]
