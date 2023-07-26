"""
# Module: contacts

The `contacts` module is part of the `lahuta` package and is designed to provide 
comprehensive functionality for computing a variety of interactions in molecular 
structures. These include but are not limited to aromatic, hydrophobic, ionic, 
covalent, and Van der Waals interactions. It also provides functionalities to 
compute interactions between atoms and planes, as well as between planes themselves.

**Author**: [Besian I. Sejdiu](../path/to/author/page)

**Institution**: St. Jude Children's Research Hospital

**Contact**: [besian.sejdiu@stjude.org](mailto:besian.sejdiu@stjude.org)

This program is licensed under [GNU General Public License](../path/to/license/file). 
For more details, see also the [GNU General Public License webpage](http://www.gnu.org/licenses/).

The module consists of several submodules, each responsible for computing a specific 
type of interaction. In each submodule, you will find the appropriate documentation 
for the functions provided, including detailed descriptions and usage examples.

"""
from .aromatic import AromaticContacts
from .atom_plane import AtomPlaneContacts, carbon_pi, cation_pi, donor_pi, sulphur_pi
from .carbon_pi import CarbonPi
from .carbonyl import CarbonylContacts
from .cation_pi import CationPi
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
from .donor_pi import DonorPi
from .hbonds import HBondContacts, PolarHBondContacts, WeakHBondContacts, WeakPolarHBondContacts
from .hydrophobic import HydrophobicContacts
from .ionic import IonicContacts
from .metal import MetalicContacts
from .plane_plane import plane_plane_neighbors
from .sulphur_pi import SulphurPi
from .vdw import VanDerWaalsContacts

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
    "donor_pi",
    "sulphur_pi",
    "carbon_pi",
    "cation_pi",
    "plane_plane_neighbors",
    "AromaticContacts",
    "CarbonylContacts",
    "CovalentContacts",
    "HBondContacts",
    "HydrophobicContacts",
    "IonicContacts",
    "MetalicContacts",
    "PolarHBondContacts",
    "VanDerWaalsContacts",
    "WeakHBondContacts",
    "WeakPolarHBondContacts",
    "CarbonPi",
    "CationPi",
    "DonorPi",
    "SulphurPi",
    "AtomPlaneContacts",
]
