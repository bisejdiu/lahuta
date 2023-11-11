"""Selections for the Lahuta package."""

from enum import Enum
from typing import Optional

import numpy as np
from MDAnalysis.core.selection import Selection
from MDAnalysis.core.topologyattrs import Resnames
from numpy.typing import NDArray

from lahuta.analysis.dssp import DSSP, DSSPParser
from lahuta.config.models.ions import IONS
from lahuta.config.models.lipids import LIPIDS
from lahuta.config.models.sugars import SUGARS
from lahuta.config.models.water import WATERS
from lahuta.lahuta_types.mdanalysis import AtomGroupType


class ResTypeSelectionTokens(Enum):
    """Selection tokens."""

    ACIDIC = ("ASP", "GLU")
    ALIPHATIC = ("ALA", "GLY", "ILE", "LEU", "VAL")
    AROMATIC = ("PHE", "TRP", "TYR", "HIS")
    AT = ("ADA", "A", "THY", "T")
    BASIC = ("ARG", "LYS", "HIS")
    BURIED = ("ALA", "LEU", "VAL", "ILE", "PHE", "CYS", "MET", "TRP")
    CG = ("CYT", "C", "GUA", "G")
    CYCLIC = ("HIS", "PHE", "PRO", "TRP", "TYR")
    HYDROPHOBIC = ("ALA", "LEU", "VAL", "ILE", "PRO", "PHE", "MET", "TRP")
    MEDIUM = ("VAL", "THR", "ASP", "ASN", "PRO", "CYS", "ASX", "PCA", "HYP")
    NEUTRAL = ("VAL", "PHE", "GLN", "TYR", "HIS", "CYS", "MET", "TRP", "ASX", "GLX", "PCA", "HYP")
    PURINE = ("ADE", "A", "GUA", "G")
    PYRIMIDINE = ("CYT", "C", "THY", "T", "URI", "U")
    SMALL = ("ALA", "GLY", "SER")

    # Components
    WATER = tuple(WATERS)
    ION = tuple(IONS)
    LIPID = tuple(LIPIDS)
    SUGAR = tuple(SUGARS)

class SecondaryStructureTokens(Enum):
    """Secondary structure tokens."""

    HELIX = ("H",)
    BETA_BRIDGE = ("B",)
    STRAND = ("E",)
    HELIX_3_10 = ("G",)
    PI_HELIX = ("I",)
    TURN = ("T",)
    BEND = ("S",)
    NO_SS = ("-",)


class BaseSelection:
    """Base class for selections."""

    token: Optional[str] = None
    prot_res: Optional[NDArray[np.str_]] = None

    def _apply(self, group: AtomGroupType) -> AtomGroupType:
        resname_attr = group.universe._topology.resnames  # noqa: SLF001

        matches = self._match_indices(resname_attr)

        # index of each atom's resname
        nmidx = resname_attr.nmidx[group.resindices]

        # intersect atom's resname index and matches to prot_res
        mask: NDArray[np.bool_] = np.in1d(nmidx, matches)
        return group[mask]

    def _match_indices(self, resname_attr: Resnames) -> NDArray[np.int32]:
        """Return indices of residues with recognized residue names."""
        # which values in resname attr are in prot_res?
        assert self.prot_res is not None
        names, indices = np.array(list(resname_attr.namedict.items())).T
        matches: NDArray[np.str_] = indices[np.isin(names, self.prot_res)]
        return matches.astype(int)


class BaseDSSPSelection(Selection):  # type: ignore
    """Base class for helix selections."""

    DTYPE = DSSPParser.DTYPES

    def get_ss_type(self) -> str:
        """Override method in subclasses to define specific secondary structure type."""
        raise NotImplementedError("This method should be overridden in subclasses")

    def _apply(self, group: AtomGroupType) -> AtomGroupType:
        dssp = DSSP(group.universe.filename)
        dssp.parse_or_calculate_dssp()

        group_resinfo = np.empty(len(group.resnames), dtype=self.DTYPE)
        group_resinfo["chain_auths"] = group.chainIDs
        group_resinfo["resname"] = group.resnames
        group_resinfo["resid"] = group.resids

        ss_indices = np.where(dssp.ss_array == self.get_ss_type())[0]
        protein_mask = np.isin(group_resinfo, dssp.resinfo_array[ss_indices])

        return group[protein_mask]


def create_restype_selection_classes() -> list[type[Selection]]:
    """Create a selection class.

    Iterate over the tokens and create a selection class for each token.

    Returns:
        type[Selection]: The selection class.
    """
    types = []
    for token_enum in ResTypeSelectionTokens:
        sel_cls = type(
            f"{token_enum.name.capitalize()}Selection",
            (
                Selection,
                BaseSelection,
            ),
            {"token": token_enum.name.lower(), "prot_res": np.array(token_enum.value, dtype=str)},
        )
        types.append(sel_cls)

    return types


def create_dssp_selection_classes() -> list[type[Selection]]:
    """Create a selection class.

    Iterate over the tokens and create a selection class for each token.

    Returns:
        type[Selection]: The selection class.
    """
    types = []
    for token in SecondaryStructureTokens:
        class_name = token.name.replace("_", " ").title().replace(" ", "")
        sel_cls = type(
            class_name,
            (BaseDSSPSelection,),
            {"token": token.name.lower(), "ss_type": token.value, "get_ss_type": lambda self: self.ss_type},
        )
        types.append(sel_cls)

    return types

class LigandSelection(Selection):
    """Selection for ligands."""

    token: str = "ligand"

    def _apply(self, group: AtomGroupType) -> AtomGroupType:

        exclusion_elements = ["protein", "nucleic", "water", "ion", "lipid", "sugar"]
        return group.select_atoms(f"not ({' or '.join(exclusion_elements)})")
