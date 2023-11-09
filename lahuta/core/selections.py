"""Selections for the Lahuta package."""

from enum import Enum, auto
from typing import Optional

import numpy as np
from MDAnalysis.core.selection import Selection
from MDAnalysis.core.topologyattrs import Resnames
from numpy.typing import NDArray

from lahuta.analysis.dssp import DSSP, DSSPParser
from lahuta.lahuta_types.mdanalysis import AtomGroupType


class SelectionTokens(Enum):
    """Selection tokens."""

    ACIDIC = auto()
    ALIPHATIC = auto()
    AROMATIC = auto()
    AT = auto()
    BASIC = auto()
    BURIED = auto()
    CG = auto()
    CYCLIC = auto()
    HYDROPHOBIC = auto()
    MEDIUM = auto()
    NEUTRAL = auto()
    PURINE = auto()
    PYRIMIDINE = auto()
    SMALL = auto()
    WATER = auto()


# Mapping of tokens to their respective prot_res values
PROT_RES_MAP = {
    SelectionTokens.ACIDIC: np.array(["ASP", "GLU"]),
    SelectionTokens.ALIPHATIC: np.array(["ALA", "GLY", "ILE", "LEU", "VAL"]),
    SelectionTokens.AROMATIC: np.array(["PHE", "TRP", "TYR", "HIS"]),
    SelectionTokens.AT: np.array(["ADA", "A", "THY", "T"]),
    SelectionTokens.BASIC: np.array(["ARG", "LYS", "HIS"]),
    SelectionTokens.BURIED: np.array(["ALA", "LEU", "VAL", "ILE", "PHE", "CYS", "MET", "TRP"]),
    SelectionTokens.CG: np.array(["CYT", "C", "GUA", "G"]),
    SelectionTokens.CYCLIC: np.array(["HIS", "PHE", "PRO", "TRP", "TYR"]),
    SelectionTokens.HYDROPHOBIC: np.array(["ALA", "LEU", "VAL", "ILE", "PRO", "PHE", "MET", "TRP"]),
    SelectionTokens.MEDIUM: np.array(["VAL", "THR", "ASP", "ASN", "PRO", "CYS", "ASX", "PCA", "HYP"]),
    SelectionTokens.NEUTRAL: np.array(
        ["VAL", "PHE", "GLN", "TYR", "HIS", "CYS", "MET", "TRP", "ASX", "GLX", "PCA", "HYP"]
    ),
    SelectionTokens.PURINE: np.array(["ADE", "A", "GUA", "G"]),
    SelectionTokens.PYRIMIDINE: np.array(["CYT", "C", "THY", "T", "URI", "U"]),
    SelectionTokens.SMALL: np.array(["ALA", "GLY", "SER"]),
    SelectionTokens.WATER: np.array(["H2O", "HH0", "OHH", "HOH", "OH2", "SOL", "WAT", "TIP", "TIP2", "TIP3", "TIP4"]),
}


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


class HELIX_3_10(Selection):  # type: ignore
    """3_10 helix selection."""

    token = "helix_3_10"
    DTYPE = DSSPParser.DTYPES

    def _apply(self, group: AtomGroupType) -> AtomGroupType:
        dssp = DSSP(group.universe.filename)
        dssp.parse_or_calculate_dssp()

        group_resinfo = np.empty(len(group.resnames), dtype=DSSPParser.DTYPES)
        group_resinfo["chain_auths"] = group.chainIDs
        group_resinfo["resname"] = group.resnames
        group_resinfo["resid"] = group.resids

        g_indices = np.where(dssp.ss_array == "G")[0]

        protein_mask = np.isin(group_resinfo, dssp.resinfo_array[g_indices])

        return group[protein_mask]


def create_selection_classes() -> list[type[Selection]]:
    """Create a selection class.

    Iterate over the tokens and create a selection class for each token.

    Returns:
        type[Selection]: The selection class.
    """
    types = []
    for token_enum, prot_res in PROT_RES_MAP.items():
        sel_cls = type(
            f"{token_enum.name.capitalize()}Selection",
            (
                Selection,
                BaseSelection,
            ),
            {"token": token_enum.name.lower(), "prot_res": prot_res},
        )
        types.append(sel_cls)

    return types
