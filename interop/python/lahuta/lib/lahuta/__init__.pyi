"""
lahuta: A Python binding for the Lahuta library
"""
from __future__ import annotations
import numpy
import typing
from . import db
from . import metrics
from . import neighbors
from . import pipeline
from . import rdkit
__all__: list[str] = ['ArpeggioContactsEngine', 'AtomRec', 'AtomType', 'AtomTypingMethod', 'Category', 'Contact', 'ContactProvider', 'ContactSet', 'Element', 'EntityID', 'EntityResolver', 'FastNS', 'FeatureGroup', 'Flavor', 'GetContactsEngine', 'GroupRec', 'GroupView', 'IR', 'IdentityAnalyzerLuni', 'InputType', 'InteractionType', 'InteractionTypeSet', 'KDIndex', 'Kind', 'LahutaSystem', 'LahutaSystemProperties', 'Logger', 'MolStarContactsEngine', 'NSResults', 'PropertyAnalyzerLuni', 'PropertyKey', 'PropertyQueryLuni', 'Residue', 'Residues', 'RingRec', 'RingView', 'SearchOptions', 'Topology', 'TopologyBuildingOptions', 'TopologyComputers', 'compute_angles', 'covalent_radius', 'db', 'decode_contacts_batch_parallel', 'decode_contacts_binary', 'decode_contacts_binary_columnar', 'decode_contacts_binary_direct', 'factorize', 'find_contacts', 'metrics', 'neighbors', 'pipeline', 'rdkit', 'vdw_radius']
class ArpeggioContactsEngine:
    def __init__(self) -> None:
        ...
    @typing.overload
    def compute(self, topology: Topology) -> ContactSet:
        ...
    @typing.overload
    def compute(self, topology: Topology, only: typing.Any = None) -> ContactSet:
        ...
class AtomRec:
    """
    Atom record with typing and RDKit atom handle
    """
    def __repr__(self) -> str:
        ...
    @property
    def idx(self) -> int:
        """
        RDKit atom index
        """
    @property
    def is_acceptor(self) -> bool:
        """
        True if atom type includes HbondAcceptor
        """
    @property
    def is_aromatic(self) -> bool:
        """
        True if atom type includes Aromatic
        """
    @property
    def is_donor(self) -> bool:
        """
        True if atom type includes HbondDonor
        """
    @property
    def is_hydrophobic(self) -> bool:
        """
        True if atom type includes Hydrophobic
        """
    @property
    def is_negative(self) -> bool:
        """
        True if atom type includes NegativeCharge
        """
    @property
    def is_positive(self) -> bool:
        """
        True if atom type includes PositiveCharge
        """
    @property
    def type(self) -> AtomType:
        """
        Atom type classification
        """
    @type.setter
    def type(self, arg1: AtomType) -> None:
        ...
class AtomType:
    Aromatic: typing.ClassVar[AtomType]  # value = <AtomType.Aromatic: 256>
    CarbonylCarbon: typing.ClassVar[AtomType]  # value = <AtomType.CarbonylCarbon: 128>
    CarbonylOxygen: typing.ClassVar[AtomType]  # value = <AtomType.CarbonylOxygen: 64>
    DativeBondPartner: typing.ClassVar[AtomType]  # value = <AtomType.DativeBondPartner: 8192>
    HbondAcceptor: typing.ClassVar[AtomType]  # value = <AtomType.HbondAcceptor: 1>
    HbondDonor: typing.ClassVar[AtomType]  # value = <AtomType.HbondDonor: 2>
    Hydrophobic: typing.ClassVar[AtomType]  # value = <AtomType.Hydrophobic: 512>
    Invalid: typing.ClassVar[AtomType]  # value = <AtomType.Invalid: 65536>
    IonicTypeMetal: typing.ClassVar[AtomType]  # value = <AtomType.IonicTypeMetal: 32768>
    IonicTypePartner: typing.ClassVar[AtomType]  # value = <AtomType.IonicTypePartner: 4096>
    NegativeCharge: typing.ClassVar[AtomType]  # value = <AtomType.NegativeCharge: 32>
    NoType: typing.ClassVar[AtomType]  # value = <AtomType.NoType: 0>
    PositiveCharge: typing.ClassVar[AtomType]  # value = <AtomType.PositiveCharge: 16>
    TransitionMetal: typing.ClassVar[AtomType]  # value = <AtomType.TransitionMetal: 16384>
    WeakHbondAcceptor: typing.ClassVar[AtomType]  # value = <AtomType.WeakHbondAcceptor: 4>
    WeakHbondDonor: typing.ClassVar[AtomType]  # value = <AtomType.WeakHbondDonor: 8>
    XBondAcceptor: typing.ClassVar[AtomType]  # value = <AtomType.XBondAcceptor: 1024>
    XbondDonor: typing.ClassVar[AtomType]  # value = <AtomType.XbondDonor: 2048>
    __members__: typing.ClassVar[dict[str, AtomType]]  # value = {'NoType': <AtomType.NoType: 0>, 'HbondAcceptor': <AtomType.HbondAcceptor: 1>, 'HbondDonor': <AtomType.HbondDonor: 2>, 'WeakHbondAcceptor': <AtomType.WeakHbondAcceptor: 4>, 'WeakHbondDonor': <AtomType.WeakHbondDonor: 8>, 'PositiveCharge': <AtomType.PositiveCharge: 16>, 'NegativeCharge': <AtomType.NegativeCharge: 32>, 'CarbonylOxygen': <AtomType.CarbonylOxygen: 64>, 'CarbonylCarbon': <AtomType.CarbonylCarbon: 128>, 'Aromatic': <AtomType.Aromatic: 256>, 'Hydrophobic': <AtomType.Hydrophobic: 512>, 'XBondAcceptor': <AtomType.XBondAcceptor: 1024>, 'XbondDonor': <AtomType.XbondDonor: 2048>, 'IonicTypePartner': <AtomType.IonicTypePartner: 4096>, 'DativeBondPartner': <AtomType.DativeBondPartner: 8192>, 'TransitionMetal': <AtomType.TransitionMetal: 16384>, 'IonicTypeMetal': <AtomType.IonicTypeMetal: 32768>, 'Invalid': <AtomType.Invalid: 65536>}
    def __and__(self, arg0: AtomType) -> AtomType:
        ...
    def __contains__(self, flag: AtomType) -> bool:
        """
        Return True if self contains all bits in flag
        """
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __iand__(self, arg0: AtomType) -> AtomType:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ior__(self, arg0: AtomType) -> AtomType:
        ...
    def __ixor__(self, arg0: AtomType) -> AtomType:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __or__(self, arg0: AtomType) -> AtomType:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    def __xor__(self, arg0: AtomType) -> AtomType:
        ...
    def all(self, flags: AtomType) -> bool:
        """
        Return True if self contains all bits in flags
        """
    def any(self, flags: AtomType) -> bool:
        """
        Return True if self contains any bits in flags
        """
    def has(self, flag: AtomType) -> bool:
        """
        Return True if all bits in flag are set in self
        """
    def has_any(self, flags: AtomType) -> bool:
        """
        Return True if any bits in flags are set in self
        """
    def none(self, flags: AtomType) -> bool:
        """
        Return True if self contains none of the bits in flags
        """
    def split(self) -> list[AtomType]:
        """
        Return list of individual flags set in self
        """
    @property
    def empty(self) -> bool:
        """
        True if self has no bits set (AtomType.NoType)
        """
    @property
    def label(self) -> str:
        """
        Human-readable label for this atom type
        """
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class AtomTypingMethod:
    """
    Atom typing backends used when classifying atoms for contacts.
    """
    Arpeggio: typing.ClassVar[AtomTypingMethod]  # value = <AtomTypingMethod.Arpeggio: 0>
    GetContacts: typing.ClassVar[AtomTypingMethod]  # value = <AtomTypingMethod.GetContacts: 2>
    MolStar: typing.ClassVar[AtomTypingMethod]  # value = <AtomTypingMethod.MolStar: 1>
    __members__: typing.ClassVar[dict[str, AtomTypingMethod]]  # value = {'Arpeggio': <AtomTypingMethod.Arpeggio: 0>, 'MolStar': <AtomTypingMethod.MolStar: 1>, 'GetContacts': <AtomTypingMethod.GetContacts: 2>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class Category:
    Aromatic: typing.ClassVar[Category]  # value = <Category.Aromatic: 8>
    CarbonPi: typing.ClassVar[Category]  # value = <Category.CarbonPi: 17>
    Carbonyl: typing.ClassVar[Category]  # value = <Category.Carbonyl: 13>
    CationPi: typing.ClassVar[Category]  # value = <Category.CationPi: 11>
    DonorPi: typing.ClassVar[Category]  # value = <Category.DonorPi: 15>
    Generic: typing.ClassVar[Category]  # value = <Category.Generic: 1>
    Halogen: typing.ClassVar[Category]  # value = <Category.Halogen: 3>
    HydrogenBond: typing.ClassVar[Category]  # value = <Category.HydrogenBond: 4>
    Hydrophobic: typing.ClassVar[Category]  # value = <Category.Hydrophobic: 2>
    Ionic: typing.ClassVar[Category]  # value = <Category.Ionic: 9>
    MetalCoordination: typing.ClassVar[Category]  # value = <Category.MetalCoordination: 10>
    PiStacking: typing.ClassVar[Category]  # value = <Category.PiStacking: 12>
    PolarHydrogenBond: typing.ClassVar[Category]  # value = <Category.PolarHydrogenBond: 6>
    SulphurPi: typing.ClassVar[Category]  # value = <Category.SulphurPi: 16>
    Unclassified: typing.ClassVar[Category]  # value = <Category.Unclassified: 0>
    VanDerWaals: typing.ClassVar[Category]  # value = <Category.VanDerWaals: 14>
    WeakHydrogenBond: typing.ClassVar[Category]  # value = <Category.WeakHydrogenBond: 5>
    WeakPolarHydrogenBond: typing.ClassVar[Category]  # value = <Category.WeakPolarHydrogenBond: 7>
    __members__: typing.ClassVar[dict[str, Category]]  # value = {'Unclassified': <Category.Unclassified: 0>, 'Generic': <Category.Generic: 1>, 'Hydrophobic': <Category.Hydrophobic: 2>, 'Halogen': <Category.Halogen: 3>, 'HydrogenBond': <Category.HydrogenBond: 4>, 'WeakHydrogenBond': <Category.WeakHydrogenBond: 5>, 'PolarHydrogenBond': <Category.PolarHydrogenBond: 6>, 'WeakPolarHydrogenBond': <Category.WeakPolarHydrogenBond: 7>, 'Aromatic': <Category.Aromatic: 8>, 'Ionic': <Category.Ionic: 9>, 'MetalCoordination': <Category.MetalCoordination: 10>, 'CationPi': <Category.CationPi: 11>, 'PiStacking': <Category.PiStacking: 12>, 'Carbonyl': <Category.Carbonyl: 13>, 'VanDerWaals': <Category.VanDerWaals: 14>, 'DonorPi': <Category.DonorPi: 15>, 'SulphurPi': <Category.SulphurPi: 16>, 'CarbonPi': <Category.CarbonPi: 17>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class Contact:
    def __eq__(self, arg0: typing.Any) -> bool:
        ...
    def __hash__(self) -> int:
        ...
    def __init__(self, lhs: EntityID, rhs: EntityID, distance_sq: float, type: InteractionType) -> None:
        ...
    def __lt__(self, arg0: Contact) -> bool:
        ...
    def __ne__(self, arg0: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __str__(self) -> str:
        ...
    def describe(self, topology: Topology) -> str:
        """
        Return a human-readable description using topology resolution.
        """
    @typing.overload
    def to_dict(self) -> dict:
        """
        Return a dict representation suitable for testing/serialization
        """
    @typing.overload
    def to_dict(self, topology: Topology) -> dict:
        """
        Return a human-readable dict representation using topology resolution
        """
    @property
    def distance_sq(self) -> float:
        """
        Distance squared between entities (A^2)
        """
    @distance_sq.setter
    def distance_sq(self, arg1: float) -> None:
        ...
    @property
    def lhs(self) -> EntityID:
        """
        Left entity (EntityID)
        """
    @lhs.setter
    def lhs(self, arg1: EntityID) -> None:
        ...
    @property
    def rhs(self) -> EntityID:
        """
        Right entity (EntityID)
        """
    @rhs.setter
    def rhs(self, arg1: EntityID) -> None:
        ...
    @property
    def type(self) -> InteractionType:
        """
        Interaction type (InteractionType)
        """
    @type.setter
    def type(self, arg1: InteractionType) -> None:
        ...
class ContactProvider:
    Arpeggio: typing.ClassVar[ContactProvider]  # value = <ContactProvider.Arpeggio: 1>
    GetContacts: typing.ClassVar[ContactProvider]  # value = <ContactProvider.GetContacts: 2>
    MolStar: typing.ClassVar[ContactProvider]  # value = <ContactProvider.MolStar: 0>
    __members__: typing.ClassVar[dict[str, ContactProvider]]  # value = {'MolStar': <ContactProvider.MolStar: 0>, 'Arpeggio': <ContactProvider.Arpeggio: 1>, 'GetContacts': <ContactProvider.GetContacts: 2>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class ContactSet:
    def __and__(self, arg0: ContactSet) -> ContactSet:
        ...
    def __getitem__(self, arg0: int) -> Contact:
        ...
    def __iand__(self, arg0: ContactSet) -> ContactSet:
        ...
    def __init__(self) -> None:
        ...
    def __ior__(self, arg0: ContactSet) -> ContactSet:
        ...
    def __isub__(self, arg0: ContactSet) -> ContactSet:
        ...
    def __iter__(self) -> typing.Iterator[Contact]:
        ...
    def __ixor__(self, arg0: ContactSet) -> ContactSet:
        ...
    def __len__(self) -> int:
        ...
    def __or__(self, arg0: ContactSet) -> ContactSet:
        ...
    def __repr__(self) -> str:
        ...
    def __sub__(self, arg0: ContactSet) -> ContactSet:
        ...
    def __xor__(self, arg0: ContactSet) -> ContactSet:
        ...
    def data(self) -> list[Contact]:
        ...
    def empty(self) -> bool:
        ...
    @typing.overload
    def insert(self, contacts: ContactSet) -> None:
        ...
    @typing.overload
    def insert(self, contact: Contact) -> None:
        ...
    def make_generic(self) -> None:
        ...
    def set_difference(self, other: ContactSet) -> ContactSet:
        ...
    def set_intersection(self, other: ContactSet) -> ContactSet:
        ...
    def set_symmetric_difference(self, other: ContactSet) -> ContactSet:
        ...
    def set_union(self, other: ContactSet) -> ContactSet:
        ...
    def size(self) -> int:
        ...
class Element:
    """
    Chemical element descriptor.
    """
    @staticmethod
    def from_symbol(symbol: str) -> Element:
        """
        Create an element from its chemical symbol (case-sensitive).
        """
    def __init__(self, atomic_number: int) -> None:
        """
        Create an element from its atomic number (Z).
        """
    def __repr__(self) -> str:
        ...
    def __str__(self) -> str:
        ...
    @property
    def atomic_number(self) -> int:
        """
        Atomic number (Z).
        """
    @property
    def covalent_radius(self) -> float:
        """
        Covalent radius (A).
        """
    @property
    def name(self) -> str:
        """
        IUPAC element name.
        """
    @property
    def symbol(self) -> str:
        """
        Chemical symbol for the element.
        """
    @property
    def vdw_radius(self) -> float:
        """
        van der Waals radius (A).
        """
class EntityID:
    """
    
    64-bit packed identifier: [ kind:8 | index:56 ].
    
    - `kind` returns the canonical `Kind`
    - `index` is a 32-bit view; `int(self)` yields the raw payload
    - Ordering and hashing use the raw 64-bit value
    
    Example
    >>> from lahuta import Kind, EntityID
    >>> eid = EntityID.make(Kind.Ring, 42)
    >>> int(eid) == (1 << 56) | 42
    True
    >>> eid.kind is Kind.Ring
    True
    >>> eid.index
    42
    >>> str(eid)
    'Ring#42'
    """
    @staticmethod
    def make(kind: Kind, index: int) -> EntityID:
        """
        Pack `(kind, index)` into an `EntityID`.
        
        The raw 56-bit index is preserved; `.index` returns its low 32 bits.
        """
    def __eq__(self, arg0: typing.Any) -> typing.Any:
        """
        Equality by raw payload; returns `NotImplemented` for non-`EntityID` operands.
        """
    def __hash__(self) -> int:
        """
        Hash of the raw 64-bit payload (consistent with equality).
        """
    def __init__(self, raw: int) -> None:
        """
        Construct from raw 64-bit payload `[ kind:8 | index:56 ]`.
        
        No validation. `str(self)` raises `ValueError` if the kind is invalid.
        """
    def __int__(self) -> int:
        """
        Exact 64-bit raw payload (`int`). Layout: `[ kind:8 | index:56 ]`.
        """
    def __lt__(self, arg0: typing.Any) -> typing.Any:
        """
        Ordering by raw payload (equivalently by `(kind, index)`). Returns `NotImplemented` for other types.
        """
    def __ne__(self, arg0: typing.Any) -> typing.Any:
        """
        Inequality by raw payload; returns `NotImplemented` for non-`EntityID` operands.
        """
    def __repr__(self) -> str:
        """
        Debug representation: `"EntityID(Kind#index)"`.
        """
    def __str__(self) -> str:
        """
        Human-readable form: `"Kind#index"`. Raises `ValueError` if the kind is invalid.
        """
    @property
    def index(self) -> int:
        """
        Low 32-bit view of the 56-bit index. The full value is preserved in `int(self)`.
        """
    @property
    def kind(self) -> typing.Any:
        """
        Kind encoded in the upper 8 bits; returns the canonical `Kind` (identity-safe).
        """
class EntityResolver:
    """
    
    Resolve EntityID to concrete records (AtomRec, RingRec, GroupRec) without copying.
    
    Notes
    - The resolver holds a reference to the provided Topology; the topology is kept
      alive for as long as the resolver exists (via keep_alive).
    - Returned records are references to internal storage; do not store them beyond
      the lifetime of the topology.
    """
    def __init__(self, topology: Topology) -> None:
        """
        Create a resolver bound to a specific topology.
        """
    def resolve(self, id: EntityID) -> typing.Any:
        """
        Resolve an EntityID to a concrete record. Returns AtomRec, RingRec, or GroupRec.
        """
    def resolve_all(self, contacts: ContactSet) -> list[tuple[AtomRec | RingRec | GroupRec, AtomRec | RingRec | GroupRec]]:
        """
        Materialize a list of resolved record pairs for a ContactSet.
        """
    def resolve_contact(self, contact: Contact) -> tuple[AtomRec, AtomRec] | tuple[AtomRec, RingRec] | tuple[AtomRec, GroupRec] | tuple[RingRec, AtomRec] | tuple[RingRec, RingRec] | tuple[RingRec, GroupRec] | tuple[GroupRec, AtomRec] | tuple[GroupRec, RingRec] | tuple[GroupRec, GroupRec]:
        """
        Resolve both sides of a contact. Returns a tuple (left_record, right_record).
        """
    @property
    def topology(self) -> Topology:
        """
        Underlying topology (borrowed reference).
        """
class FastNS:
    @staticmethod
    def dist(arg0: float, arg1: float) -> float:
        """
        Return the Euclidean distance between two contiguous float[3] coordinates.
        """
    @staticmethod
    def dist_sq(arg0: float, arg1: float) -> float:
        """
        Return the squared Euclidean distance between two contiguous float[3] coordinates.
        """
    @typing.overload
    def __init__(self, coords: list[rdkit.Point3D]) -> None:
        """
        Construct from a sequence of 3D points.
        
        Args:
            coords: Sequence[Point3D]
        """
    @typing.overload
    def __init__(self, coords: list[list[float]]) -> None:
        """
        Construct from an iterable of (x, y, z) triples.
        
        Args:
            coords: Sequence[Sequence[float]]
        
        Raises:
            ValueError: If any inner sequence does not have length 3.
        """
    @typing.overload
    def __init__(self, coords: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]]) -> None:
        """
        Construct from a NumPy array of 3D coordinates.
        
        Args:
            coords: ndarray of shape (n, 3), dtype=float64. Input coordinates.
        
        Raises:
            ValueError: If coords does not have shape (n, 3).
        """
    def build(self, cutoff: float, brute_force_fallback: bool = True) -> bool:
        """
        Build the neighbor-search grid with the given cutoff.
        
        Args:
            cutoff: Distance threshold used for neighbor detection.
            brute_force_fallback: When True (default), fall back to an internal brute-force search for small systems if the grid cannot be configured.
        
        Returns:
            bool: True if the grid (or fallback) was successfully prepared, False otherwise.
        """
    def get_cutoff(self) -> float:
        """
        Return the current cutoff distance used for neighbor searches.
        """
    def self_search(self) -> NSResults:
        """
        Find all neighbor pairs within the cutoff among the stored coordinates.
        
        Returns:
            NSResults: Neighbor pairs with squared distances.
        """
class FeatureGroup:
    """
    Functional group classification used by the topology
    """
    Acetamidine: typing.ClassVar[FeatureGroup]  # value = <FeatureGroup.Acetamidine: 9>
    Carboxylate: typing.ClassVar[FeatureGroup]  # value = <FeatureGroup.Carboxylate: 10>
    Guanidine: typing.ClassVar[FeatureGroup]  # value = <FeatureGroup.Guanidine: 8>
    Halocarbon: typing.ClassVar[FeatureGroup]  # value = <FeatureGroup.Halocarbon: 7>
    NoGroup: typing.ClassVar[FeatureGroup]  # value = <FeatureGroup.NoGroup: 0>
    Phosphate: typing.ClassVar[FeatureGroup]  # value = <FeatureGroup.Phosphate: 6>
    QuaternaryAmine: typing.ClassVar[FeatureGroup]  # value = <FeatureGroup.QuaternaryAmine: 1>
    Sulfate: typing.ClassVar[FeatureGroup]  # value = <FeatureGroup.Sulfate: 5>
    SulfonicAcid: typing.ClassVar[FeatureGroup]  # value = <FeatureGroup.SulfonicAcid: 4>
    Sulfonium: typing.ClassVar[FeatureGroup]  # value = <FeatureGroup.Sulfonium: 3>
    TertiaryAmine: typing.ClassVar[FeatureGroup]  # value = <FeatureGroup.TertiaryAmine: 2>
    __members__: typing.ClassVar[dict[str, FeatureGroup]]  # value = {'NoGroup': <FeatureGroup.NoGroup: 0>, 'QuaternaryAmine': <FeatureGroup.QuaternaryAmine: 1>, 'TertiaryAmine': <FeatureGroup.TertiaryAmine: 2>, 'Sulfonium': <FeatureGroup.Sulfonium: 3>, 'SulfonicAcid': <FeatureGroup.SulfonicAcid: 4>, 'Sulfate': <FeatureGroup.Sulfate: 5>, 'Phosphate': <FeatureGroup.Phosphate: 6>, 'Halocarbon': <FeatureGroup.Halocarbon: 7>, 'Guanidine': <FeatureGroup.Guanidine: 8>, 'Acetamidine': <FeatureGroup.Acetamidine: 9>, 'Carboxylate': <FeatureGroup.Carboxylate: 10>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class Flavor:
    Default: typing.ClassVar[Flavor]  # value = <Flavor.Default: 0>
    Parallel: typing.ClassVar[Flavor]  # value = <Flavor.Parallel: 1>
    TShape: typing.ClassVar[Flavor]  # value = <Flavor.TShape: 2>
    __members__: typing.ClassVar[dict[str, Flavor]]  # value = {'Default': <Flavor.Default: 0>, 'Parallel': <Flavor.Parallel: 1>, 'TShape': <Flavor.TShape: 2>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class GetContactsEngine:
    def __init__(self) -> None:
        ...
    @typing.overload
    def compute(self, topology: Topology) -> ContactSet:
        ...
    @typing.overload
    def compute(self, topology: Topology, only: typing.Any = None) -> ContactSet:
        ...
class GroupRec:
    """
    Functional group record with atoms
    """
    def __init__(self) -> None:
        """
        Create empty group record
        """
    def __repr__(self) -> str:
        ...
    @property
    def a_type(self) -> AtomType:
        """
        Atom type associated with this group
        """
    @a_type.setter
    def a_type(self, arg1: AtomType) -> None:
        ...
    @property
    def atoms(self) -> list[rdkit.Atom]:
        """
        Atoms participating in the group
        """
    @atoms.setter
    def atoms(self, arg0: list[rdkit.Atom]) -> None:
        ...
    @property
    def type(self) -> FeatureGroup:
        """
        Functional group classification
        """
    @type.setter
    def type(self, arg1: FeatureGroup) -> None:
        ...
class GroupView:
    """
    Functional group view with dynamic geometry accessors
    """
    def __repr__(self) -> str:
        ...
    @property
    def a_type(self) -> AtomType:
        """
        Atom type associated with this group
        """
    @property
    def atoms(self) -> list:
        """
        Atoms participating in the group
        """
    @property
    def center(self) -> rdkit.Point3D:
        """
        Geometric center of the group (A)
        """
    @property
    def type(self) -> FeatureGroup:
        """
        Functional group classification
        """
class IR:
    """
    Intermediate representation for molecular data. All arrays must be 0-based and aligned by atom index.
    """
    @typing.overload
    def __init__(self) -> None:
        """
        Create empty IR
        """
    @typing.overload
    def __init__(self, atom_indices: list[int], atomic_numbers: list[int], atom_names: list[str], resids: list[int], resnames: list[str], chainlabels: list[str], positions: list[list[float]]) -> None:
        """
        Create IR from molecular data arrays.
        """
    @property
    def atom_indices(self) -> list[int]:
        """
        0-based atom indices (0..N-1)
        """
    @atom_indices.setter
    def atom_indices(self, arg0: list[int]) -> None:
        ...
    @property
    def atom_names(self) -> list[str]:
        """
        PDB atom names (e.g., 'CA')
        """
    @atom_names.setter
    def atom_names(self, arg0: list[str]) -> None:
        ...
    @property
    def atomic_numbers(self) -> list[int]:
        """
        Atomic numbers (Z)
        """
    @atomic_numbers.setter
    def atomic_numbers(self, arg0: list[int]) -> None:
        ...
    @property
    def chainlabels(self) -> list[str]:
        """
        Chain labels (e.g., 'A')
        """
    @chainlabels.setter
    def chainlabels(self, arg0: list[str]) -> None:
        ...
    @property
    def positions(self) -> list[list[float]]:
        """
        Cartesian coordinates in A, shape (N, 3)
        """
    @positions.setter
    def positions(self, arg0: list[list[float]]) -> None:
        ...
    @property
    def resids(self) -> list[int]:
        """
        Residue sequence identifiers
        """
    @resids.setter
    def resids(self, arg0: list[int]) -> None:
        ...
    @property
    def resnames(self) -> list[str]:
        """
        Residue names (e.g., 'ALA')
        """
    @resnames.setter
    def resnames(self, arg0: list[str]) -> None:
        ...
class IdentityAnalyzerLuni:
    def __init__(self) -> None:
        ...
class InputType:
    """
    Input categories for LahutaSystem file parsing.
    """
    AlphaFold: typing.ClassVar[InputType]  # value = <InputType.AlphaFold: 1>
    Generic: typing.ClassVar[InputType]  # value = <InputType.Generic: 0>
    __members__: typing.ClassVar[dict[str, InputType]]  # value = {'Generic': <InputType.Generic: 0>, 'AlphaFold': <InputType.AlphaFold: 1>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class InteractionType:
    All: typing.ClassVar[InteractionType]  # value = InteractionType(All)
    Aromatic: typing.ClassVar[InteractionType]  # value = InteractionType(Aromatic)
    CarbonPi: typing.ClassVar[InteractionType]  # value = InteractionType(CarbonPi)
    Carbonyl: typing.ClassVar[InteractionType]  # value = InteractionType(Carbonyl)
    CationPi: typing.ClassVar[InteractionType]  # value = InteractionType(CationPi)
    DonorPi: typing.ClassVar[InteractionType]  # value = InteractionType(DonorPi)
    Generic: typing.ClassVar[InteractionType]  # value = InteractionType(Generic)
    Halogen: typing.ClassVar[InteractionType]  # value = InteractionType(Halogen)
    HydrogenBond: typing.ClassVar[InteractionType]  # value = InteractionType(HydrogenBond)
    Hydrophobic: typing.ClassVar[InteractionType]  # value = InteractionType(Hydrophobic)
    Ionic: typing.ClassVar[InteractionType]  # value = InteractionType(Ionic)
    MetalCoordination: typing.ClassVar[InteractionType]  # value = InteractionType(MetalCoordination)
    NoInteraction: typing.ClassVar[InteractionType]  # value = InteractionType(None)
    PiStacking: typing.ClassVar[InteractionType]  # value = InteractionType(PiStacking)
    PiStackingP: typing.ClassVar[InteractionType]  # value = InteractionType(PiStackingP)
    PiStackingT: typing.ClassVar[InteractionType]  # value = InteractionType(PiStackingT)
    PolarHydrogenBond: typing.ClassVar[InteractionType]  # value = InteractionType(PolarHydrogenBond)
    SulphurPi: typing.ClassVar[InteractionType]  # value = InteractionType(SulphurPi)
    VanDerWaals: typing.ClassVar[InteractionType]  # value = InteractionType(VanDerWaals)
    WeakHydrogenBond: typing.ClassVar[InteractionType]  # value = InteractionType(WeakHydrogenBond)
    WeakPolarHydrogenBond: typing.ClassVar[InteractionType]  # value = InteractionType(WeakPolarHydrogenBond)
    def __eq__(self, arg0: typing.Any) -> typing.Any:
        ...
    def __hash__(self) -> int:
        ...
    def __init__(self, category: Category = ..., flavor: Flavor = ...) -> None:
        ...
    def __int__(self) -> int:
        """
        Return the 32-bit packed code for this InteractionType.
        
        Layout: [ flavor:16 | category:16 ]. Both Category and Flavor are 16-bit enums.
        """
    def __ne__(self, arg0: typing.Any) -> typing.Any:
        ...
    @typing.overload
    def __or__(self, arg0: InteractionType) -> ...:
        ...
    @typing.overload
    def __or__(self, arg0: ...) -> ...:
        ...
    def __repr__(self) -> str:
        ...
    @typing.overload
    def __ror__(self, arg0: InteractionType) -> ...:
        ...
    @typing.overload
    def __ror__(self, arg0: ...) -> ...:
        ...
    def __str__(self) -> str:
        ...
    @property
    def category(self) -> Category:
        ...
    @property
    def flavor(self) -> Flavor:
        ...
class InteractionTypeSet:
    @staticmethod
    def all() -> InteractionTypeSet:
        ...
    def __contains__(self, arg0: InteractionType) -> bool:
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: InteractionType) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: ...) -> None:
        ...
    @typing.overload
    def __ior__(self, arg0: InteractionTypeSet) -> InteractionTypeSet:
        ...
    @typing.overload
    def __ior__(self, arg0: InteractionType) -> InteractionTypeSet:
        ...
    def __len__(self) -> int:
        ...
    @typing.overload
    def __or__(self, arg0: InteractionTypeSet) -> InteractionTypeSet:
        ...
    @typing.overload
    def __or__(self, arg0: InteractionType) -> InteractionTypeSet:
        ...
    def __repr__(self) -> str:
        ...
    @typing.overload
    def __ror__(self, arg0: InteractionTypeSet) -> InteractionTypeSet:
        ...
    @typing.overload
    def __ror__(self, arg0: InteractionType) -> InteractionTypeSet:
        ...
    def __str__(self) -> str:
        ...
    def empty(self) -> bool:
        ...
    def is_all(self) -> bool:
        ...
    def members(self) -> list:
        ...
class KDIndex:
    def __init__(self) -> None:
        """
        Create an empty KDIndex. Call build() before searching.
        """
    @typing.overload
    def build(self, coords: list[rdkit.Point3D]) -> bool:
        """
        Build KD index from a sequence of Point3D.
        """
    @typing.overload
    def build(self, coords: list[list[float]]) -> bool:
        """
        Build KD index from a sequence of (x,y,z) triples.
        """
    @typing.overload
    def build(self, coords: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]]) -> bool:
        """
        Build KD index from an ndarray of shape (n, 3), dtype=float64.
        """
    @typing.overload
    def build_view(self, coords: numpy.ndarray[typing.Any, numpy.dtype[numpy.float32]], leaf_size: int = 40) -> bool:
        """
        Build KD index by viewing a float32 or float64 NumPy array (n,3) without copying. The array must remain alive.
        """
    @typing.overload
    def build_view(self, coords: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]], leaf_size: int = 40) -> bool:
        """
        Build KD index by viewing a float32 or float64 NumPy array (n,3) without copying. The array must remain alive.
        """
    def radius_neighbors_flat(self, queries: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]], radius: float, return_distance: bool = False, sort_results: bool = False) -> tuple:
        """
        Return neighbors in CSR-like flat arrays: (indices, indptr) or (distances, indices, indptr).
        """
    @typing.overload
    def radius_search(self, queries: list[rdkit.Point3D], radius: float, grouped: bool = False, return_distance: bool = False, sort_results: bool = False) -> typing.Any:
        """
        Search neighbors within radius.
        
        Args:
            queries: Query coordinates.
            radius: Search radius.
            grouped: If True, return per-query lists.
            return_distance: If grouped is True, include distances.
            sort_results: If grouped is True, sort neighbors per query.
        
        Returns:
            NSResults if grouped is False.
            list or (distances, indices) if grouped is True.
        """
    @typing.overload
    def radius_search(self, queries: numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]], radius: float, grouped: bool = False, return_distance: bool = False, sort_results: bool = False) -> typing.Any:
        """
        Search neighbors within radius.
        
        Args:
            queries: Query coordinates.
            radius: Search radius.
            grouped: If True, return per-query lists.
            return_distance: If grouped is True, include distances.
            sort_results: If grouped is True, sort neighbors per query.
        
        Returns:
            NSResults if grouped is False.
            list or (distances, indices) if grouped is True.
        """
    @property
    def ready(self) -> bool:
        """
        Return True if the KD index is built and ready.
        """
class Kind:
    """
    
    Entity category (Atom=0, Ring=1, Group=2).
    
    Stored in the upper 8 bits of an EntityID and used for ordering.
    """
    Atom: typing.ClassVar[Kind]  # value = <Kind.Atom: 0>
    Group: typing.ClassVar[Kind]  # value = <Kind.Group: 2>
    Ring: typing.ClassVar[Kind]  # value = <Kind.Ring: 1>
    __members__: typing.ClassVar[dict[str, Kind]]  # value = {'Atom': <Kind.Atom: 0>, 'Ring': <Kind.Ring: 1>, 'Group': <Kind.Group: 2>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class LahutaSystem:
    """
    Main class storing the parsed molecular structure
    """
    @staticmethod
    @typing.overload
    def create(ir: IR) -> LahutaSystem:
        """
        Create from an intermediate representation
        """
    @staticmethod
    @typing.overload
    def create(mol: rdkit.RWMol) -> LahutaSystem:
        """
        Create from an RDKit molecule
        """
    def __init__(self, file_name: str, *, input_type: InputType = ...) -> None:
        """
        Create a LahutaSystem object from a molecular structure file. Set input_type for model inputs.
        """
    @typing.overload
    def build_topology(self, t_opts: TopologyBuildingOptions | None = None) -> bool:
        """
        Build the topology with optional configuration
        """
    @typing.overload
    def build_topology(self, t_opts: TopologyBuildingOptions, include: TopologyComputers) -> bool:
        """
        Build topology and ensure requested computations
        """
    @typing.overload
    def build_topology(self, include: TopologyComputers) -> bool:
        """
        Build topology with default options and ensure requested computations
        """
    @typing.overload
    def build_topology(self, t_opts: TopologyBuildingOptions, include: typing.Iterable) -> bool:
        """
        Build topology and ensure requested computations (list of flags)
        """
    @typing.overload
    def build_topology(self, include: typing.Iterable) -> bool:
        """
        Build topology with default options and ensure requested computations (list of flags)
        """
    def filter(self, atom_indices: list[int]) -> LahutaSystem:
        """
        Create a filtered copy of the molecule with specified atom indices
        """
    def find_neighbors(self, cutoff: float, residue_difference: int = 0) -> NSResults:
        """
        Find neighboring atoms within cutoff distance. If residue_difference > 0, exclude neighbors from the same residue or within residue_difference residues.
        """
    def get_atom(self, idx: int) -> rdkit.Atom:
        """
        Get atom by index
        """
    def get_or_build_topology(self, t_opts: TopologyBuildingOptions | None = None) -> Topology:
        """
        Return topology if built, otherwise build with optional options
        """
    def get_topology(self) -> Topology:
        """
        Get the topology object (keeps the parent LahutaSystem alive)
        """
    def has_topology_built(self) -> bool:
        """
        Check if the topology has been successfully built
        """
    def reset_topology(self) -> LahutaSystem:
        """
        Return a fresh system by reloading the original input file
        """
    def set_search_cutoff_for_bonds(self, cutoff: float) -> None:
        """
        Set the cutoff distance for neighbor search in bond perception
        """
    @property
    def file_name(self) -> str:
        """
        Source file name
        """
    @property
    def is_model(self) -> bool:
        """
        Whether the system originated from a model input
        """
    @property
    def n_atoms(self) -> int:
        """
        Number of atoms
        """
    @property
    def props(self) -> LahutaSystemProperties:
        """
        Molecular properties
        """
    @property
    def residues(self) -> Residues:
        """
        Residues container
        """
class LahutaSystemProperties:
    """
    Wrapper for accessing molecular properties
    """
    @property
    def atom_nums(self) -> numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]]:
        """
        Atomic numbers
        """
    @property
    def chainlabels(self) -> numpy.ndarray[typing.Any, numpy.dtype[typing.Any]]:
        """
        Chain labels
        """
    @property
    def elements(self) -> numpy.ndarray[typing.Any, numpy.dtype[typing.Any]]:
        """
        Element names
        """
    @property
    def indices(self) -> numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]]:
        """
        Atom indices
        """
    @property
    def names(self) -> numpy.ndarray[typing.Any, numpy.dtype[typing.Any]]:
        """
        Atom names
        """
    @property
    def positions(self) -> numpy.ndarray[typing.Any, numpy.dtype[numpy.float64]]:
        """
        Atom coordinates (copy, float64, shape (n,3))
        """
    @property
    def positions_view(self) -> numpy.ndarray[typing.Any, numpy.dtype[typing.Any]]:
        """
        Atom coordinates view (zero-copy, float64, shape (n,3))
        """
    @property
    def resids(self) -> numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]]:
        """
        Residue IDs
        """
    @property
    def resindices(self) -> numpy.ndarray[typing.Any, numpy.dtype[typing.Any]]:
        """
        Residue indices (AtomPDBResidueInfo.residueIndex)
        """
    @property
    def resnames(self) -> numpy.ndarray[typing.Any, numpy.dtype[typing.Any]]:
        """
        Residue names
        """
    @property
    def symbols(self) -> numpy.ndarray[typing.Any, numpy.dtype[typing.Any]]:
        """
        Element symbols
        """
class Logger:
    class FormatStyle:
        Detailed: typing.ClassVar[Logger.FormatStyle]  # value = <FormatStyle.Detailed: 1>
        Simple: typing.ClassVar[Logger.FormatStyle]  # value = <FormatStyle.Simple: 0>
        __members__: typing.ClassVar[dict[str, Logger.FormatStyle]]  # value = {'Simple': <FormatStyle.Simple: 0>, 'Detailed': <FormatStyle.Detailed: 1>}
        def __eq__(self, other: typing.Any) -> bool:
            ...
        def __getstate__(self) -> int:
            ...
        def __hash__(self) -> int:
            ...
        def __index__(self) -> int:
            ...
        def __init__(self, value: int) -> None:
            ...
        def __int__(self) -> int:
            ...
        def __ne__(self, other: typing.Any) -> bool:
            ...
        def __repr__(self) -> str:
            ...
        def __setstate__(self, state: int) -> None:
            ...
        def __str__(self) -> str:
            ...
        @property
        def name(self) -> str:
            ...
        @property
        def value(self) -> int:
            ...
    class LogLevel:
        Critical: typing.ClassVar[Logger.LogLevel]  # value = <LogLevel.Critical: 5>
        Debug: typing.ClassVar[Logger.LogLevel]  # value = <LogLevel.Debug: 1>
        Error: typing.ClassVar[Logger.LogLevel]  # value = <LogLevel.Error: 4>
        Info: typing.ClassVar[Logger.LogLevel]  # value = <LogLevel.Info: 2>
        Off: typing.ClassVar[Logger.LogLevel]  # value = <LogLevel.Off: 6>
        Trace: typing.ClassVar[Logger.LogLevel]  # value = <LogLevel.Trace: 0>
        Warn: typing.ClassVar[Logger.LogLevel]  # value = <LogLevel.Warn: 3>
        __members__: typing.ClassVar[dict[str, Logger.LogLevel]]  # value = {'Trace': <LogLevel.Trace: 0>, 'Debug': <LogLevel.Debug: 1>, 'Info': <LogLevel.Info: 2>, 'Warn': <LogLevel.Warn: 3>, 'Error': <LogLevel.Error: 4>, 'Critical': <LogLevel.Critical: 5>, 'Off': <LogLevel.Off: 6>}
        def __eq__(self, other: typing.Any) -> bool:
            ...
        def __getstate__(self) -> int:
            ...
        def __hash__(self) -> int:
            ...
        def __index__(self) -> int:
            ...
        def __init__(self, value: int) -> None:
            ...
        def __int__(self) -> int:
            ...
        def __ne__(self, other: typing.Any) -> bool:
            ...
        def __repr__(self) -> str:
            ...
        def __setstate__(self, state: int) -> None:
            ...
        def __str__(self) -> str:
            ...
        @property
        def name(self) -> str:
            ...
        @property
        def value(self) -> int:
            ...
    @staticmethod
    def get_instance() -> Logger:
        ...
    def log(self, arg0: Logger.LogLevel, arg1: str) -> None:
        ...
    def set_format(self, arg0: Logger.FormatStyle) -> None:
        ...
    def set_log_level(self, arg0: Logger.LogLevel) -> None:
        ...
class MolStarContactsEngine:
    def __init__(self) -> None:
        ...
    @typing.overload
    def compute(self, topology: Topology) -> ContactSet:
        ...
    @typing.overload
    def compute(self, topology: Topology, only: typing.Any = None) -> ContactSet:
        ...
class NSResults:
    @typing.overload
    def __init__(self) -> None:
        """
        Create an empty NSResults object.
        """
    @typing.overload
    def __init__(self, other: NSResults) -> None:
        """
        Create a new NSResults as a copy of another.
        """
    @typing.overload
    def __init__(self, pairs: numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]], distances: numpy.ndarray[typing.Any, numpy.dtype[numpy.float32]]) -> None:
        """
        Create an NSResults from NumPy arrays.
        
        Args:
            pairs: ndarray of shape (n, 2), dtype=int. Each row is a neighbor pair (i, j).
            distances: ndarray of shape (n,), dtype=float32. Squared distances corresponding to each pair.
        
        Raises:
            ValueError: If shapes do not match or arrays are malformed.
        """
    def __iter__(self) -> typing.Iterator[tuple[tuple[int, int], float]]:
        """
        Iterate over ((i, j), distance_sq) tuples.
        """
    def __len__(self) -> int:
        """
        Number of stored neighbor pairs.
        """
    def add_neighbors(self, i: int, j: int, distance_sq: float) -> None:
        """
        Append a neighbor pair (i, j) with its squared distance to this results container.
        
        Args:
            i: Index of the first item.
            j: Index of the second item.
            distance_sq: Squared distance between the two items.
        """
    def clear(self) -> None:
        """
        Remove all stored pairs and distances.
        """
    @typing.overload
    def filter(self, cutoff: float) -> NSResults:
        """
        Return a new NSResults containing only pairs with squared distance <= cutoff**2.
        
        Args:
            cutoff: Distance threshold.
        
        Returns:
            NSResults: Filtered results (copy).
        """
    @typing.overload
    def filter(self, indices: list[int]) -> NSResults:
        """
        Return a new NSResults containing only pairs where either index is in indices.
        
        Args:
            indices: Set of indices to keep.
        
        Returns:
            NSResults: Filtered results (copy).
        """
    @typing.overload
    def filter(self, indices: list[int], offset: int) -> NSResults:
        """
        Return a new NSResults filtered by a specific column.
        
        Args:
            indices: Indices to keep.
            offset: Column selector: 0 to match the first element of each pair, 1 to match the second.
        
        Returns:
            NSResults: Filtered results (copy).
        
        Raises:
            ValueError: If offset is not 0 or 1.
        """
    def get_distances(self) -> numpy.ndarray[typing.Any, numpy.dtype[numpy.float32]]:
        """
        Alias of 'distances': squared distances (copy)
        """
    def get_sqrt_distances(self) -> numpy.ndarray[typing.Any, numpy.dtype[numpy.float32]]:
        """
        Return the square root of the stored squared distances as a new (n,) float32 NumPy array.
        
        Notes:
          - This computes the square root, but does NOT modify internal data.
          - Very small negative values (< -1e-8f) caused by floating-point round-off are clamped to 0.
        """
    def size(self) -> int:
        """
        Return the number of stored neighbor pairs.
        """
    @property
    def distances(self) -> numpy.ndarray[typing.Any, numpy.dtype[numpy.float32]]:
        """
        Squared distances as a new (n,) float array (copy)
        """
    @property
    def distances_sq(self) -> numpy.ndarray[typing.Any, numpy.dtype[numpy.float32]]:
        """
        Alias of 'distances': squared distances (copy)
        """
    @property
    def distances_view(self) -> numpy.ndarray[typing.Any, numpy.dtype[typing.Any]]:
        """
        Zero-copy read-only view of squared distances with shape (n,)
        """
    @property
    def pairs(self) -> numpy.ndarray[typing.Any, numpy.dtype[numpy.int32]]:
        """
        Neighbor index pairs as a new (n,2) int array (copy)
        """
    @property
    def pairs_view(self) -> numpy.ndarray[typing.Any, numpy.dtype[typing.Any]]:
        """
        Zero-copy read-only view of pairs with shape (n,2) and custom strides
        """
class PropertyAnalyzerLuni:
    def __init__(self, arg0: PropertyQueryLuni) -> None:
        ...
class PropertyKey:
    Elements: typing.ClassVar[PropertyKey]  # value = <PropertyKey.Elements: 2>
    Indices: typing.ClassVar[PropertyKey]  # value = <PropertyKey.Indices: 1>
    Names: typing.ClassVar[PropertyKey]  # value = <PropertyKey.Names: 0>
    Positions: typing.ClassVar[PropertyKey]  # value = <PropertyKey.Positions: 3>
    __members__: typing.ClassVar[dict[str, PropertyKey]]  # value = {'Names': <PropertyKey.Names: 0>, 'Indices': <PropertyKey.Indices: 1>, 'Elements': <PropertyKey.Elements: 2>, 'Positions': <PropertyKey.Positions: 3>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class PropertyQueryLuni:
    def __init__(self) -> None:
        ...
    def properties(self) -> list[PropertyKey]:
        ...
    @typing.overload
    def select(self, arg0: PropertyKey) -> PropertyQueryLuni:
        ...
    @typing.overload
    def select(self, arg0: list[PropertyKey]) -> PropertyQueryLuni:
        ...
class Residue:
    """
    Residue view with identity and atom membership. Comparable and hashable.
    """
    def __eq__(self, arg0: typing.Any) -> typing.Any:
        ...
    def __hash__(self) -> int:
        ...
    @typing.overload
    def __init__(self) -> None:
        """
        Create empty residue
        """
    @typing.overload
    def __init__(self, chain_id: str, number: int, name: str, alt_loc: str) -> None:
        """
        Create residue with specified properties
        """
    def __ne__(self, arg0: typing.Any) -> typing.Any:
        ...
    def __repr__(self) -> str:
        ...
    def __str__(self) -> str:
        ...
    @property
    def alt_loc(self) -> str:
        """
        Alternative location identifier
        """
    @alt_loc.setter
    def alt_loc(self, arg0: str) -> None:
        ...
    @property
    def atoms(self) -> list[rdkit.Atom]:
        """
        Atom indices in this residue
        """
    @atoms.setter
    def atoms(self, arg0: list[rdkit.Atom]) -> None:
        ...
    @property
    def chain_id(self) -> str:
        """
        Chain identifier
        """
    @chain_id.setter
    def chain_id(self, arg0: str) -> None:
        ...
    @property
    def idx(self) -> int:
        """
        0-based index of this residue in the parent molecule
        """
    @idx.setter
    def idx(self, arg0: int) -> None:
        ...
    @property
    def name(self) -> str:
        """
        Residue name
        """
    @name.setter
    def name(self, arg0: str) -> None:
        ...
    @property
    def number(self) -> int:
        """
        Residue number
        """
    @number.setter
    def number(self, arg0: int) -> None:
        ...
class Residues:
    """
    Random-access container over residues with functional helpers.
    """
    def __getitem__(self, arg0: int) -> Residue:
        """
        Get residue by index (view kept alive by the container)
        """
    def __init__(self, mol: rdkit.RWMol) -> None:
        """
        Create residues from RDKit molecule
        """
    def __iter__(self) -> typing.Iterator[Residue]:
        """
        Iterate over residues
        """
    def __len__(self) -> int:
        """
        Get number of residues
        """
    def filter(self, func: typing.Callable) -> Residues:
        """
        Return a new Residues filtered by predicate
        """
    def get_atom_ids(self) -> list[int]:
        """
        Concatenated atom indices of all residues (order-preserving)
        """
    def get_residue_names(self) -> list[str]:
        """
        Residue names aligned to container order
        """
    def map(self, func: typing.Callable) -> list[typing.Any]:
        """
        Apply function to each residue, returning a Python list
        """
    def residue_index_of_atom(self, atom_idx: int) -> int:
        """
        Return residue index for an atom index
        """
    def residue_of_atom(self, atom_idx: int) -> Residue:
        """
        Return residue object for an atom index
        """
    @property
    def atom_to_residue_indices(self) -> list[int]:
        """
        Atom-aligned residue index mapping (-1 for unmapped atoms)
        """
    @property
    def residues(self) -> list[Residue]:
        """
        Snapshot list of residues
        """
class RingRec:
    """
    Ring record with atom membership
    """
    def __init__(self) -> None:
        """
        Create empty ring record
        """
    def __repr__(self) -> str:
        ...
    @property
    def aromatic(self) -> bool:
        """
        Whether the ring is aromatic
        """
    @aromatic.setter
    def aromatic(self, arg0: bool) -> None:
        ...
    @property
    def atoms(self) -> list[rdkit.Atom]:
        """
        Atoms participating in the ring
        """
    @atoms.setter
    def atoms(self, arg0: list[rdkit.Atom]) -> None:
        ...
    @property
    def size(self) -> int:
        """
        Number of atoms in the ring
        """
class RingView:
    """
    Ring view with dynamic geometry accessors
    """
    def __repr__(self) -> str:
        ...
    @property
    def aromatic(self) -> bool:
        """
        Whether the ring is aromatic
        """
    @property
    def atoms(self) -> list:
        """
        Atoms participating in the ring
        """
    @property
    def center(self) -> rdkit.Point3D:
        """
        Ring geometric center (A)
        """
    @property
    def normal(self) -> rdkit.Point3D:
        """
        Normal vector to ring plane
        """
    @property
    def size(self) -> int:
        """
        Number of atoms in the ring
        """
class SearchOptions:
    distance_max: float
    def __init__(self) -> None:
        ...
class Topology:
    """
    Topology manager built over an RDKit molecule and conformer.
    """
    def assign_typing(self, method: AtomTypingMethod) -> None:
        """
        Assign per-atom types using specified method; populates atom_records
        """
    def atom_types_filter_by_fn(self, func: typing.Callable) -> list[AtomRec]:
        """
        Return atom records for which the predicate returns True
        """
    def atoms_with_type(self, type: AtomType) -> list[AtomRec]:
        """
        Return atom records containing the specified type flag
        """
    def build(self, options: TopologyBuildingOptions) -> bool:
        """
        Build all enabled stages. Returns a boolean
        """
    def conformer(self, idx: int = 0) -> rdkit.Conformer:
        """
        RDKit conformer by id (borrowed reference)
        """
    def get_atom(self, idx: int) -> AtomRec:
        """
        Atom record at atom index
        """
    def get_atom_ids(self) -> list[int]:
        """
        0-based atom indices present in topology
        """
    def get_group(self, idx: int) -> GroupRec:
        """
        Group record at index (0..n_groups-1)
        """
    def get_ring(self, idx: int) -> RingRec:
        """
        Ring record at index (0..n_rings-1)
        """
    def has_computed(self, comp: TopologyComputers) -> bool:
        """
        Whether a stage completed successfully
        """
    def molecule(self) -> rdkit.RWMol:
        """
        Underlying RDKit molecule (borrowed reference)
        """
    def residue_index_of_atom(self, atom_idx: int) -> int:
        """
        Return residue index for an atom index
        """
    def residue_of_atom(self, atom_idx: int) -> Residue:
        """
        Return residue object for an atom index
        """
    def resolve_atom(self, id: ...) -> AtomRec:
        """
        Resolve an Atom EntityID to AtomRec (reference, lifetime tied to Topology)
        """
    def resolve_group(self, id: ...) -> GroupRec:
        """
        Resolve a Group EntityID to GroupRec (reference, lifetime tied to Topology)
        """
    def resolve_ring(self, id: ...) -> RingRec:
        """
        Resolve a Ring EntityID to RingRec (reference, lifetime tied to Topology)
        """
    def set_compute_nonstandard_bonds(self, compute: bool) -> None:
        """
        Include metal/coordination bonds if True
        """
    def set_cutoff(self, cutoff: float) -> None:
        """
        Neighbor cutoff used by bond perception (A)
        """
    @property
    def atom_records(self) -> list[AtomRec]:
        """
        Per-atom typing records (size N_atoms)
        """
    @property
    def atom_to_residue_indices(self) -> list[int]:
        """
        Atom-aligned residue index mapping (-1 for unmapped atoms)
        """
    @property
    def groups(self) -> list[GroupView]:
        """
        Detected functional groups as dynamic views
        """
    @property
    def residues(self) -> Residues:
        """
        Get residues container (kept alive by the Topology)
        """
    @property
    def rings(self) -> list[RingView]:
        """
        Detected rings as dynamic views
        """
class TopologyBuildingOptions:
    """
    Options controlling topology construction.
    """
    def __init__(self) -> None:
        """
        Create default options
        """
    @property
    def atom_typing_method(self) -> AtomTypingMethod:
        """
        Backend used for atom typing
        """
    @atom_typing_method.setter
    def atom_typing_method(self, arg0: AtomTypingMethod) -> None:
        ...
    @property
    def auto_heal(self) -> bool:
        """
        Resolve required dependencies automatically during execution
        """
    @auto_heal.setter
    def auto_heal(self, arg0: bool) -> None:
        ...
    @property
    def compute_nonstandard_bonds(self) -> bool:
        """
        Include non-standard/metal coordination where applicable
        """
    @compute_nonstandard_bonds.setter
    def compute_nonstandard_bonds(self, arg0: bool) -> None:
        ...
    @property
    def cutoff(self) -> float:
        """
        Neighbor-search cutoff used in bond perception and neighbors (A)
        """
    @cutoff.setter
    def cutoff(self, arg0: float) -> None:
        ...
class TopologyComputers:
    """
    Bitmask flags representing topology stages. Combine with | and &.
    """
    All: typing.ClassVar[TopologyComputers]  # value = <TopologyComputers.All: 4294967295>
    AtomTyping: typing.ClassVar[TopologyComputers]  # value = <TopologyComputers.AtomTyping: 32>
    Basic: typing.ClassVar[TopologyComputers]  # value = <TopologyComputers.Basic: 7>
    Bonds: typing.ClassVar[TopologyComputers]  # value = <TopologyComputers.Bonds: 2>
    Complete: typing.ClassVar[TopologyComputers]  # value = <TopologyComputers.Complete: 63>
    Extended: typing.ClassVar[TopologyComputers]  # value = <TopologyComputers.Extended: 31>
    Neighbors: typing.ClassVar[TopologyComputers]  # value = <TopologyComputers.Neighbors: 1>
    NoComp: typing.ClassVar[TopologyComputers]  # value = <TopologyComputers.NoComp: 0>
    NonStandardBonds: typing.ClassVar[TopologyComputers]  # value = <TopologyComputers.NonStandardBonds: 4>
    Residues: typing.ClassVar[TopologyComputers]  # value = <TopologyComputers.Residues: 8>
    Rings: typing.ClassVar[TopologyComputers]  # value = <TopologyComputers.Rings: 16>
    Standard: typing.ClassVar[TopologyComputers]  # value = <TopologyComputers.Standard: 15>
    __members__: typing.ClassVar[dict[str, TopologyComputers]]  # value = {'NoComp': <TopologyComputers.NoComp: 0>, 'Neighbors': <TopologyComputers.Neighbors: 1>, 'Bonds': <TopologyComputers.Bonds: 2>, 'NonStandardBonds': <TopologyComputers.NonStandardBonds: 4>, 'Residues': <TopologyComputers.Residues: 8>, 'Rings': <TopologyComputers.Rings: 16>, 'AtomTyping': <TopologyComputers.AtomTyping: 32>, 'Basic': <TopologyComputers.Basic: 7>, 'Standard': <TopologyComputers.Standard: 15>, 'Extended': <TopologyComputers.Extended: 31>, 'Complete': <TopologyComputers.Complete: 63>, 'All': <TopologyComputers.All: 4294967295>}
    def __and__(self, arg0: TopologyComputers) -> TopologyComputers:
        """
        Bitwise AND: intersect flags
        """
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __or__(self, arg0: TopologyComputers) -> TopologyComputers:
        """
        Bitwise OR: combine flags
        """
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
def compute_angles(topology: Topology, ring_indices: list[int], points: list[list[float]]) -> list[float]:
    ...
def covalent_radius(atomic_number: int) -> float:
    """
    Covalent radius (A) for the element with the given atomic number.
    """
def decode_contacts_batch_parallel(payloads: list, columnar: bool = True) -> list:
    ...
def decode_contacts_binary(payload: bytes) -> dict:
    ...
def decode_contacts_binary_columnar(payload: bytes) -> dict:
    ...
def decode_contacts_binary_direct(payload: bytes) -> dict:
    ...
def factorize(arg0: list[str], arg1: list[int], arg2: list[str]) -> tuple:
    """
    Factorize
    """
@typing.overload
def find_contacts(topology: Topology, kind_a: Kind, pred_a: typing.Callable, kind_b: Kind, pred_b: typing.Callable, opts: SearchOptions = ...) -> ContactSet:
    ...
@typing.overload
def find_contacts(topology: Topology, kind_a: Kind, pred_a: typing.Callable, kind_b: Kind, pred_b: typing.Callable, tester: typing.Callable, opts: SearchOptions = ...) -> ContactSet:
    ...
@typing.overload
def find_contacts(topology: Topology, kind: Kind, pred: typing.Callable, opts: SearchOptions = ...) -> ContactSet:
    ...
@typing.overload
def find_contacts(topology: Topology, kind: Kind, pred: typing.Callable, tester: typing.Callable, opts: SearchOptions = ...) -> ContactSet:
    ...
def vdw_radius(atomic_number: int) -> float:
    """
    van der Waals radius (A) for the element with the given atomic number.
    """
__version__: str = '2.0.3'
__version_info__: tuple = (2, 0, 3, '')
