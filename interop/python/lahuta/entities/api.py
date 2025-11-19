from __future__ import annotations

from typing import Callable, Iterable, Literal, TypeAlias, overload

from .._interaction_utils import normalize_interaction_selection
from ..lib import lahuta as lxx
from .selectors import Selector

# fmt: off
Tester:       TypeAlias = Callable[[int, int, float], lxx.InteractionType]
Provider:     TypeAlias = Literal["molstar", "arpeggio"]
SelectorLike: TypeAlias = Selector[lxx.AtomRec] | Selector[lxx.RingRec] | Selector[lxx.GroupRec]

DEFAULT_PROVIDER: Provider = "molstar"
SUPPORTED_PROVIDERS = frozenset(["molstar", "arpeggio", "getcontacts"])

def _create_search_options(distance_max: float | None = None, opts: lxx.SearchOptions | None = None) -> lxx.SearchOptions:
    """Create or modify search options.

    Args:
        distance_max: Maximum search distance in Angstroms
        opts: Existing search options to modify, or None to create new

    Returns:
        Configured search options object
    """
    search_opts = opts or lxx.SearchOptions()
    if distance_max is not None:
        if distance_max <= 0:
            raise ValueError(f"distance_max must be positive, got {distance_max}")
        search_opts.distance_max = float(distance_max)
    return search_opts

def _get_selector(selector: SelectorLike) -> Callable[..., bool]:
    """Extract selector function from selector, using a default true function if None.

    Args:
        selector: Entity selector with optional selector

    Returns:
        Callable selector function
    """
    if selector.selector is not None:
        return selector.selector

    def _always_true(_: SelectorLike) -> bool:
        return True

    return _always_true

@overload
def find_contacts(
    topology: lxx.Topology,
    a: SelectorLike,
    *,
    tester: Tester | None = ...,
    distance_max: float | None = ...,
    opts: lxx.SearchOptions | None = ...,
) -> lxx.ContactSet:
    """Find contacts within a single entity type."""
    ...


@overload
def find_contacts(
    topology: lxx.Topology,
    a: SelectorLike,
    b: SelectorLike,
    *,
    tester: Tester | None = ...,
    distance_max: float | None = ...,
    opts: lxx.SearchOptions | None = ...,
) -> lxx.ContactSet:
    """Find contacts between two different entity types."""
    ...


def find_contacts(
    topology: lxx.Topology,
    a: SelectorLike,
    b: SelectorLike | None = None,
    *,
    tester: Tester | None = None,
    distance_max: float | None = None,
    opts: lxx.SearchOptions | None = None,
) -> lxx.ContactSet:
    """Find contacts between molecular entities using flexible selectors.

    This function provides a unified interface for contact detection between
    different types of molecular entities (atoms, rings, functional groups).

    Args:
        topology: Molecular topology containing entity information
        a: Primary entity selector (atoms, rings, or groups)
        b: Secondary entity selector. If None, finds contacts within 'a'
        tester: Optional custom interaction tester function. Takes (index1, index2, distance^2)
               and returns InteractionType or bool
        distance_max: Maximum distance for contact detection (Angstroms)
        opts: Advanced search options. If provided, distance_max is ignored

    Returns:
        ContactSet containing all detected molecular contacts

    Raises:
        TypeError: If topology is not a valid topology type
        ValueError: If distance_max is not positive

    Examples:
        >>> # Find all atom-atom contacts within 5A
        >>> contacts = find_contacts(topology, atoms(), distance_max=5.0)

        >>> # Find hydrogen bonds between donors and acceptors
        >>> donors = atoms(lambda a: a.type.has(AtomType.HbondDonor))
        >>> acceptors = atoms(lambda a: a.type.has(AtomType.HbondAcceptor))
        >>> hbonds = find_contacts(topology, donors, acceptors, tester=hbond_tester, distance_max=3.5)
    """
    # Input validation and preparation
    search_opts = _create_search_options(distance_max, opts)

    # Extract selectors
    sel_a = _get_selector(a)

    # Handle single entity type (self-contact)
    if b is None:
        if tester is None:
            return lxx.find_contacts(topology, a.kind, sel_a, search_opts)
        return lxx.find_contacts(topology, a.kind, sel_a, tester, search_opts)

    # Handle dual entity types (cross-contact)
    sel_b = _get_selector(b)
    if tester is None:
        return lxx.find_contacts(topology, a.kind, sel_a, b.kind, sel_b, search_opts)
    return lxx.find_contacts(topology, a.kind, sel_a, b.kind, sel_b, tester, search_opts)


def compute_contacts(
    topology: lxx.Topology,
    provider: Provider = DEFAULT_PROVIDER,
    only: lxx.InteractionType | lxx.InteractionTypeSet | Iterable[lxx.InteractionType] | None = None,
) -> lxx.ContactSet:
    """Compute molecular contacts using predefined interactions types and providers.

    Args:
        topology: Molecular topology containing structural information
        provider: Contact detection algorithm ('molstar' or 'arpeggio')
        only: Optional filter to detect only specific interaction types

    Returns:
        ContactSet containing all detected molecular interactions

    Raises:
        TypeError: If topology is not a valid topology type
        ValueError: If provider is not supported

    Examples:
        >>> all_contacts = compute_contacts(topology)
        >>> hbonds = compute_contacts(topology, "arpeggio", lxx.InteractionType.HydrogenBond)
    """
    if provider not in SUPPORTED_PROVIDERS:
        raise ValueError(f"Unsupported provider '{provider}'. Supported providers: {', '.join(SUPPORTED_PROVIDERS)}")

    engine = lxx.MolStarContactsEngine() if provider == "molstar" else lxx.ArpeggioContactsEngine()

    # Compute contacts with optional filtering
    if only is None:
        return engine.compute(topology)

    interaction_filter = normalize_interaction_selection(only, arg_name="only")
    if interaction_filter.is_all():
        return engine.compute(topology)
    return engine.compute(topology, interaction_filter)

__all__ = [
    "find_contacts",
    "compute_contacts",
    "SUPPORTED_PROVIDERS",
]
