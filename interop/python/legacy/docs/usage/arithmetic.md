This section follows up directly on the previous section and builds on the concepts introduced there, so please make sure you have read the [Contacts](ns/neighbors.md) section before continuing. 

In this section, we will discuss how `NeighborPairs` supports operator overloading and how this enables the creation of complex contact types using simple arithmetic-like expressions. This becomes quite powerful when applied to evolutionary related proteins. 

## Arithmetic Operators

For two `NeighborPairs` objects `a` and `b`, the following arithmetic operators are supported:

| Operator | Example | Description | Analogous set operation |
| --- | --- | --- | --- |
| `+` | `a + b` | Union | `a.union(b)` | 
| `-` | `a - b` | Difference | `a.difference(b)` |
| `&` | `a & b` | Intersection | `a.intersection(b)` |
| `|` | `a | b` | Symmetric difference | `a.symmetric_difference(b)` |
| `^` | `a ^ b` | XOR | `a.symmetric_difference(b)` |
| `>` | `a > b` | Strict superset | `a.issuperset(b)` |
| `<` | `a < b` | Strict subset | `a.issubset(b)` |
| `>=` | `a >= b` | Superset | `a.issuperset(b)` |
| `<=` | `a <= b` | Subset | `a.issubset(b)` |
| `==` | `a == b` | Equality | `a.isequal(b)` |
| `!=` | `a != b` | Inequality | `not a.isequal(b)` |

The following example illustrates the use of these operators:

???+ example "Example - Arithmetic Operators"
    ```py 
    from lahuta import Luni
    from lahuta.contacts import F # (1)!

    luni = Luni("path/to/file.pdb")
    ns = luni.compute_neighbors()

    # Get aromatic neighbors
    aromatic = F.aromatic_neighbors(ns)

    # Get carbonyl neighbors
    carbonyl = F.carbonyl_neighbors(ns)

    # Get aromatic or carbonyl neighbors
    aromatic_or_carbonyl = aromatic + carbonyl

    # Get aromatic and carbonyl neighbors
    aromatic_and_carbonyl = aromatic & carbonyl

    # Get aromatic but not carbonyl neighbors
    aromatic_not_carbonyl = aromatic - carbonyl # (2)!

    # Get aromatic or carbonyl neighbors but not both 
    aromatic_xor_carbonyl = aromatic | carbonyl # (3)!
    ```

    1. `F` contains the functional forms of pre-defined contact types. See the [F documentation]() for more information.
    2. For mutually exclusive contacts, `a - b = a` and `b - a = b`.
    3. For mutually exclusive contacts, `a | b = a + b`.

!!! tip "Note"
    Arithmetic or set operations between `NeighborPairs` objects will always return a new `NeighborPairs` object. The original `NeighborPairs` objects are not modified. 

!!! warning "Warning"
    While such set operations are very powerful and will almost always be valid syntax for `NeighborPairs` objects, you still have to make sure that the biological meaning and interpretation of the results make sense. You have to think about the meaning of the operations and then interpret it in the context of atom neighbors and contacts. For example, `a - b` will always be a subset of `a`. And as we saw above, for mutually exclusive contact definitions, `a - b = a` and `b - a = b`. This means that for such cases, the `difference` operations is not very useful.

## Set Operations
In addition to the arithmetic operators, `NeighborPairs` also supports the following set operations:

| Method | Example | Description |
| --- | --- | --- |
| `isdjoint` | `a.isdisjoint(b)` | Returns `True` if `a` and `b` have no common elements |
| `issubset` | `a.issubset(b)` | Returns `True` if `a` is a subset of `b` |
| `issuperset` | `a.issuperset(b)` | Returns `True` if `a` is a superset of `b` |
| `isequal` | `a.isequal(b)` | Returns `True` if `a` and `b` are equal |
| `isunique` | `a.isunique()` | Returns `True` if all pairs in `a` are unique |
| `is_strict_subset` | `a.is_strict_subset(b)` | Returns `True` if `a` is a strict subset of `b` |
| `is_strict_superset` | `a.is_strict_superset(b)` | Returns `True` if `a` is a strict superset of `b` |

The following example illustrates how to use these methods:

??? example "Example - Set Operations"
    ```py 
    from lahuta import Luni
    from lahuta.contacts import F

    luni = Luni("path/to/file.pdb")
    ns = luni.compute_neighbors()

    # Get aromatic neighbors
    aromatic = F.aromatic_neighbors(ns)

    # Get carbonyl neighbors
    carbonyl = F.carbonyl_neighbors(ns)

    aromatic.isdisjoint(carbonyl) # (1)!
    #> True

    aromatic.issubset(carbonyl) # (2)!
    #> False

    aromatic.isunique() # (3)!
    #> True

    ```

    1. Aromatic and carbonyl neighbors are disjoint.
    2. Aromatic neighbors are not a subset of carbonyl neighbors.
    3. All pairs in aromatic neighbors are unique.

!!! tip "Note"
    Uniqueness should be a guarantee and you will never encounter duplicate entries. If you find that a `NeighborPairs` object contains duplicate pairs, please report this as a bug on the [issue tracker]().
## Evolutionary Related Proteins

TBW