## Selection and Filtering API
As mentioned, contacts are simply neighbor atoms that satisfy some additional criteria. In this section, we will discuss how to use Lahuta to extract contacts from neighbor atoms. We will also discuss how to customize the criteria that are used to extract contacts. It's best to start with an example:

???+ example "Example - Applying Filters to Neihboring Atoms"
    _Applying filters to neighboring atoms._
    ```py
    from lahuta import Luni

    # Load the structure and compute neighbors
    luni = Luni("path/to/file.pdb")
    ns = luni.compute_neighbors()
    print (ns.pairs.shape) # (1)!
    #> (3024, 2)

    # Apply a simple type filter
    aromatic1_ns = ns.type_filter("aromatic", 1) # (2)!
    print (aromatic1_ns.pairs.shape) # (3)!
    #> (287, 2)

    ```
    
    1. After computing neighbors, we have 3024 pairs of neighboring atoms.
    2. We apply a type filter to the neighbors. This filter will only keep neighbors where the first atom is aromatic. 
    3. After applying the filter, we are left with 287 pairs of aromatic neighbors.

!!! tip "Keep in mind"
    `pairs` and `distances` are uniquely sorted such that the first atom in a pair is the one with a lower index. This means that if we have a pair of atoms with indices `(i, j)`, then `i < j`. This way, we can uniquely specify the first and second atom in a pair. `1` is used to specify the first atom and `2` is used to specify the second atom. This is true for all methods that apply filters to `NeighborPairs` objects.

!!! tip "Also note"
    The `type_filter` method returns a new `NeighborPairs` object. The original `NeighborPairs` object is not modified. This Fluent API design pattern is used throughout Lahuta. It allows you to chain multiple filters together and apply them in a single step. For example, we can apply two filters in a single step as follows:

    ```py  
    aromatic_neighbors = ns.type_filter("aromatic", 1).type_filter("aromatic", 2)
    ```

    This will only keep neighbors where both atoms are aromatic.

This selection syntax is quite powerful and allows you to capture complicated relationships between atoms in an intuitive way. It is also very flexible and extensible. To further demonstrate this, we provide the full code for how aromatic contacts are computed in Lahuta:

???+ example "Example - Aromatic Contacts"
    _Computing aromatic contacts in Lahuta._
    ```py

    aromatic_ns = ns.type_filter("aromatic", 1).type_filter("aromatic", 2).distance_filter(4.0) # (1)!
    print (aromatic_ns.pairs.shape)

    ```

    1. This is all that is required to compute aromatic contacts!

Here is another example that demonstrates how to compute ionic contacts:

???+ example "Example - Ionic Contacts"
    _Computing ionic contacts in Lahuta._
    ```py

    contacts_atom12 = ns.type_filter("pos_ionisable", 1).type_filter("neg_ionisable", 2).distance_filter(4.0) # (1)!
    contacts_atom21 = ns.type_filter("neg_ionisable", 1).type_filter("pos_ionisable", 2).distance_filter(4.0) # (2)!

    contacts = contacts_atom12 + contacts_atom21 # (3)!

    ```

    1. This will only keep neighbors where the first atom is positively charged and the second atom is negatively charged.
    2. This will only keep neighbors where the first atom is negatively charged and the second atom is positively charged.
    3. This will combine the two `NeighborPairs` objects into a single `NeighborPairs` object.

!!! tip "Learn more"
    See the API documentation on [contacts](api/contacts.md) for more information.

I hope that these examples show how easy and intuitive it is to use Lahuta to extract contacts from neighboring atoms. I hope you also see how `NeighborPairs` objects can be combined to create more complex `NeighborPairs` objects. This is the core idea behind Lahuta's selection and filtering API. It is designed to be intuitive, flexible, and extensible.

## Computed Atom Types

Below are the atom types that are computed by Lahuta:

- `aromatic`: Aromatic atoms
- `pos_ionisable`: Positively charged atoms
- `neg_ionisable`: Negatively charged atoms
- `hbond_acceptor`: Hydrogen bond acceptor atoms
- `hbond_donor`: Hydrogen bond donor atoms
- `weak_hbond_acceptor`: Weak hydrogen bond acceptor atoms
- `weak_hbond_donor`: Weak hydrogen bond donor atoms
- `xbond_acceptor`: Halogen bond acceptor atoms
- `xbond_donor`:  Halogen bond donor atoms
- `hydrophobic`: Hydrophobic atoms
- `carbonyl_oxygen`: Carbonyl oxygen atoms
- `carbonyl_carbon`: Carbonyl carbon atoms

## Types of Filters

`NeighborPairs` objects can be filtered using the following filters (methods):

### Type Filter: `NeighborPairs.type_filter(...)`
:::lahuta.core.neighbors.NeighborPairs.type_filter
### Distance Filter: `NeighborPairs.distance_filter(...)`
:::lahuta.core.neighbors.NeighborPairs.distance_filter
### Index Filter: `NeighborPairs.index_filter(...)`
:::lahuta.core.neighbors.NeighborPairs.index_filter
### Numeric Filter `NeighborPairs.numeric_filter(...)`
:::lahuta.core.neighbors.NeighborPairs.numeric_filter
### Radius Filter: `NeighborPairs.radius_filter(...)`
:::lahuta.core.neighbors.NeighborPairs.radius_filter

There are also hbond-specific filters:

- `hbond_distance_filter`: Filters based on the distance between the hbonded atoms
- `hbond_angle_filter`: Filters based on the angle between the hydrogen bonded atoms

The last two filters operate on three atoms at a time and Lahuta implements vectorized versions of these filters. This means that they are very fast and efficient, but also that the code is a bit more complicated. For this reason, we will not discuss them here. Please see the API documentation for more information.

!!! tip "Learn more"
    See the API documentation on [NeighborPairs](api/neighborpairs.md) for more information.