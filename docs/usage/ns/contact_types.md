# Available Contact Types

Lahuta comes with several built-in contact types. The majority of contacts are atom-atom contacts, but there are also atom-plane and plane-plane contacts. We provide two different API's to access built-in contact types: functional and object-oriented. The functional API is more concise and easier to use, but the object-oriented API is more flexible and allows for more customization.

## Functional API

It's best to show how the functional API works with an example. The following code snippet shows how to extract all atom-atom contacts from a structure:

???+ example "Example - Getting Contacts using the Functional API"
    ```py

    from lahuta import Luni
    from lahuta.contacts import F # (1)!

    luni = Luni("path/to/file.pdb")
    ns = luni.compute_neighbors()

    # several atom-atom contacts
    aromatic = F.aromatic_neighbors(ns) # (2)!
    ionic = F.ionic_neighbors(ns)
    hydrophobic = F.hydrophobic_neighbors(ns)
    vdw = F.vdw_neighbors(ns)
    hbond = F.hbond_neighbors(ns)

    # atom-plane contacts
    carbonyl_pi = F.carbon_pi(ns)
    cation_pi = F.cation_pi(ns)
    sulphur_pi = F.sulphur_pi(ns)
    donor_pi = F.donor_pi(ns)

    # plane-plane contacts
    plane_plane = F.plane_plane_neighbors(ns)

    ```

    1. The `F` object provides access to all built-in contact types.
    2. May take optional arguments. For example, `distance` is an optional argument for `aromatic_neighbors`. If not provided, we read the default distance from `lahuta.config.defaults` module. 

The cool thing about the functional API is the uniformity of the API. All contact types are accessed using the same API. 

Note that atom-plane contacts are a bit special. This is reflected in an important performance consideration when using the Functional API for atom-plane contacts. Specifically, atom-plane contacts require the computation of several properties about the system, and the functional API forces the re-computation of these properties for each contact type. It is possible to circumvent this by caching the results of the computation. 

??? example "Example - Caching Results"
    ```py
    
    ns = ...

    # caching atom-plane contacts
    carbonyl_pi = F.carbon_pi(ns, cache=True) # (1)!
    cation_pi = F.cation_pi(ns, cache=True) # (2)!
    sulphur_pi = F.sulphur_pi(ns, cache=True)
    donor_pi = F.donor_pi(ns, cache=True)
    ```

    1. General atom-plane data will be computed once and cached for all atom-plane contacts.
    2. The cached data will be used. No re-computation will be done.

!!! warning Warning
    The `atom-plane` functional API is not thread-safe. This means that you cannot use it in a multi-threaded environment. If you need to use Lahuta in a multi-threaded environment, then you need to use the object-oriented API for atom-plane contacts.

!!! note
    Caching is not without trade-offs. It requires more memory and consumes time to store and retrieve the cached data. For this reason, caching is disabled by default, and you need to test your system to see if caching is beneficial or not.

## Object-Oriented API

The object-oriented API is more flexible and allows for more customization. It is also more verbose and requires more code to use. The following code snippet shows how to extract all atom-atom contacts from a structure:

???+ example "Example - Getting Contacts using the Object-Oriented API"
    ```py

    from lahuta import Luni
    from lahuta.contacts import ( # (1)!
        AromaticContacts, 
        IonicContacts, 
        HydrophobicContacts, 
        HBondContacts
    )

    luni = Luni("path/to/file.pdb")
    ns = luni.compute_neighbors()

    # aromatic contacts
    ac = AromaticContacts(ns)
    print (ac.distance) # (2)!
    ac.compute()
    ac.results

    # ionic contacts
    ic = IonicContacts(ns)
    ic.compute()
    ic.results

    # hydrophobic contacts
    hc = HydrophobicContacts(ns)
    hc.compute()
    hc.results

    # hydrogen bond contacts
    hbc = HBondContacts(ns)
    hbc.compute()
    hbc.results

    ```

    1. Each contact type has its own class. 
    2. The `distance` class attribute is used to access the distance cutoff.

With the object-oriented API, atom-plane contacts avoid the performance issues of the functional API:

???+ example "Example - Object-Oriented API for Atom-Plane Contacts"
    ```py

    from lahuta import Luni
    from lahuta.contacts import AtomPlaneContacts, PlanePlaneContacts

    luni = Luni("path/to/file.pdb")
    ns = luni.compute_neighbors()

    # atom-plane contacts
    apc = AtomPlaneContacts(ns)

    # plane-plane contacts
    ppc = PlanePlaneContacts(ns)

    # carbon-pi contacts
    apc.carbon_pi() # (1)!
    apc.cation_pi()

    # plane-plane contacts
    results = ppc.compute() # (2)!

    ```

    1. The AtomPlaneContacts instance holds the data for all atom-plane contacts.
    2. The `compute` method for PlanePlaneContacts returns the results.

## Custom Contact Types

The Functional API simply uses the selection/filtering language. As such, it is trivial to add new contact types. We highly recommend you check the API for the [functional API]() to see how to add new contact types.

To add new contact types using the object-oriented API, you can inherit from the `ContactAnalysis` base class and then implement one of two methods:

1. `compute` - This method is used to compute the contacts, similar to the functional API. It is fast because it is vectorized.
2. `compute_elementwise` - This method is used to compute the contacts element-wise. It is slower than `compute` because it iterates over all pairs of neighbor atoms.

The following code snippet shows how to implement a new contact type using the object-oriented API:

???+ example "Example - Implementing a New Contact Type"
    ```py

    from lahuta import Luni
    from lahuta.contacts import ContactAnalysis

    luni = Luni("path/to/file.pdb")
    ns = luni.compute_neighbors()

    class MyContact(ContactAnalysis): # (1)!
        def compute(self): -> "NeighborPairs" # (2)!
            # compute must return a NeighborPairs object

            # For example, return stored neighbors without any filtering
            return self.ns

    # my contact type
    mc = MyContact(ns) # (3)!
    mc.results

    ```

    1. Inherit from the `ContactAnalysis` base class.
    2. Implement the `compute` method.
    3. The `compute` method is automatically invoked when the `MyContact` object is created.    

???+ example "Example - New Contact Type using Element-Wise Computation"
    ```py

    from lahuta import Luni
    from lahuta.contacts import ContactAnalysis

    luni = Luni("path/to/file.pdb")
    ns = luni.compute_neighbors()

    class MyContact(ContactAnalysis): # (1)!
        def compute_elementwise(self, atom1, atom2, distance) -> Any: # (2)!
            return (atom1, atom2, distance)

    # my contact type
    mc = MyContact(ns) # (3)!
    mc.results [:3]
    #    [(<Atom 2: CA of type C of resname ALA, resid -3 and segid 1>,
    #    <Atom 170: NZ of type N of resname LYS, resid 9 and segid 1>,
    #    4.691295052545275),
    #    (<Atom 3: C of type C of resname ALA, resid -3 and segid 1>,
    #    <Atom 25: N of type N of resname LEU, resid -1 and segid 1>,
    #    4.526018623575078),
    #    (<Atom 4: O of type O of resname ALA, resid -3 and segid 1>,
    #    <Atom 25: N of type N of resname LEU, resid -1 and segid 1>,
    #    4.711188940149044)]

    ```

    1. Inherit from the `ContactAnalysis` base class.
    2. Implement the `compute_elementwise` method. See the `Note` below for more information.
    3. The `compute_elementwise` method is automatically invoked when the `MyContact` object is created.

!!! note
    `compute_elementwise` takes three arguments: `atom1`, `atom2`, and `distance`. The `atom1` and `atom2` arguments are `MDAnalysis.core.groups.Atom` objects representing the neighboring atoms. The `distance` argument is a `float` value representing the distance between the two atoms. The method must return a value that will be stored in the `results` attribute, as such it can be anything. In the example above we return a tuple of the form `(atom1, atom2, distance)`. 
    