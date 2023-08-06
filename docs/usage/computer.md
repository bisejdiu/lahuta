# A Uniform Interface for Computing Contacts

Both the Functional and Object-Oriented API described earlier offer a lot of flexibility, but they are also verbose and very expressive. Lahuta implements two different interfaces that are easier to use. The first is for structural data and the second is for MD simulation data. 

## Structural Data

We can use the `LahutaContacts` class to pool and compute contacts for structural data. With this class, you do not need to worry about performance. 

???+ example "Example - Computing Contacts for Structural Data"
    ```py
    from lahuta.contacts.computer import LahutaContacts

    luni = Luni("path/to/file.pdb")
    ns = luni.compute_neighbors()

    contacts = LahutaContacts()
    contacts.compute(ns)
    contacts.results # (1)!
    ```

    1. The `results` attribute contains a dictionary of all computed contacts. The keys are the contact types and the values are the computed contacts.

The following example provides a more detailed example of how to use the `LahutaContacts` class:

???+ example "Example - More involved Example"
    ```py

    contacts = LahutaContacts(contact_type='atom-plane') # (1)!
    contacts.compute(ns)

    contacts = LahutaContacts(contact_type=None) # (2)!
    custom_contact_func = ...
    contacts.register(custom_contact_func) # (3)!
    contacts.register([F.hydrophobic_neighbors, F.vdw_neighbors]) # (4)!
    
    contacts.compute(ns)
    contacts.results # (5)!
    ```

    1. `contact_type` can be any of None, `all`, `atom-atom`, `atom-plane`, or `plane-plane`. The default is `all`.
    2. Initialize with no contact type.
    3. `custom_contact_func` is our own custom function. We register it with the `LahutaContacts` class.
    4. We register two functions from the `F` module.
    5. `results` will only contain the results of the registered functions.

## MD Simulation Data

Computing contacts for MD simulation data is slightly different because we need to compute contacts for each frame. Lahuta has a special class, `LahutaTrajectoryContacts` that handles this for us. It will iterate over all frames in the trajectory and compute contacts for each frame. It can also do this in parallel, as the following example shows:

???+ example "Example - Computing Contacts for MD Simulation Data"
    ```py
    import MDAnalysis as mda
    from lahuta.contacts.computer import LahutaContacts, LahutaTrajectoryContacts

    universe = mda.Universe("path/to/file.pdb", "path/to/file.xtc")
    luni = Luni(universe.atoms)
    ns = luni.compute_neighbors()

    contacts = LahutaContacts(contact_type='atom-atom') # (1)!
    trj_contacts = LahutaTrajectoryContacts(res_dif=2, radius=5) # (2)!
    trj_contacts.compute(luni, contacts, n_jobs=4) # (3)!
    trj_contacts.results
    ```

    1. We initialize the `LahutaContacts` class with the `contact_type` argument set to `atom-atom`. This means that we will only compute atom-atom contacts.
    2. We initialize the `LahutaTrajectoryContacts` class with the `res_dif` and `radius` arguments set to 2 and 5, respectively. This means that we will only compute contacts between residues that are at most 2 residues apart and that are within 5 Å of each other.
    3. We compute contacts for each frame in the trajectory. We use 4 processes to do this in parallel.
    4. The `results` attribute contains a dictionary of all computed `atom-atom` contacts. It stores results for each frame in the trajectory. 

