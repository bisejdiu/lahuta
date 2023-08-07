The entry point to Lahuta is provided via its `Luni` class from its `core` module. Think of this class as the main interface to Lahuta, similar to `MDAnalysis.Universe`. In fact the name `Luni` stands for `Lahuta Universe`. 

There are two ways to initialize the `Luni` class. The first way is to pass a path to a PDB or PDBx/mmCIF file to it. The second way is to pass an `MDAnalysis.AtomGroup` object to it. 

???+ example "Example - Initializing the `Luni` Class Using File Paths"
    ```py 
    from lahuta import Luni

    luni_pdb = Luni("path/to/file.pdb") # (1)!
    luni_cif = Luni("path/to/file.cif") # (2)!
    luni_pdb_gz = Luni("path/to/file.pdb.gz") # (3)!
    luni_cif_gz = Luni("path/to/file.cif.gz") # (4)!
    luni_gro = Luni("path/to/file.gro") # (5)!
    
    ```

    1. The `Luni` class is initialized by passing a path to a PDB file to it.
    2. The `Luni` class is initialized by passing a path to a PDBx/mmCIF file to it.
    3. The `Luni` class is initialized by passing a path to a compressed PDB or PDBx/mmCIF file to it.
    4. The `Luni` class is initialized by passing a path to a compressed PDB or PDBx/mmCIF file to it.
    5. The `Luni` class is initialized by passing a path to a GRO file to it.

???+ example "Example - Initializing the `Luni` Class Using `MDAnalysis.AtomGroup` Objects"
    ```py 
    import MDAnalysis as mda
    from lahuta import Luni

    u = mda.Universe("path/to/file.pdb") # (1)!
    u_traj = mda.Universe("path/to/file.pdb", "path/to/file.xtc") # (2)!

    luni = Luni(u.atoms) # (3)!
    luni_traj = Luni(u_traj.atoms) # (4)!
    ```

    1. The `MDAnalysis.Universe` class is used to load the PDB file.
    2. The `MDAnalysis.Universe` class is used to load the PDB file and the XTC trajectory.
    3. The `Luni` class is initialized by passing an `MDAnalysis.AtomGroup` object to it.
    4. The `Luni` class is initialized by passing an `MDAnalysis.AtomGroup` object to it.

When using Lahuta with MDAnalysis, the primary adjustment is to supply the `Luni` class with an `MDAnalysis.AtomGroup` object rather than a file path. When an `MDAnalysis.Universe` object has an associated trajectory, Lahuta recognizes and utilizes it. In fact, as we will discuss later, Lahuta provides specialized features for extracting contacts from MD trajectories.

In essence, PDB files as well as PDBx/mmCIF files are supported. Compressed files are also supported. Further, by supporting MDAnalysis, Lahuta supports all file formats supported by MDAnalysis. This includes PDB, GRO, XTC, TRR, DCD, and many more. See the [MDAnalysis documentation]() for more information.

## Accessing Information 

Upon initializing the `Luni` class, the contained information is accessible through the `arc` attribute, specifically its `atoms`, `residues`, and `chains` sub-attributes. These sub-attributes are instantiated from the `Atoms`, `Residues`, and `Chains` classes, respectively, that adeptly encapsulate the pertinent data. Integral to these classes is their foundation on NumPy's [Structured Arrays](https://numpy.org/doc/stable/user/basics.rec.html). Structured Arrays in NumPy are efficient because they allow for the flexible indexing of complex data types within a contiguous memory block, facilitating rapid access and operations on the data. This choice of data structure not only offers streamlined storage but also optimizes the retrieval and manipulation of the loaded biomolecular details.

??? example "Example - Accessing Loaded Information"
    ```py 
    from lahuta import Luni

    luni = Luni("path/to/file.pdb")

    # Accessing Atoms 
    atoms = luni.arc.atoms # (1)!
    ids = atoms.ids
    elements = atoms.elements
    names = atoms.names
    coords = atoms.coordinates

    # Accessing Residues
    residues = luni.arc.residues # (2)!
    resids = residues.resids
    resnames = residues.resnames

    # Accessing Chains
    chains = luni.arc.chains # (3)!
    ids = chains.ids 
    auths = chains.auths
    labels = chains.labels

    print (names)
    print (resids)
    print (labels)

    #> array(['N', 'C', 'C', ..., 'H', 'H', 'H'], dtype='<U10')
    #> array([-3, -3, -3, ..., 90, 90, 90])
    #> array(['A', 'A', 'A', ..., 'B', 'B', 'B'], dtype='<U10')
    ```

    1. The `atoms` attribute contains an `Atoms` object that encapsulates the atom information.
    2. The `residues` attribute contains a `Residues` object that encapsulates the residue information.
    3. The `chains` attribute contains a `Chains` object that encapsulates the chain information.

??? example "Example - Accessing Information using Indexing" 
    ```py 
    from lahuta import Luni

    luni = Luni("path/to/file.pdb")

    print (luni.arc.atoms.names[:3]) # (1)!
    #> array(['N', 'CA', 'C'], dtype='<U10')

    atom_ix0 = luni.arc[0] # (2)!
    print (atom_ix0)
    #> Atom(name=N, id=0, element=N, type=N, resname=ALA, resid=-3, chain_label=A, chain_id=1)
    print (atom_ix0.resname)
    #> ALA

    atoms_ix20_ix22 = luni.arc[20:22] # (3)!
    #> [
    #>      Atom(name=H, id=20, element=H, type=H, resname=ASP, resid=-2, chain_label=A, chain_id=1),
    #>      Atom(name=HA, id=21, element=H, type=H, resname=ASP, resid=-2, chain_label=A, chain_id=1)
    #> ]
    ```

    1. The `Atoms` object can be indexed using NumPy's indexing syntax.
    2. Indexing the `ARC` object returns an `Atom` object that encapsulates that atom's information.
    3. Indexing the `ARC` object with a slice returns a list of `Atom` objects that encapsulate the atoms' information.
## Third-Party Libraries

Lahuta relies specifically on three libraries: [MDAnalysis](), [Gemmi](), and [OpenBabel](). MDAnalysis is used for loading data from MD simulations. Gemmi is used for parsing PDBx/mmCIF files. OpenBabel is used for SMARTS pattern matching and perception of chemical properties (bonds, aromaticity, etc.). Lahuta has internal API that allows it to use these libraries very efficiently and initialize the required objects depending on the data source. For example, if the data is loaded from a PDB file, then the `Gemmi.Structure` object is used for parsing the file. The `MDAnalysis.Universe`, and `OpenBabel.OBMol` objects are created from the `Gemmi.Structure` object. 

???+ example "Example - Accessing the Underlying Objects"
    ```py 
    from lahuta import Luni

    luni = Luni("path/to/file.pdb")

    # Accessing Underlying Objects
    mda_ag = luni.to("mda") # (1)!
    ob_obmol = luni.to("mol") # (2)!
    ```

    1. The `to` method is used to convert the `Luni` object to an `MDAnalysis.AtomGroup` object.
    2. The `to` method is used to convert the `Luni` object to an `OpenBabel.OBMol` object.

!!! warning "Warning"
    Calling the `to` method with `"mol"` as an argument will trigger the perception of chemical properties. This is done for performance reasons, because OpenBabel does not support vectorized OBMol creation and we have to iterate over all atoms instead. To ensure we only do this once, all required chemical perceptions and SMARTS pattern matching are done the first time the `to` method is called with `"mol"` as an argument. The resulting `OpenBabel.OBMol` object is then cached and returned on subsequent calls.

## Chemical Perception and SMARTS Pattern Matching

Lahuta relies on OpenBabel for chemical perception and SMARTS pattern matching. Depending on the system and the type of hardware, this process may be the most time consuming step in the analysis pipeline. We have made a few attempts to speed this process up, but challenges using OpenBabel's API and inherent difficulties in parallelizing unpicklable objects have limited our success. Regardless, for the vast majority of structures and systems, this process is very fast and should not be a bottleneck. 

For MD simulations, which can get very large, and if the user wishes to compute contacts for all atoms without any pre-processing or filtering, and if this analyis is done on old hardware, then this process may take a few seconds. Even then, this process is only done once and **not** for all frames in the trajectory.

??? example "Example - Getting the system ready"
    ```py 
    from lahuta import Luni

    luni = Luni("path/to/file.pdb")

    # Getting the system ready (compute perception and SMARTS pattern matching)
    luni.ready() # (1)!

    # Computing Neighbors
    ns = luni.compute_neighbors()
    ```

    1. This method is auto-invoked by `compute_neighbors` or when the `to` method uses `"mol"` as an argument. While it's usually automatic, you can call it directly if desired to make the computation explicit.

