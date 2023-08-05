# Welcome to Lahuta

{{ version }}.

Lahuta is a sophisticated and high-performance Python library designed for the analysis of non-covalent interactions within biomolecular systems. Its core functionality encompasses a wide range of applications, focusing on both speed and efficiency.

The library is architected with a user-friendly approach, promoting both intuitive usage and ease of extension. Lahuta also strives to provide comprehensive support for the most common file formats utilized in biomolecular simulations (including Molecular Dynamics (MD) data). 

Its development has been influenced by [arpeggio](https://github.com/harryjubb/arpeggio), which we consider as its predecessor. 

{% raw %}
## Wetting Your Appetite {#quickstart}
{% endraw %}

The following example demonstrates how to compute non-covalent interactions using Lahuta.
For a more detailed introduction, please refer to the [User Guide](user_guide.md).

???+ example "Example - First dive into Lahuta"
    _The following code snippet shows a simple example of how to use Lahuta._
    ```py requires="3.9"
    from lahuta import Luni

    luni = Luni("path/to/file.pdb") # (1)!

    ns = luni.compute_neighbors(
        radius=4.5, # (2)!
        res_dif=2,  # (3)!
    )

    print (ns.pairs) # (4)!
    ```

    1. The `Luni` class is used as the entry point. You can supply a PDB, PDBx/mmCIF, or even a GRO file.
    2. The `radius` parameter specifies the cutoff distance for the neighbor search.
    3. The `res_dif` parameter specifies the minimum residue difference between two atoms (i.e. how far apart they are in the sequence).
    4. The `pairs` attribute contains an array of all pairs of atom indices that are within the specified cutoff distance and satisfy the residue difference constraint.

!!! tip "Learn more"
    See the documentation on [supported file formats](usage/loading.md#supported-file-formats) for more information.

{% raw %}
## Performance {#performance}
{% endraw %}

Lahuta is designed to be fast and efficient. Two overarching design principles are at the core of its architecture and implementation:

1. **Vectorization**: Lahuta relies heavily on vectorized operations to perform numerical computations. This is achieved by using NumPy arrays and SciPy sparse matrices.
2. **Ease of Use & Intuitive Usage**: Written in Python, Lahuta is not only easy to use but also easy to extend. Its high-level API is designed to be intuitive and user-friendly.

{% raw %}
## Scalability {#scalability}
{% endraw %}

Lahuta is designed to be scalable. It is capable of handling large biomolecular systems with millions of atoms, or directories containing thousands of files.
It is also designed to be able to handle large-scale MD simulations. This is important in the age of high-accuracty protein structure prediction, and parallel advances in workflows for quickly generating and simulating MD systems. 

Analyzing proteome-wide structure prediction databases on one side and large-scale MD simulations on the other side is a very challenging task. To ensure that Lahuta is capable of being used in such scenarios, we have designed it to be scalable and resource-efficient.  

{% raw %}
## Support for Molecular Dynamics (MD) Simulations {#md-support}
{% endraw %}

Lahuta provides support for MD simulations by relying on MDAnalysis for reading and writing structures and trajectories. The way this is done is that the entry point to Lahuta can be triggered by passing an `MDAnalysis.AtomGroup` object to the `Luni` class. This is demonstrated in the following example:

???+ example "Example - Using Lahuta with MDAnalysis"
    _The following code snippet shows how to use Lahuta with MDAnalysis._
    ```py requires="3.9"
    import MDAnalysis as mda
    from lahuta import Luni

    u = mda.Universe("path/to/file.pdb") # (1)!
    luni = Luni(u.atoms) # (2)!

    ns = luni.compute_neighbors(
        radius=4.5,
        res_dif=2,
    )

    print (ns.pairs)
    ```

    1. The `MDAnalysis.Universe` class is used to load the PDB file.
    2. The `Luni` class is initialized by passing an `MDAnalysis.AtomGroup` object to it.

As you can see, the only difference is that we pass an `MDAnalysis.AtomGroup` object to the `Luni` class instead of a file path. This is the only change that is required to use Lahuta with MDAnalysis. If a trajectory is loaded to the `MDAnalysis.Universe` object, then it will be also available and used by Lahuta. In fact, Lahuta provides dedicated support for extracting contacts from MD trajectories. 

!!! tip "Learn more"
    See the documentation on [working with MD Simulation data](usage/loading.md#md-contacts) for more information.

{% raw %}

## Command Line Interface (CLI) {#cli}
{% endraw %}

For many common tasks, a dedicated API, regardless of how intuitive it is, can be cumbersome to use. For this reason, Lahuta also comes with a command line interface (CLI) that provides a convenient way to perform common tasks. The following example demonstrates how to use the CLI to compute contacts:

???+ example "Example - Using the CLI"
    _The following code snippet shows how to use the CLI to compute contacts._
    ```bash
    $ lahuta contacts -h
    usage: lahuta contacts [-h] [-o OUTPUT] [-r RADIUS] [-d RES_DIF] [-v] input

    positional arguments:
    input                 Path to the input file.

    optional arguments:
    -h, --help            show this help message and exit
    -o OUTPUT, --output OUTPUT
                            Path to the output file.
    -r RADIUS, --radius RADIUS
                            Cutoff distance for the neighbor search.
    -d RES_DIF, --res-dif RES_DIF
                            Minimum residue difference between two atoms.
    -v, --verbose         Increase verbosity.
    ```

    ```bash
    $ lahuta contacts path/to/file.pdb -o contacts.csv -r 4.5 -d 2
    ```

!!! tip "Learn more"
    See the documentation on the [Command Line Interface (CLI) ](usage/cli.md) for more information.


{% raw %}
## Tracing Contacts Through Evolutionary Time {#evolutionary-time}
{% endraw %}

We can use Lahuta to compute contacts from related proteins and trace them through evolutionary time. Lahuta is capable of parsing MSA files and generate unique identifiers for each residue. It can then quickly map between real atom & residue indices and MSA-based indices. 

This opens the door to a whole class of applications with a wide range of use cases. We discuss this more in the [Tutorials](tutorials.md) and [Usage](usage.md) sections. We also have a [Blog Post]() that discusses this in more detail in the context of GPCRs.

{% raw %}
## Standing on the Shoulders of Giants {#dependencies}
{% endraw %}

Lahuta is written from the ground up in Python. Nevertheless, it relies on a few external libraries to offload tasks such as I/O and numerical operations. The following list provides a brief overview of the libraries that are used by Lahuta and their respective purposes:

### Bio- and Cheminformatics Libraries:
- [MDAnalysis](): Used for reading and writing biomolecular structures and trajectories, as well as facilitating rapid neighbor search.
- [Gemmi](): A specialized library for handling the reading and writing of structures in the PDBx/mmCIF format.
- [OpenBabel](): Essential for SMARTS pattern matching and advanced chemical perception.

### General-Purpose Libraries:
- [NumPy](): Used for fast numerical operations on arrays by relying extensively on vectorized operations.
- [SciPy](): We employ sparse matrices to efficiently store and quickly access atom indices. 


{% raw %}
## Active Development {#active-development}
{% endraw %}

Lahuta is under active development. We are constantly adding new features and improving existing ones. We are also working on improving the documentation and adding more tutorials. If you have any questions or suggestions, please feel free to [open an issue]() or [contact us](). 