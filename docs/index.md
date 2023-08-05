# Lahuta

{{ version }}.

Lahuta is a sophisticated and high-performance Python library designed for the analysis of non-covalent interactions within biomolecular systems. Its core functionality encompasses a wide range of applications, focusing on both speed and efficiency.

The library is architected with a user-friendly approach, promoting both intuitive usage and ease of extension. Lahuta also strives to provide comprehensive support for the most common file formats utilized in biomolecular simulations. 

Its development has been influenced by [arpeggio](https://github.com/harryjubb/arpeggio), which we consider as its predecessor. 


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
    See the [documentation on supported file formats](usage/loading.md#supported-file-formats) for more information.

{% raw %}
## Performance {#performance}
{% endraw %}

Lahuta is designed to be fast and efficient. Two overarching design principles are at the core of its architecture:

1. **Vectorization**: Lahuta relies heavily on vectorized operations to perform numerical computations. This is achieved by using NumPy arrays and SciPy sparse matrices.
2. **Ease of Use & Intuitive Usage**: Written in Python, Lahuta is not only easy to use but also easy to extend. Its high-level API is designed to be intuitive and user-friendly.


