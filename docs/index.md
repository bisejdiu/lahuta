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

The following code snippet provides a brief overview of the library's capabilities. For a more detailed introduction, please refer to the [User Guide](user_guide.md).

```python
from lahuta import Luni

# Create a Luni object from a PDB file
luni = Luni.from_pdb("path/to/file.pdb")

# Create a Luni object from a PDBx/mmCIF file
luni = Luni.from_mmcif("path/to/file.cif")


```
