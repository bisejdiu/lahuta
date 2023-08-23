# Lahuta
[![Github-CI][github-ci]][github-link]
[![Code style: black][black-badge]][black-link]
[![Conda][conda-badge]][conda-link]
[![Update Date][Anaconda-Update Badge]][update-link]
[![Downloads][Anaconda-Server Badge]][down-link]


<!-- [![Coverage Status][codecov-badge]][codecov-link] -->


<p align="center">
  <a href="#"><img src="https://github.com/bisejdiu/lahuta/blob/dev/docs/usage/big-logo.svg" alt="Lahuta"></a>
</p>

---

**Documentation**: <a href="https://bisejdiu.github.io/lahuta/" target="_blank">https://bisejdiu.github.io/lahuta/</a>

**Source Code**: <a href="https://github.com/bisejdiu/lahuta" target="_blank">https://github.com/bisejdiu/lahuta</a>

---

🧬 **Lahuta** is a Python package designed for chemically-aware contact analysis in biomolecular systems. 

Developed with a focus on performance and scalability, it enables the detailed analysis of large-scale molecular dynamics simulations, in addition to handling structural databases arising from current high-accuracy protein structure prediction pipelines such as AlphaFold, ESMFold, and others. 

**Key Features**:
- Comprehensive support for a variety of input data formats, including PDB, PDBx/mmCIF, GRO, XTC, DCD, and more.
- Ultra-responsive with sub- to low-millisecond computation times. (See benchmarks: [TBA]())
- A biologically-guided and data-driven API
- Capability to trace and juxtapose contacts between evolutionarily cognate proteins.
- First-class support for MD simulations
- Intuitive interface that is both user-friendly and extensible.
- Built-in exporters to common biological formats (including PyMOL, VMD). 
- Fully typed and comprehensively tested

## Motivation
Recent advances in high-accuracy structure prediction have resulted in an explosion in available structural data and dedicated databases to store and access the data. On the experimental side, improvements in sample preparation and, in large part thanks to the cryo-EM resolution revolution, the number of deposited structures has continued increasing along with the size of the deposited structures. Data files with hundreds of thousands of atoms are quite common. Parallel to all of this, there has been a great concerted effort to make Molecular Dynamics (MD) simulation data available and accessible. Such data is compositionally heterogeneous and has a time component along which the system evolves (trajectories). 

Gaining structural insight that is robust and consistent is very challenging. `Lahuta` aims to provide a uniform interface that works for both structural data and MD simulations, in addition to being performant and scalable. 

`Lahuta` is inspired by `arpeggio`, a widely-used tool for protein structure contact analysis, which remains in use even without active maintenance.

## Installation
We strongly recommend using `conda` for installation. Note that `lahuta` requires **Python 3.10** or higher to work.

```bash
conda create -n lahuta python=3.10 -y && conda activate lahuta
conda install -c bisejdiu lahuta
```

## Documentation
`Lahuta` provides extensive internal code documentation, an understandable API, and dedicated documentation pages. The [documentation](https://bisejdiu.github.io/lahuta/) provides a detailed usage guide with many examples, a few tutorials, and a detailed overview of the API. 
We highly recommend you go through the examples and usage guide in the documentation to learn how to use `Lahuta`. 

## Reporting Issues
If you encounter any issues, please report them in the [issues](https://github.com/bisejdiu/lahuta/issues) section of this project. 

[github-ci]: https://github.com/bisejdiu/lahuta/actions/workflows/ci.yml/badge.svg?branch=dev
[github-link]: https://github.com/bisejdiu/lahuta
[conda-badge]: https://anaconda.org/bisejdiu/lahuta/badges/version.svg
[conda-link]: https://anaconda.org/bisejdiu/lahuta
<!-- [codecov-badge]: https://codecov.io/gh/bisejdiu/lahuta/branch/main/graph/badge.svg -->
[black-badge]: https://img.shields.io/badge/code%20style-black-000000.svg
[black-link]: https://github.com/ambv/black
[Anaconda-Server Badge]: https://anaconda.org/bisejdiu/lahuta/badges/downloads.svg
[down-link]: https://anaconda.org/bisejdiu/lahuta
[Anaconda-Update Badge]: https://anaconda.org/bisejdiu/lahuta/badges/latest_release_date.svg
[update-link]: https://anaconda.org/bisejdiu/lahuta
