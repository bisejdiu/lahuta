# Release Notes for Lahuta v0.9.1 (2023-08-25)
[GitHub release v0.9.1](https://github.com/bisejdiu/lahuta/releases/tag/v0.9.1)

**New Features**:
- Mapping of Indices Based on Multiple Sequence Alignment ([#41](https://github.com/bisejdiu/lahuta/pull/41)): Introducing a new mechanism to map neighbor information from residues in MSA. This initial support lays the foundation for further development.
- Basic Plotting Interface ([#45](https://github.com/bisejdiu/lahuta/pull/45)): A new plotting interface provides users with straightforward visualization tools.
- CLI Interface ([#47](https://github.com/bisejdiu/lahuta/pull/47)): Access and control Lahuta functionalities via a new command-line interface.

**Improvements**:
- Documentation Overhaul ([#39](https://github.com/bisejdiu/lahuta/pull/39)): Comprehensive updates and improvements to documentation enhance usability and understanding.
- Extending and Improving Tests ([#50](https://github.com/bisejdiu/lahuta/pull/50)): Robustness of the codebase has been increased with significant testing enhancements. 
- Refactoring & Python 3.10 Support ([#44](https://github.com/bisejdiu/lahuta/pull/44), [#56](https://github.com/bisejdiu/lahuta/pull/56)): Codebase refactoring for better performance and compatibility with Python 3.10 as the minimum supported version.

**Set Operations & Testing**:
- Set Operations Testing ([#43](https://github.com/bisejdiu/lahuta/pull/43)): Rigorous first-principles testing on set operations ensures reliability and consistency.
- Refactoring & Test Improvement ([#44](https://github.com/bisejdiu/lahuta/pull/44)): Enhancing the code quality and performance through systematic refactoring and improved testing.

**CI/CD Enhancements**:
- Fixed Building and Deployment of Documentation ([#51](https://github.com/bisejdiu/lahuta/pull/51)): Streamlined building process and deployment for accessible and up-to-date documentation.
- Linting & Type Checking ([#53](https://github.com/bisejdiu/lahuta/pull/53)): Code maintainability and correctness have been strengthened through rigorous linting and type checking.
- Build Testing Workflow ([#55](https://github.com/bisejdiu/lahuta/pull/55)): Continuous integration now includes build testing for regular validation.
- Publish Conda Package ([#57](https://github.com/bisejdiu/lahuta/pull/57)): Automated publishing of conda packages for every new release for seamless distribution.

Lahuta v0.9.1 is quite a significant update in both features and overall quality. 


# Release Notes for Lahuta v0.8.0 (2023-08-03)
[GitHub release v0.8.0](https://github.com/bisejdiu/lahuta/releases/tag/v0.8.0)

This release represents a significant move from a proof-of-concept to a fully-fledged library. Here's a brief overview of everything that's new in this release.

## Key Highlights

- **Functionality:** Lahuta now supports a wide range of features and contact types. 
- **Dual API Design:** Lahuta is versatile, with both functional and OOP APIs.
- **Precise Data Handling:** The Selection & Filtering Syntax allows for a much better handling of input data.
- **Fluent Design:** Makes it possible to do method chaining for a more intuitive and streamlined workflow.
- **Strict Typing:** Almost the entire codebase is now fully typed in the strictest mode.
- **Comprehensive Documentation:** Thorough documentation is now available for all modules, classes, and methods.
- **Testing:** Lahuta now has a comprehensive test suite that covers all major features and functionalities.
- **Trajectory Support:** Lahuta now supports the efficient handling of MD trajectories.
- **Performance:** Lahuta is now very fast and efficient, with a focus on memory usage and speed. 

## New Features & Major Changes

### CIF File Support [#14]
- **Added support for CIF file formats.** Now, users can work with CIF files, making the data input process more straightforward.
    - [Commit 604ddeb](https://github.com/bisejdiu/Lahuta/commit/604ddeb) introduced this feature.
    - Bugs related to this feature were addressed and support was extended in [Commit 65acc71](https://github.com/bisejdiu/Lahuta/commit/65acc71).

### Enhanced Molecular Dynamics (MD) Trajectory Support [#19]
- **Integrated MDAnalysis for Trajectory Handling.** This enables Lahuta to compute atom-atom contacts over simulation trajectories for a more dynamic view of interatomic interactions.
    - **Abstraction for Data Loading Improved:** Efficiency and versatility at their best. Check out the changes in [Commit 85e2d0f](https://github.com/bisejdiu/Lahuta/commit/85e2d0f).
    - **Optimizations in Atom Matching:** Methods for atom matchers and assigners are now optimized for faster computation. Technical details in [Commit eca4c37](https://github.com/bisejdiu/Lahuta/commit/eca4c37) and [Commit 949826e](https://github.com/bisejdiu/Lahuta/commit/949826e).
    - **Enhanced Index Handling:** Proper handling of indices for different atom selections now has better support. More in [Commit 8a77d77](https://github.com/bisejdiu/Lahuta/commit/8a77d77).
    - **Refined Trajectory Loading:** We now support the efficient loading and accessing of trajectories. More in [Commit 685ae8d](https://github.com/bisejdiu/Lahuta/commit/685ae8d).
    - **Improved Testing:** Our tests are now more aligned with the refactored codebase, ensuring that everything works just as expected. More in commits [9a1979a](https://github.com/bisejdiu/Lahuta/commit/9a1979a), [606a049](https://github.com/bisejdiu/Lahuta/commit/606a049), and [db668ae](https://github.com/bisejdiu/Lahuta/commit/db668ae).

### Updated Tests for Refactored Code [#21]
- **Tests Overhaul:** Aligned tests with the entirely refactored codebase, ensuring consistency and robustness. This addresses issue #20.

### Improved Covalent Contacts Computation [#23]
- **Covalent Contacts Enhancement:** Refined the computation of covalent contacts, which directly impacts issue #22. 

### New Contact Interface [#24]
- **Example Base Class for Contacts:** Added an example base class for contacts that serves as a prototype. Introduced in [Commit 3db583e](https://github.com/bisejdiu/Lahuta/commit/3db583e).
- **Optimized `mda.Universe` Creation:** Refactored the `mda.Universe` from scratch, adopting a faster and superior approach. Check out the changes in [Commit ff29ac2](https://github.com/bisejdiu/Lahuta/commit/ff29ac2).

### Implementation of ARC (Atoms, Residues, Chains) [#25]
- **Leveraging Structured Numpy Arrays:** Refactored the Atoms, Chains, and Residues (ARC) classes to use structured numpy arrays, enhancing code readability, reducing memory usage, and standardizing atomic data handling.
    - **ARC Classes Overview:** Atoms store atom information; Chains handle chain labels, authentication codes, and ids; Residues manage residue names and ids.
    - **Enhanced Indexing & Iteration:** Updated `__getitem__` and `__iter__` methods for improved efficiency and cleaner code.
    - **Improved Data Loading:** Updated methods `from_gemmi` and `from_mda` for creating structured arrays from their respective inputs.
    - **Streamlined Creation of Unique Structured Arrays:** Simplified the creation process of structured arrays for unique residue names, ids, and chain ids.
    - **Introduction of Atom Class:** A new Atom class represents individual atoms in the structure.
- **Notable Commits for ARC:** Initial ARC implementation with working Atoms class in [Commit 205b9ba](https://github.com/bisejdiu/Lahuta/commit/205b9ba), Chains and Residues classes using structured arrays in [Commit 32691f9](https://github.com/bisejdiu/Lahuta/commit/32691f9), renaming from CRA to ARC in [Commit 2295b77](https://github.com/bisejdiu/Lahuta/commit/2295b77), and ensuring consistent ARC behavior in [Commit 15eeb44](https://github.com/bisejdiu/Lahuta/commit/15eeb44).


### Enhanced API and Expanded Functionality for Interaction Calculations [#26]
#### **Introduction of Functional and Class-based Approaches:** 

The new release introduces two major paradigms for computing atom-atom interactions.
  - **Functional Approach:** Allows efficient calculation of specific contact types by passing relevant data to corresponding functions.
    - Sample usage:
      ```python
      universe = Universe(...)
      ns = universe.compute_neighbors()

      # functional
      from lahuta.contacts import F

      carb = F.carbonyl_neighbors(ns)
      aromatic = F.aromatic_neighbors(n)
      ...
      ```
  - **Class-Based Approach:** Predefined classes corresponding to various contact types are provided. Instantiating these classes with relevant data provides results through the `results` attribute.
    - Sample usage:
      ```python
      from lahuta.contacts import AromaticContacts, AtomPlaneContacts

      # class-based
      ac = AromaticContacts(ns)
      print(ac.results)
      print(AtomPlaneContacts(ns).donor_pi())
      ```

#### **API Consistency and Feature Expansion** 

Maintained uniformity in API usage across various contact types. Atom-plane and plane-plane computation logic have been significantly improved. 

  - **Major Commit Highlights:**
    - Fixed various bugs related to `mda AtomGroup`, `NeighborPairs`, and atom-plane contacts, as seen in commits [d42edda](https://github.com/bisejdiu/Lahuta/commit/d42edda), [fa0854d](https://github.com/bisejdiu/Lahuta/commit/fa0854d), and [10c2989](https://github.com/bisejdiu/Lahuta/commit/10c2989).
    - Enhancements to ARC Atoms to store positions/coordinates simplifying various processes [84335d7](https://github.com/bisejdiu/Lahuta/commit/84335d7).
    - Implementation of a dedicated module for writing results and integrating df and dict writers to `NeighborPairs` [62267f3](https://github.com/bisejdiu/Lahuta/commit/62267f3).
    - Refinement of class-based contact definitions and the introduction of new tests [50c028f](https://github.com/bisejdiu/Lahuta/commit/50c028f) and [6f15fa1](https://github.com/bisejdiu/Lahuta/commit/6f15fa1).
    - Overhaul of plane-plane interactions, offering enhanced functionality [c699850](https://github.com/bisejdiu/Lahuta/commit/c699850).
    - Addition of support for annotations to `NeighborPairs` and `DataFrameWriter`, improving data contextuality [34285ab](https://github.com/bisejdiu/Lahuta/commit/34285ab).


#### **Typing Support**
   - Added extensive type annotations to Lahuta.
   - Introduced types module which defines basic wrapper types.
   - Integrated tools like mypy, Pylance, and pyright for static type checking and linting.

#### **Extensive Documentation**
   - Module docstrings were added to various modules like config, contacts, types, array utils, arc.py, _loaders.py, assigners, matchers, universe class, writers, and neighbors.py.
   - Detailed documentation added for specific functions, methods, and classes.
   - Documentation-related configs updated (e.g., pylint and vscode settings).

#### **Code Cleanup & Refactoring**
   - Unnecessary files removed: old protocol, loaders files, main.py, types/base.py, and outdated utilities.
   - Code formatting adjustments: increased line length to 100, updated files to reflect this, and integrated the black formatter.
   - Refactoring and updating code: replaced lambda calls with `getattr` and dynamic method calls, refactored vdw radii and hbonded atom calculations.
   - Enhanced configurations: updates to lahuta.toml, pyright config, and pylint config.

#### Commits of Note
- **Code Refinement**: Shifted hbond_array creation responsibility to the NeighborPairs class and decoupled the neighbor class from the universe class.
- **Typing**: Added tentative typing annotation for NeighborPairs and fully typed modules achieved.
- **Code Cleanup**: Removed redundant/unnecessary files including the old protocol, loaders files, utils module files, and main.py.
- **Documentation**: Extensive documentation updates across various modules and functions.


### Testing & Bug Fixes
#### **Array Utilities Testing**
   - Rigorous testing of the functions within `array_utils.py`, which consists mainly of pure functions.
   - Utilization of randomly generated data for covering a vast range of test cases and edge conditions.

#### **NeighborPairs Testing & Functionality**
   - Validation of operations supported by NeighborPairs using tested functions from `array_utils.py`.
   - Integration of Python's built-in operator functionality to provide an intuitive user interface. Specific operators supported include:
     - `in` or `operator`: Tests if all pairs in one object are in another (issubset).
     - `+`: Union operation that combines pairs from two objects.
     - `-`: Difference operation returning pairs from the first object that aren't in the second.
     - `|` and `^`: Both represent the symmetric difference, returning unique pairs to each object.
     - `==`: Tests equality of two objects based on pairs and distances.
     - `&`: Intersection operation returning pairs common to both objects.

## Set Operations & Sparse Arrays
#### **Operator Decisions**
   - Ongoing decisions on the usage of `|` and `^` operators. One will represent symmetric difference, while the purpose of the other is still under consideration.


#### **Sparse Arrays Implementation**
   - Transitioned from using mapping dictionaries and NumPy arrays to sparse arrays for efficiency.
   - Implemented column-wise slicing, row-wise slicing, and DOK format for data population, ensuring reduced memory consumption.

#### **Trajectory Support**
   - Made robust improvements in trajectory support.
   - Added tests specifically for validating trajectory functionalities.

#### **Type Hints & Code Readability**
   - Updated type hints across various modules, bolstering code clarity and maintenance.

#### **OBMol Generation**
   - Refactoring to use `atom_id` rather than `idx` for generating OBMol, ensuring accurate atom identification.

#### **Vectorized Atom Assigner & Sparse Matrices**
   - Adaptation of vectorized atom assigner to be compatible with sparse matrices.

#### **Atom Types Correction**
   - Addressed issues of wrongly concatenated atom types.
   - Updated relevant types for both loaders and assigners.

#### **Code Cleanup**
   - Removed redundant testing code and incorporated incremental changes to improve overall code quality.

#### Other Noteworthy Changes
- Enabled `lahuta` to support various alternative residue names.
- Streamlined handling of indices ensuring compatibility across different libraries.
- Eliminated the necessity of providing `n_atoms` to the atom assigner.


### Trajectories: Testing & Bug Fixes
#### **Trajectory Compatibility**
   - Fixed issues with how trajectories are supported, specifically as they are read by MDAnalysis.

#### **Trajectory Testing**
   - Introduced tests aimed at validating universes with attached trajectories.
   - New data has been added specifically for trajectory testing, including both small and large MD files. The smaller datasets allow for faster test execution.

#### **System Dimension Check**
   - Before the computation of neighbors, the system now checks if dimensions are lacking. In such cases, a coordinate shift is executed to correct the problem.

#### **Conftest File Integration**
   - Introduced a `conftest` file to provide more granular control over testing procedures.
   - Implemented a flag `--large-files` for pytest which facilitates testing with sizable MD files.

#### Commit Details
- **6b88203**: Fixed contacts to be compatible with trajectory files.
- **ba4e5c4**: Replaced the older 'indices' term with 'ix'.
- **896e07e**: Added pertinent data specifically for trajectory testing.
- **4bd5f5f**: Included smaller systems for trajectory testing which hastens the testing process.
- **3a5aa6a**: Integrated code to validate trajectory support thoroughly.
- **de8f383**: Incorporated the `conftest` file to enhance test control.

