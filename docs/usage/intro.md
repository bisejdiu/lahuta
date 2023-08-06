# Introduction 

<img src="../big-logo.svg" alt="Alternative text" width="650" height="200">

## Background and Motivation
Lahuta is a computational tool that builds on the foundation laid by CREDO, a well-regarded database for understanding complex protein interactions. CREDO categorizes chemical entities according to their interaction with amino acid residues. It calculates these entities' chemical properties and consolidates this data in an easily accessible and comprehensive database.

Specifically, Lahuta leverages the refined and advanced SMARTS-based atom-typing system that was initially introduced by CREDO and further refined by Arpeggio. SMARTS, an acronym for SMiles ARbitrary Target Specification, is a language developed for the accurate description of molecular structures and patterns. By utilizing SMARTS, Lahuta is equipped to identify and categorize atoms based on defined criteria. 

The concept of atom typing is central to Lahuta's functionality. This process involves assigning a specific 'atom type' to each atom in a molecule, a step that is essential for accurately identifying and analyzing interatomic interactions. Atom typing, as it is implemented in Lahuta, forgoes the need to hard-code interactions based on exclusively predefined atom types. Instead, it opts for a more flexible, inclusive, and accurate approach rooted in the CREDO atom-typing system. As a result, Lahuta can recognize and analyze complex atomic interactions in protein structures with increased accuracy and thoroughness.

Lahuta supports almost all interaction types initially recognized by Arpeggio. These include hydrogen contacts, carbonyl contacts, π-π interactions, π-cation interactions, cation-π interactions, interactions with metal ions, van der Waals interactions, and more. Such interactions frequently occur in protein structures, and their recognition and analysis by Lahuta enables an in-depth understanding of complex atomic interactions.

## Neighbors and Contacts

Lahuta's core functionality is based on the concept of neighbors. In Lahuta, a neighbor is defined as an atom that is within a specified distance of another atom. This distance cutoff is user-defined. Lahuta stores neighbor information in the `NeighborPairs` class. It is the central class in Lahuta and everything else either builds on it or uses it. The next section in this documentation, [Contacts](), discusses the `NeighborPairs` class in detail.

Here we only outline its importance and to introduce it to you. A contact is then simply a pair of atoms that satisfy the user-defined distance criterion and any (none or several) additional criteria. For example, contacts with aromatic residues require that in addition to the distance criterion, one of the atoms must be aromatic.

### Built-in Contact Types

By default, Lahuta supports several commonly occurring contact types across three different non-exclusive contact categories. We list them below, but please note that, as will be described in the next section, Lahuta is designed to be customizable and easily extensible. This means that you can easily add your own contact types using the selection language described in the [Selection Language]() section.

#### Atom-Atom Contacts

1. `hbond` - Hydrogen bond interactions
2. `weak_hbond` - Weak hydrogen bond interactions
3. `polar_hbond` - Polar hydrogen bond interactions
4. `weak_polar_hbond` - Weak polar hydrogen bond interactions
5. `ionic` - Ionic interactions
6. `hydrophobic` - Hydrophobic interactions
7. `aromatic` - Aromatic interactions
8. `carbonyl` - Carbonyl interactions
9. `covalent` - Covalent interactions
10. `metal` - Metal interactions

#### Atom-Plane Contacts
1. `carbonyl-π` - Carbonyl-π interactions
2. `cation-π` - Cation-π interactions
3. `sulphur-π` - Sulphur-π interactions
4. `donor-π` - Donor-π interactions

#### Plane-Plane Contacts
1. `π-π` - π-π interactions
2. `π-π_t-shaped` - T-shaped π-π interactions
3. `π-π_parallel_displaced` - Parallel-displaced π-π interactions
4. `π-π_stacked` - Stacked π-π interactions

While the above list is not exhaustive, it covers the most common contact types. Adding new contact types is easy and straightforward and is described later in this documentation.

## Avoiding Re-inventing The Wheel

Writing a parser for even a common file format such as PDB is non-trivial and time-consuming. It is also very error-prone due to the many edge cases that need to be considered (e.g. missing information, negative residue numbers, etc.). For this reason, Lahuta relies on excellent third-party libraries to handle the parsing of file formats. This allows us to focus on the core functionality of the libraray and ensure that this core functionality is as robust and efficient as possible and available to as many users as possible.

Lahuta supports both PDB and PDBx/mmCIF file formats. It also supports loading data from MDAnalysis, a popular Python library for the analysis of MD simulations. This means that you can use Lahuta to analyze MD simulations without having to convert them to PDB or PDBx/mmCIF files. See the [Loading Data]() section for more information.


