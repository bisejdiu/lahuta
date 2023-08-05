# Neighbors and Contacts in Lahuta

## The Concept of Neighbor Atoms  
In the terminology that we employ, a "neighbor" refers to an atom placed within a defined proximity to another atom. In Lahuta, a neighbor is determined by a user-defined distance criterion from another atom. If one atom is within this specified distance, it's considered a neighbor. Lahuta captures this information in its `NeighborPairs` class, which is the main concept in Lahuta as will be detailed here and in the following sections.

 
## Distance cutoffs and Residue Differences  
Atomic distances in Lahuta are calculated using the standard Euclidean norm. The second parameter, the residue difference, determines adjacency at the residue level. For instance, setting the residue difference to zero will result in the inclusions of contacts between same and adjacent residues which likely is not indicative of functional importance. Contiguous residues are by definition adjacent, and thus, their inclusion increases the level of noise in the data. 

## Rapid Neighbor Search  
Lahuta uses a highly optimized neighbor search algorithm provided by MDAnalysis. The implementation is in CPython and it can be used for both periodic and non-periodic systems. The algorithm is based on the [Cell Lists](https://en.wikipedia.org/wiki/Cell_lists) algorithm and it is highly optimized for speed and efficiency. It is also capable of handling large systems with millions of atoms. This is usually more favorable than using a KDTree-based algorithm. 

## Transition from Neighbors to Contacts  
In a sense, neighbor atoms can be considered to already nr in contact with each other. Usually, however, we are interested in specific types of contacts or contacts that satisfy certain criteria. These criteria can be geometric, chemical, or both. For example, we might be interested in contacts between ionic or aromatic atoms or contacts between atoms that are within a certain distance of each other and satisfy a certain angle constraint. We may also want to apply restrictions on the difference between the Van der Waals radii of the atoms (e.g. avoid contacts between atoms that are too close to each other - this is useful for avoiding clashes). The `NeighborPairs` class in Lahuta is designed to facilitate the transition from neighbors to contacts. It provides a set of built-in filters and criteria that can be used to filter the neighbors and extract contacts.
