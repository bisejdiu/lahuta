
<!-- 1. Add chainID to labels perhaps? -->
<!-- 2. Add filter method to Luni  -->
<!-- 3. Add convenience methods to Luni to access system information -->
<!-- 4. Now that we have a filter method to Luni, we may want to add a method to add and/or substract Luni objects
   together. The problem here is that combining multiple Luni objects together may not provide any benefit 
   from just using the selection syntax to the filter. There is not really anything to gain from this. However,
   what we can do is provide a Luni instance as a parameter to `Luni.neighbors` calculation. The meaning of this
   is that the resulting neighbors will be filtered to only include the atoms with the provided Luni instance.
   That is, a neighbor pair of atoms must have at least one atom with the provided Luni instance. 
 -->

- [x] Luni now accepts Path objects. Update tests to reflect this (avoid casting to str). 
- [x] For custom DSSP selection, ensure we also take into account the chainID.
- [x] Update MSA mappings to include chainID information.
- [x] Tests for DSSP API 
- [ ] Add docs to the filter method 
- [x] Add tests for the filter method
- [x] How does dssp handle MD files?
  - It seems to only consider protein residues and ignore everything else. 
  - Test that the DSSP API works with MD files.

- [ ] Implement a `pairwise_distances` method for Luni.