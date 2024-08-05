#include <GraphMol/RWMol.h>
#include "nsgrid.hpp"


bool connectOBMol(RDKit::Atom *p, RDKit::Atom *q, double dist_sq, double tolerance);
RDKit::RWMol lahutaBondAssignment(RDKit::RWMol &mol, const NSResults &results, std::vector<int> &non_protein_indices);
