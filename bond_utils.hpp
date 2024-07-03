#include "GraphMol/Atom.h"
#include "GraphMol/RWMol.h"

double computeLengthSq(const RDKit::Conformer &conf, const RDKit::Bond *bond);

unsigned int GetExplicitValenceFromAtomBonds(const RDKit::RWMol &mol,
                                             const RDKit::Atom *atom);

bool isAngleLessThan45Degrees(const RDKit::RWMol &mol,
                              const RDKit::Conformer &conf,
                              const RDKit::Atom *atom1,
                              const RDKit::Atom *atom2,
                              const RDKit::Atom *atom3);

double computeAngle(const RDKit::RWMol &mol, const RDKit::Conformer &conf,
                    const RDKit::Atom *atom1, const RDKit::Atom *atom2,
                    RDKit::Atom *atom3);

double ComputeSmallestBondAngle(const RDKit::RWMol &mol,
                                const RDKit::Conformer &conf,
                                const RDKit::Atom *atom);

bool SmallestBondAngle(const RDKit::RWMol &mol, const RDKit::Conformer &conf,
                       const RDKit::Atom *atom);

void CleanUpMolecule(RDKit::RWMol &mol, RDKit::Conformer &conf);
