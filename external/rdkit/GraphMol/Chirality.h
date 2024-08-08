//
//  Copyright (C) 2008-2022 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
/*! \file Chirality.h

*/
#ifndef RD_CHIRALITY_20AUG2008_H
#define RD_CHIRALITY_20AUG2008_H
#include <RDGeneral/types.h>
#include <GraphMol/Bond.h>
#include <boost/dynamic_bitset.hpp>
#include <limits>

namespace RDKit {
class Atom;
class Bond;
class ROMol;
class Conformer;

namespace Chirality {

//! double bond stereo will be ignored/removed for rings smaller than this:
constexpr unsigned int minRingSizeForDoubleBondStereo = 8;

constexpr auto nonTetrahedralStereoEnvVar = "RDK_ENABLE_NONTETRAHEDRAL_STEREO";
constexpr auto useLegacyStereoEnvVar = "RDK_USE_LEGACY_STEREO_PERCEPTION";
constexpr bool nonTetrahedralStereoDefaultVal =
    true;  //!< whether or not nontetrahedral stereo is perceived by default
constexpr bool useLegacyStereoDefaultVal =
    true;  //!< whether or not the legacy stereo perception code is used by
           //!< default

extern void setAllowNontetrahedralChirality(bool val);
extern bool getAllowNontetrahedralChirality();

extern void setUseLegacyStereoPerception(bool val);
extern bool getUseLegacyStereoPerception();

void removeNonExplicit3DChirality(ROMol &mol);

extern bool
    useLegacyStereoPerception;  //!< Toggle usage of the legacy stereo
                                //!< perception code

extern bool
    useLegacyStereoPerception;  //!< Toggle usage of the legacy stereo
                                //!< perception code

/// @cond
/*!
  \param mol the molecule to be altered
  \param ranks  used to return the set of ranks.
                Should be at least mol.getNumAtoms() long.

  <b>Notes:</b>
     - All atoms gain a property common_properties::_CIPRank with their overall
       CIP ranking.

*/
void assignAtomCIPRanks(const ROMol &mol,
                                              UINT_VECT &ranks);

bool hasStereoBondDir(const Bond *bond);

/**
 *  Returns the first neighboring bond that can be found which has a stereo
 * bond direction set. If no such bond can be found, it returns null. No
 * checks are made to ensure there aren't any other conflicting directed bonds.
 */
const Bond *getNeighboringDirectedBond(const ROMol &mol,
                                                             const Atom *atom);

/**
 *  This just translates the labels, setting/translating StereoAtoms or the
 * label is not the responsibility of this function. If the passed label is not
 * E/Z, it will be returned unchanged.
 */
Bond::BondStereo translateEZLabelToCisTrans(
    Bond::BondStereo label);
/// @endcond

enum class StereoType {
  Unspecified,
  Atom_Tetrahedral,
  Atom_SquarePlanar,
  Atom_TrigonalBipyramidal,
  Atom_Octahedral,
  Bond_Double,         // single double bond and odd-numbered cumulenes
  Bond_Cumulene_Even,  // even-numbered cumulenes
  Bond_Atropisomer
};

enum class StereoDescriptor {
  None,
  Tet_CW,
  Tet_CCW,
  Bond_Cis,
  Bond_Trans,
  Bond_AtropCW,
  Bond_AtropCCW
};

enum class StereoSpecified {
  Unspecified,  // no information provided
  Specified,
  Unknown  // deliberately marked as unknown
};

struct StereoInfo {
  // REVIEW: absolute stereo data member?
#ifdef _MSC_VER
  static const unsigned NOATOM =
      std::numeric_limits<unsigned>::max();  // used to mark missing atoms
#else
  static const unsigned NOATOM;  // used to mark missing atoms
#endif
  StereoType type = StereoType::Unspecified;
  StereoSpecified specified = StereoSpecified::Unspecified;
  unsigned centeredOn = NOATOM;
  StereoDescriptor descriptor = StereoDescriptor::None;
  unsigned permutation = 0;  // for the non-tetrahedral stereo cases
  std::vector<unsigned> controllingAtoms;  // all atoms around the atom or bond.
  // Order is important
  bool operator==(const StereoInfo &other) const {
    return type == other.type && specified == other.specified &&
           centeredOn == other.centeredOn && descriptor == other.descriptor &&
           permutation == other.permutation &&
           controllingAtoms == other.controllingAtoms;
  }
  bool operator!=(const StereoInfo &other) const { return !(*this == other); }
};

//! identifies potential stereoatoms and stereobonds in a molecule
/*!
  Note that this function is still somewhat experimental and the API
  and results may change in a future release.

  \param mol the molecule to look for stereo in
  \param cleanIt remove chirality/stereo specifications from atoms/bonds that
     cannot be chiral/stereo
*/
std::vector<StereoInfo> findPotentialStereo(
    ROMol &mol, bool cleanIt, bool flagPossible = true);
//! overload
std::vector<StereoInfo> findPotentialStereo(
    const ROMol &mol);

//! removes atoms without specified chirality from stereo groups
void cleanupStereoGroups(ROMol &mol);

//! calls the approximate legacy code for assigning CIP labels
void assignLegacyCIPLabels(
    ROMol &mol, bool flagPossibleStereoCenters = false);

/// @cond
namespace detail {
bool isAtomPotentialNontetrahedralCenter(
    const Atom *atom);
bool isAtomPotentialTetrahedralCenter(const Atom *atom);
bool isAtomPotentialStereoAtom(const Atom *atom);
bool isBondPotentialStereoBond(const Bond *bond);
StereoInfo getStereoInfo(const Bond *bond);
StereoInfo getStereoInfo(const Atom *atom);
bool bondAffectsAtomChirality(const Bond *bond,
                                                    const Atom *atom);
unsigned int getAtomNonzeroDegree(const Atom *atom);

bool has_protium_neighbor(const ROMol &mol,
                                                const Atom *atom);

}  // namespace detail
/// @endcond

INT_VECT findStereoAtoms(const Bond *bond);

//! \name Non-tetrahedral stereochemistry
//! @{
bool hasNonTetrahedralStereo(const Atom *center);
Bond *getChiralAcrossBond(const Atom *center,
                                                const Bond *qry);
Bond *getChiralAcrossBond(const Atom *center,
                                                const Atom *qry);
Atom *getChiralAcrossAtom(const Atom *center,
                                                const Bond *qry);
Atom *getChiralAcrossAtom(const Atom *center,
                                                const Atom *qry);
//! \param which: if this is -1 then the second axial bond will be returned,
//! otherwise the first
Bond *getTrigonalBipyramidalAxialBond(const Atom *center,
                                                            int which = 0);
Atom *getTrigonalBipyramidalAxialAtom(const Atom *center,
                                                            int which = 0);

//! \returns 1 if it's the first axial atom, -1 if it's the second
int isTrigonalBipyramidalAxialBond(const Atom *center,
                                                         const Bond *qry);
int isTrigonalBipyramidalAxialAtom(const Atom *center,
                                                         const Atom *qry);

double getIdealAngleBetweenLigands(const Atom *center,
                                                         const Atom *lig1,
                                                         const Atom *lig2);

unsigned int getChiralPermutation(const Atom *center,
                                                        const INT_LIST &probe);
//! @}

std::ostream &operator<<(std::ostream &oss,
                                               const StereoSpecified &s);
std::ostream &operator<<(std::ostream &oss,
                                               const StereoType &s);

struct BondWedgingParameters {
  bool wedgeTwoBondsIfPossible =
      false;  //!< If this is enabled then two bonds will be wedged at chiral
              //!< centers subject to the following constraints:
              //!<   1. ring bonds will not be wedged
              //!<   2. bonds to chiral centers will not be wedged
              //!<   3. bonds separated by more than 120 degrees will not be
              //!<      wedged
};

enum class WedgeInfoType {
  WedgeInfoTypeChiral,
  WedgeInfoTypeAtropisomer,
};

class WedgeInfoBase {
 public:
  WedgeInfoBase(int idxInit) : idx(idxInit){};
  virtual ~WedgeInfoBase(){};

  virtual WedgeInfoType getType() const = 0;
  virtual Bond::BondDir getDir() const = 0;

  int getIdx() const { return idx; }

 private:
  int idx = -1;
};

class WedgeInfoChiral : public WedgeInfoBase {
 public:
  WedgeInfoChiral(int atomId) : WedgeInfoBase(atomId){};
  ~WedgeInfoChiral(){};

  WedgeInfoType getType() const override {
    return Chirality::WedgeInfoType::WedgeInfoTypeChiral;
  }
  Bond::BondDir getDir() const override {
    throw std::runtime_error(
        "BondDir is not stored/used in Chiral type WedgInfos");
  }
};

class WedgeInfoAtropisomer : public WedgeInfoBase {
 public:
  WedgeInfoAtropisomer(int bondId, RDKit::Bond::BondDir dirInit)
      : WedgeInfoBase(bondId) {
    dir = dirInit;
  };
  ~WedgeInfoAtropisomer(){};

  RDKit::Bond::BondDir dir = RDKit::Bond::BondDir::NONE;

  WedgeInfoType getType() const override {
    return Chirality::WedgeInfoType::WedgeInfoTypeAtropisomer;
  }

  Bond::BondDir getDir() const override { return dir; }
};

namespace detail {
Bond::BondDir determineBondWedgeState(
    const Bond *bond, unsigned int fromAtomIdx, const Conformer *conf);
Bond::BondDir determineBondWedgeState(
    const Bond *bond,
    const std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
        &wedgeBonds,
    const Conformer *conf);

std::pair<bool, INT_VECT> countChiralNbrs(
    const ROMol &mol, int noNbrs);
int pickBondToWedge(
    const Atom *atom, const ROMol &mol, const INT_VECT &nChiralNbrs,
    const std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
        &resSoFar,
    int noNbrs);
void setStereoForBond(ROMol &mol, Bond *bond,
                                            Bond::BondStereo stereo);
}  // namespace detail

//! picks the bonds which should be wedged
/// returns a map from bond idx -> controlling atom idx
std::map<int, std::unique_ptr<Chirality::WedgeInfoBase>> pickBondsToWedge(
    const ROMol &mol, const BondWedgingParameters *params = nullptr);

std::map<int, std::unique_ptr<Chirality::WedgeInfoBase>> pickBondsToWedge(
    const ROMol &mol, const BondWedgingParameters *params,
    const Conformer *conf);

void wedgeMolBonds(
    ROMol &mol, const Conformer *conf = nullptr,
    const BondWedgingParameters *params = nullptr);
void wedgeBond(Bond *bond, unsigned int fromAtomIdx,
                                     const Conformer *conf);

//! Returns true for double bonds which should be shown as a crossed bonds.
// It always returns false if any adjacent bond is a squiggle bond.
bool shouldBeACrossedBond(const Bond *bond);

//! Clears existing bond wedging and forces use of atom wedging from MolBlock.
/*!
 \param mol: molecule to have its wedges altered
 \param allBondTypes: reapply the wedging also on bonds other than single and aromatic ones
 */
void reapplyMolBlockWedging(ROMol &mol, bool allBondTypes=true);
//! Remove MolBlock bond wedging information from molecule.
/*!
 \param mol: molecule to modify
 */
void clearMolBlockWedgingInfo(ROMol &mol);
//! Invert bond wedging information read from a mol block (if present).
/*!
 \param mol: molecule to modify
 */
void invertMolBlockWedgingInfo(ROMol &mol);

//! gets stereo info for a bond
/*!
 \param bond: bond to check
 \param wedgeBonds - the list of bonds to have wedges
 \param conf -  Conformer to use
 \param dirCode - receives the dircode for the bond
 \param reverse - receives the reverse flag
 only returned if it was exlicility set witha wiggle bond
 */

void GetMolFileBondStereoInfo(
    const Bond *bond,
    const std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
        &wedgeBonds,
    const Conformer *conf, int &dirCode, bool &reverse);

void GetMolFileBondStereoInfo(
    const Bond *bond,
    const std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
        &wedgeBonds,
    const Conformer *conf, Bond::BondDir &dir, bool &reverse);

//! add R/S, relative stereo, and E/Z annotations to atoms and bonds
/*!
 \param mol: molecule to modify
 \param absLabel: label for atoms in an ABS stereo group
 \param orLabel: label for atoms in an OR stereo group
 \param andLabel: label for atoms in an AND stereo group
 \param cipLabel: label for chiral atoms that aren't in a stereo group.
 \param bondLabel: label for CIP stereochemistry on bonds

 If any label is empty, the corresponding annotations will not be added.

 The labels can contain the following placeholders:
   {id} - the stereo group's index
   {cip} - the atom or bond's CIP stereochemistry

 Note that CIP labels will only be added if CIP stereochemistry has been
 assigned to the molecule.

 */
void addStereoAnnotations(
    ROMol &mol, std::string absLabel = "abs ({cip})",
    std::string orLabel = "or{id}", std::string andLabel = "and{id}",
    std::string cipLabel = "({cip})", std::string bondLabel = "({cip})");

}  // namespace Chirality
}  // namespace RDKit
#endif
