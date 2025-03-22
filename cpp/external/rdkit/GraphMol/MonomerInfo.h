//
//  Copyright (C) 2013-2024 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
/*! \file MonomerInfo.h

  \brief Defines Monomer information classes

*/
#ifndef RD_MONOMERINFO_H
#define RD_MONOMERINFO_H

#include <string>
#include <utility>
#include <boost/shared_ptr.hpp>

namespace RDKit {

//! The abstract base class for atom-level monomer info
class AtomMonomerInfo {
 public:
  typedef enum { UNKNOWN = 0, PDBRESIDUE, OTHER } AtomMonomerType;

  virtual ~AtomMonomerInfo() {}
  virtual void destroy() { delete this; }

  AtomMonomerInfo() = default;
  AtomMonomerInfo(AtomMonomerType typ, std::string nm = "")
      : d_monomerType(typ), d_name(std::move(nm)) {}
  AtomMonomerInfo(const AtomMonomerInfo &other) = default;

  const std::string &getName() const { return d_name; }
  void setName(const std::string &nm) { d_name = nm; }
  AtomMonomerType getMonomerType() const { return d_monomerType; }
  void setMonomerType(AtomMonomerType typ) { d_monomerType = typ; }

  virtual AtomMonomerInfo *copy() const { return new AtomMonomerInfo(*this); }

 private:
  AtomMonomerType d_monomerType{UNKNOWN};
  std::string d_name{""};
};

//! Captures atom-level information about peptide residues
class AtomPDBResidueInfo : public AtomMonomerInfo {
 public:
  AtomPDBResidueInfo() : AtomMonomerInfo(PDBRESIDUE) {}
  AtomPDBResidueInfo(const AtomPDBResidueInfo &other) = default;

  AtomPDBResidueInfo(const std::string &atomName, int serialNumber = 0,
                     std::string altLoc = "", std::string residueName = "",
                     int residueNumber = 0, std::string chainId = "",
                     std::string insertionCode = "", double occupancy = 1.0,
                     double tempFactor = 0.0, bool isHeteroAtom = false,
                     unsigned int secondaryStructure = 0,
                     unsigned int segmentNumber = 0)
      : AtomMonomerInfo(PDBRESIDUE, atomName),
        d_serialNumber(serialNumber),
        d_altLoc(std::move(altLoc)),
        d_residueName(std::move(residueName)),
        d_residueNumber(residueNumber),
        d_chainId(std::move(chainId)),
        d_insertionCode(std::move(insertionCode)),
        d_occupancy(occupancy),
        d_tempFactor(tempFactor),
        df_heteroAtom(isHeteroAtom),
        d_secondaryStructure(secondaryStructure),
        d_segmentNumber(segmentNumber) {}

  int getSerialNumber() const { return d_serialNumber; }
  void setSerialNumber(int val) { d_serialNumber = val; }
  const std::string &getAltLoc() const { return d_altLoc; }
  void setAltLoc(const std::string &val) { d_altLoc = val; }
  const std::string &getResidueName() const { return d_residueName; }
  void setResidueName(const std::string &val) { d_residueName = val; }
  int getResidueNumber() const { return d_residueNumber; }
  void setResidueNumber(int val) { d_residueNumber = val; }
  const std::string &getChainId() const { return d_chainId; }
  void setChainId(const std::string &val) { d_chainId = val; }
  const std::string &getInsertionCode() const { return d_insertionCode; }
  void setInsertionCode(const std::string &val) { d_insertionCode = val; }
  double getOccupancy() const { return d_occupancy; }
  void setOccupancy(double val) { d_occupancy = val; }
  double getTempFactor() const { return d_tempFactor; }
  void setTempFactor(double val) { d_tempFactor = val; }
  bool getIsHeteroAtom() const { return df_heteroAtom; }
  void setIsHeteroAtom(bool val) { df_heteroAtom = val; }
  unsigned int getSecondaryStructure() const { return d_secondaryStructure; }
  void setSecondaryStructure(unsigned int val) { d_secondaryStructure = val; }
  unsigned int getSegmentNumber() const { return d_segmentNumber; }
  void setSegmentNumber(unsigned int val) { d_segmentNumber = val; }

  AtomMonomerInfo *copy() const override {
    return static_cast<AtomMonomerInfo *>(new AtomPDBResidueInfo(*this));
  }

  unsigned int getResidueIndex() const { return d_residueIndex; }
  void setResidueIndex(unsigned int idx) { d_residueIndex = idx; }

  friend class LeanAtomPDBResidueInfo;
  friend class AtomPDBResidueInfoStatic;

 private:
  // the fields here are from the PDB definition
  // (http://www.wwpdb.org/documentation/format33/sect9.html#ATOM) [9 Aug, 2013]
  // element and charge are not present since the atom itself stores that
  // information
  unsigned int d_serialNumber = 0;
  std::string d_altLoc = "";
  std::string d_residueName = "";
  int d_residueNumber = 0;
  std::string d_chainId = "";
  std::string d_insertionCode = "";
  double d_occupancy = 1.0;
  double d_tempFactor = 0.0;
  // additional, non-PDB fields:
  bool df_heteroAtom = false;  // is this from a HETATM record?
  unsigned int d_secondaryStructure = 0;
  unsigned int d_segmentNumber = 0;

  unsigned int d_residueIndex = 0;
};

class LeanAtomPDBResidueInfo : public AtomPDBResidueInfo {
public:
  LeanAtomPDBResidueInfo(const char* atomName, int serialNumber,
                        const char* residueName, int residueNumber)
      : AtomPDBResidueInfo() {
        d_atomNamePtr = atomName;
        d_residueNamePtr = residueName;
        d_serialNumber = serialNumber;
        d_residueNumber = residueNumber;
  }

  void destroy() override {
    // explicit
    this->~LeanAtomPDBResidueInfo();
    // This object was constructed using placement-new and memory is managed externally.
  }


  const std::string& getName() const {
    if (d_cachedName.empty() && d_atomNamePtr) {
      d_cachedName = d_atomNamePtr;
    }
    return d_cachedName;
  }

  const std::string& getResidueName() const {
    if (d_cachedResidueName.empty() && d_residueNamePtr) {
      d_cachedResidueName = d_residueNamePtr;
    }
    return d_cachedResidueName;
  }

  const std::string& getChainId() const {
    static const std::string defaultChain = "A";
    return defaultChain;
  }

  const std::string& getAltLoc() const {
    static const std::string emptyString = "";
    return emptyString;
  }

  const std::string& getInsertionCode() const {
    static const std::string emptyString = "";
    return emptyString;
  }

  bool getIsHeteroAtom() const { return false; }
  AtomMonomerType getMonomerType() const { return AtomPDBResidueInfo::PDBRESIDUE; }

  RDKit::AtomMonomerInfo* copy() const override {
    return static_cast<RDKit::AtomMonomerInfo*>(new LeanAtomPDBResidueInfo(*this));
  }

private:
  // pointers to strings, no copies
  const char* d_atomNamePtr = nullptr;
  const char* d_residueNamePtr = nullptr;
  // optional cache stuff
  mutable std::string d_cachedName;
  mutable std::string d_cachedResidueName;
};


struct AtomStaticInfo {
  const char *atom_name;
  const char *residue_name;
  int local_atom_index; // Atom position within residue
};

class AtomPDBResidueInfoStatic : public RDKit::AtomMonomerInfo {
public:
  const AtomStaticInfo* static_info; // pointer to compile-time static data
  int residue_number;
  char chain_id;

  AtomPDBResidueInfoStatic(const AtomStaticInfo* staticInfo, int residueNum, char chainId = 'A')
    : AtomMonomerInfo(PDBRESIDUE, staticInfo->atom_name),
      static_info(staticInfo),
      residue_number(residueNum),
      chain_id(chainId) {}

  AtomMonomerInfo* copy() const override {
    return new AtomPDBResidueInfoStatic(*this);
  }

  /*const std::string& getResidueName() const { return static_info->residue_name; }*/
  const std::string getResidueName() const { return static_info->residue_name; }
  int getResidueNumber() const { return residue_number; }

  const std::string& getChainId() const {
    static const std::string defaultChain = "A";
    return defaultChain;
  }

  const std::string& getAltLoc() const {
    static const std::string emptyString = "";
    return emptyString;
  }

  const std::string& getInsertionCode() const {
    static const std::string emptyString = "";
    return emptyString;
  }

  bool getIsHeteroAtom() const { return false; }
  AtomMonomerType getMonomerType() const { return AtomPDBResidueInfo::PDBRESIDUE; }

  void setResidueNumber(int val) { residue_number = val; }
  void setChainId(char val) { chain_id = val; }
};



};  // namespace RDKit
//! allows AtomPDBResidueInfo objects to be dumped to streams
std::ostream &operator<<(
    std::ostream &target, const RDKit::AtomPDBResidueInfo &apri);

#endif
