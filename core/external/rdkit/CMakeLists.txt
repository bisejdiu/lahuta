cmake_minimum_required(VERSION 3.15...3.22)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(rdkit STATIC
            RDGeneral/Invariant.cpp RDGeneral/LocaleSwitcher.cpp
            RDGeneral/RDLog.cpp RDGeneral/types.cpp
            RDGeneral/utils.cpp

            GraphMol/atomic_data.cpp
            GraphMol/PeriodicTable.cpp

            DataStructs/base64.cpp DataStructs/BitOps.cpp
            DataStructs/BitVect.cpp DataStructs/ExplicitBitVect.cpp
            DataStructs/SparseBitVect.cpp 
            Geometry/point.cpp
            Geometry/Transform3D.cpp
            GraphMol/Atom.cpp 

            GraphMol/AtomIterators.cpp GraphMol/Atropisomers.cpp
            GraphMol/Bond.cpp GraphMol/BondIterators.cpp
            GraphMol/Canon.cpp GraphMol/Chirality.cpp
            GraphMol/Conformer.cpp
            GraphMol/Kekulize.cpp
            GraphMol/NontetrahedralStereo.cpp
            GraphMol/AddHs.cpp
            GraphMol/ConjugHybrid.cpp
            GraphMol/Aromaticity.cpp
            GraphMol/FindStereo.cpp
            GraphMol/Matrices.cpp
            GraphMol/AdjustQuery.cpp
            GraphMol/WedgeBonds.cpp
            GraphMol/FindRings.cpp

            GraphMol/GraphMol.cpp GraphMol/MolOps.cpp
            GraphMol/MonomerInfo.cpp GraphMol/new_canon.cpp
            GraphMol/QueryAtom.cpp
            GraphMol/QueryBond.cpp GraphMol/QueryOps.cpp
            GraphMol/Resonance.cpp GraphMol/RingInfo.cpp
            GraphMol/ROMol.cpp GraphMol/RWMol.cpp
            GraphMol/StereoGroup.cpp GraphMol/SubstanceGroup.cpp

            GraphMol/CIPLabeler/CIPLabeler.cpp GraphMol/CIPLabeler/CIPMol.cpp
            GraphMol/CIPLabeler/Digraph.cpp GraphMol/CIPLabeler/Edge.cpp
            GraphMol/CIPLabeler/Mancude.cpp GraphMol/CIPLabeler/Node.cpp
            GraphMol/CIPLabeler/Sort.cpp 
            GraphMol/CIPLabeler/configs/AtropisomerBond.cpp
            GraphMol/CIPLabeler/configs/Configuration.cpp
            GraphMol/CIPLabeler/configs/Sp2Bond.cpp
            GraphMol/CIPLabeler/configs/Tetrahedral.cpp

            GraphMol/CIPLabeler/rules/Rule1a.cpp
            GraphMol/CIPLabeler/rules/Rule1b.cpp
            GraphMol/CIPLabeler/rules/Rule2.cpp
            GraphMol/CIPLabeler/rules/Rule3.cpp
            GraphMol/CIPLabeler/rules/Rule4a.cpp
            GraphMol/CIPLabeler/rules/Rule4b.cpp
            GraphMol/CIPLabeler/rules/Rule4c.cpp
            GraphMol/CIPLabeler/rules/Rule5.cpp
            GraphMol/CIPLabeler/rules/Rule5New.cpp
            GraphMol/CIPLabeler/rules/Rule6.cpp
            GraphMol/CIPLabeler/rules/SequenceRule.cpp

            GraphMol/GenericGroups/GenericGroups.cpp
            GraphMol/GraphMol.cpp GraphMol/CIPLabeler/Mancude.cpp
            GraphMol/FileParsers/MolFileStereochem.cpp
            GraphMol/MolOps.cpp GraphMol/MonomerInfo.cpp
            GraphMol/new_canon.cpp GraphMol/CIPLabeler/Node.cpp
            GraphMol/QueryAtom.cpp
            GraphMol/QueryBond.cpp GraphMol/QueryOps.cpp
            GraphMol/Resonance.cpp GraphMol/RingInfo.cpp
            GraphMol/ROMol.cpp GraphMol/RWMol.cpp

            GraphMol/SmilesParse/SmilesParse.cpp
            GraphMol/SmilesParse/SmilesParseOps.cpp
            GraphMol/SmilesParse/CXSmilesOps.cpp
            GraphMol/SmilesParse/SmilesWrite.cpp

            GraphMol/SmilesParse/smarts.tab.cpp
            GraphMol/SmilesParse/smiles.tab.cpp

            GraphMol/SmilesParse/lex.yysmiles.cpp
            GraphMol/SmilesParse/lex.yysmarts.cpp

            GraphMol/CIPLabeler/Sort.cpp
            GraphMol/StereoGroup.cpp GraphMol/SubstanceGroup.cpp
            GraphMol/Substruct/SubstructMatch.cpp
            GraphMol/Substruct/SubstructUtils.cpp
)

set_property(TARGET rdkit PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(rdkit PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR} # rdkit headers
)
# Consume Boost from parent: either vendored (Boost::headers created by parent)
# or system (Boost::headers provided by find_package)
target_link_libraries(rdkit PUBLIC Boost::headers)
