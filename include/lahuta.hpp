#include <vector>
#include "bonds.hpp"
#include "conv.hpp"
#include "nsgrid.hpp"
#include "bond_order.hpp"

namespace MolecularStructure {

class IMoleculeSource {
public:
    virtual RDKit::RWMol& getMolecule() = 0;
    virtual const RDKit::RWMol& getMolecule() const = 0;
    virtual RDKit::Conformer& getConformer() = 0;
    virtual const RDKit::Conformer& getConformer() const = 0;
    virtual ~IMoleculeSource() = default;

    virtual void process(const gemmi::Structure& st) = 0;
};

class GemmiStructureSource : public IMoleculeSource {
private:
    RDKit::RWMol &mol;
    // RDKit::Conformer conformer;

public:
    explicit GemmiStructureSource(RDKit::RWMol &mol) : mol(mol) {}

    void process(const gemmi::Structure& st) override {
        RDKit::Conformer *conformer = new RDKit::Conformer();
        gemmiStructureToRDKit(mol, st, *conformer, false);
        std::cout << "NUmber of coordinates: " << conformer->getPositions().size() << std::endl;
        mol.addConformer(conformer, true);
        std::cout << "xx Number of coordinates: \n";
        std::cout << mol.getConformer().getPositions().size() << std::endl;
    }

    RDKit::RWMol& getMolecule() override {
        return mol;
    }

    const RDKit::RWMol& getMolecule() const override {
        return mol;
    }

    RDKit::Conformer& getConformer() override {
        return mol.getConformer();
    }

    const RDKit::Conformer& getConformer() const override {
        return mol.getConformer();
    }

};

class NeighborSearch {
private:
    double cutoff;

public:
    explicit NeighborSearch(double cut = 4.5) : cutoff(cut) {}

    NSResults findNeighbors(const std::vector<RDGeom::Point3D>& coords) const {
        FastNS grid(coords, cutoff);
        return grid.selfSearch();
    }
};

class BondComputation {
public:
    static RDKit::RWMol computeProteinBonds(RDKit::RWMol& mol, const NSResults& neighbors) {
        std::vector<int> non_protein_indices;
        return lahutaBondAssignment(mol, neighbors, non_protein_indices);
    }

    static void perceiveBondOrders(RDKit::RWMol& mol) {
        mol.updatePropertyCache(false);
        PerceiveBondOrders(mol);
    }

    static void mergeBonds(RDKit::RWMol& targetMol, RDKit::RWMol& sourceMol, const std::vector<int>& indexMap) {
        for (auto bondIt = sourceMol.beginBonds(); bondIt != sourceMol.endBonds(); ++bondIt) {
            const RDKit::Bond* bond = *bondIt;
            int bIdx = indexMap[bond->getBeginAtomIdx()];
            int eIdx = indexMap[bond->getEndAtomIdx()];
            if (targetMol.getBondBetweenAtoms(bIdx, eIdx) == nullptr) {
                targetMol.addBond(bIdx, eIdx, bond->getBondType());
            }
        }
    }

    static void computeBonds(RDKit::RWMol& mol, const NSResults& neighborResults) {

        // RDKit::RWMol resultMol = mol;
        // auto proteinMol = computeProteinBonds(resultMol, neighborResults);

        std::cout << "Computing bonds..." << mol.getNumAtoms() << " " << mol.getNumBonds() << std::endl;
        std::vector<int> non_protein_indices;
        std::cout << "A\n";
        auto newMol = lahutaBondAssignment(mol, neighborResults, non_protein_indices);
        std::cout << "B\n";
        newMol.updatePropertyCache(false);
        std::cout << "C\n";
        perceiveBondOrders(newMol);
        std::cout << "D\n";

        for (auto bondIt = newMol.beginBonds(); bondIt != newMol.endBonds();
             ++bondIt) {
          RDKit::Bond *bond = *bondIt;
          auto bAtomIdx = bond->getBeginAtomIdx();
          auto eAtomIdx = bond->getEndAtomIdx();
          int bIdx = non_protein_indices[bAtomIdx];
          int eIdx = non_protein_indices[eAtomIdx];
          if (mol.getBondBetweenAtoms(bIdx, eIdx) == nullptr) {
            mol.addBond(bIdx, eIdx, bond->getBondType());
          }
        }
        std::cout << "E\n";

        // std::vector<int> non_protein_indices(resultMol.getNumAtoms());
        // std::iota(non_protein_indices.begin(), non_protein_indices.end(), 0);

        // mergeBonds(resultMol, proteinMol, non_protein_indices);
        // return resultMol;
    }
};

class MoleculeProcessor {
private:
    IMoleculeSource &moleculeSource;
    NeighborSearch neighborSearch;
    // RDKit::RWMol &molecule;
    std::vector<RDGeom::Point3D> atomCoords;
    NSResults neighborResults;

public:
    explicit MoleculeProcessor(IMoleculeSource &source, double cutoff = 4.5)
        : moleculeSource(source), neighborSearch(cutoff) {
        initializeMolecule();
    }

    void initializeMolecule() {
        // molecule = moleculeSource.getMolecule();
        std::cout << "Processing molecule..." << std::endl;
        std::cout << "Num atoms: " << moleculeSource.getMolecule().getNumAtoms() << std::endl;
        const auto& conf = moleculeSource.getConformer();
        std::cout << "-> Number of coordinates: \n";
        std::cout << conf.getPositions().size() << std::endl;
        atomCoords = conf.getPositions();
        neighborResults = neighborSearch.findNeighbors(atomCoords);
    }

    void processStructure(RDKit::RWMol& molecule) {
        BondComputation::computeBonds(molecule, neighborResults);
    }

    // RDKit::RWMol& getMolecule() { return molecule; }
    // const RDKit::RWMol& getMolecule() const { return molecule; }
    const std::vector<RDGeom::Point3D>& getAtomCoords() const { return atomCoords; }
    const NSResults& getNeighborResults() const { return neighborResults; }
};

} // namespace MolecularStructure

// Usage example:
// auto source = std::make_unique<MolecularStructure::GemmiStructureSource>(gemmiStructure);
// MolecularStructure::MoleculeProcessor processor(std::move(source));
// processor.processStructure();
// auto processedMolecule = processor.getMolecule();
