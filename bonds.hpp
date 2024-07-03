#include <GraphMol/PeriodicTable.h>
#include <GraphMol/RWMol.h>
#include <gemmi/neighbor.hpp>

#include "nsgrid.hpp"

using namespace gemmi;

struct BondInfo {
  const_CRA cra1, cra2;
  int image_idx;
  double dist_sq;
  int order;

  void print(const UnitCell &cell) const {
    NearestImage im =
        cell.find_nearest_pbc_image(cra1.atom->pos, cra2.atom->pos, image_idx);
    assert(fabs(dist_sq - im.dist_sq) < 1e-3);
    std::printf("%s - %s  im:%s  %.3f   %d\n", atom_str(cra1).c_str(),
                atom_str(cra2).c_str(), im.symmetry_code(true).c_str(),
                std::sqrt(dist_sq), order);
  }
};
void perceiveBonds(RDKit::RWMol &mol, const NSResults &results,
                   const float tolerance);

void findBondsDeconstructed(Structure &st, Model &model, RDKit::RWMol &mol,
                            double maxRadius, double covFactor);
