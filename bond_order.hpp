#include "GraphMol/RDKitBase.h"
#include <GraphMol/Substruct/SubstructMatch.h>

using namespace RDKit;

using HybridizationType = RDKit::Atom::HybridizationType;
using SubStrMatches = std::vector<RDKit::MatchVectType>;

int GetExpVal(const RDKit::Atom *atom);

constexpr std::pair<const char *, HybridizationType> smartsList[] = {
    {"[D4]", HybridizationType::SP3},
    {"[D5]", HybridizationType::SP3D},
    {"[D6]", HybridizationType::SP3D2},
    {"[C]", HybridizationType::SP3},
    {"[c,$(C=*)]", HybridizationType::SP2},
    {"[$(C#*),$(C(=*)=*)]", HybridizationType::SP},

    {"[N]", HybridizationType::SP3},
    {"[n,$(N=*),$(N[#6,#7,#8]=,:,#*)]", HybridizationType::SP2},
    {"[ND1,ND2,ND3]a", HybridizationType::SP2},
    {"[$(N#*),$([ND2](=*)=*)]", HybridizationType::SP},

    {"[O]", HybridizationType::SP3},
    {"[o,$(O=*),$(O[#6,#7,#8]=,:*)]", HybridizationType::SP2},
    {"[$([#8D1][#6][#8D1])]", HybridizationType::SP2},
    {"[$(O#*)]", HybridizationType::SP},

    {"[S]", RDKit::Atom::HybridizationType::SP3},
    {"[#16;s,$([SD1]=*)]", RDKit::Atom::HybridizationType::SP2},
    {"[SD6]", RDKit::Atom::HybridizationType::SP3D2},
};

using SmartsPatternPair = std::pair<const char *, std::vector<int>>;
const std::vector<SmartsPatternPair> bondSmarts = {
    {"[x2,x3]1[#6]([#7D3]2)[#6][#6][#6]2[x2,x3][#6]([#7D3]3)[#6][#6][#6]3["
     "x2,"
     "x3][#6]([#7D3]4)[#6][#6][#6]4[x2,x3][#6]([#7D3]5)[#6][#6][#6]51",
     {0,  1,  2,  1,  2,  1,  1,  3,  1,  3,  4,  2,  4,  5,  1,  5,  2,
      1,  5,  6,  2,  6,  7,  1,  7,  8,  2,  7,  9,  1,  9,  10, 2,  10,
      11, 1,  11, 8,  1,  11, 12, 2,  12, 13, 1,  13, 14, 1,  13, 15, 2,
      15, 16, 1,  16, 17, 2,  17, 14, 1,  17, 18, 1,  18, 19, 2,  19, 20,
      1,  19, 21, 1,  21, 22, 2,  22, 23, 1,  23, 20, 2}},
    {"[x2,x3]1[#6]([#7D3]2)[#6][#6][#6]2[x2,x3][#6]([#7]3)[#6][#6][#6]3[x2,"
     "x3][#6]([#7D3]4)[#6][#6][#6]4[x2,x3][#6]([#7]5)[#6][#6][#6]51",
     {0,  1,  2,  1,  2,  1,  1,  3,  1,  3,  4,  2,  4,  5,  1,  5,  2,
      1,  5,  6,  2,  6,  7,  1,  7,  8,  2,  7,  9,  1,  9,  10, 2,  10,
      11, 1,  11, 8,  1,  11, 12, 2,  12, 13, 1,  13, 14, 1,  13, 15, 2,
      15, 16, 1,  16, 17, 2,  17, 14, 1,  17, 18, 1,  18, 19, 2,  19, 20,
      1,  19, 21, 1,  21, 22, 2,  22, 23, 1,  23, 20, 2}},
    {"[x2,x3]1[#6]([#7]2)[#6][#6][#6]2[x2,x3][#6]([#7]3)[#6][#6][#6]3[x2,"
     "x3]["
     "#6]([#7]4)[#6][#6][#6]4[x2,x3][#6]([#7]5)[#6][#6][#6]51",
     {0,  1,  2,  1,  2,  1,  1,  3,  1,  3,  4,  2,  4,  5,  1,  5,  2,
      1,  5,  6,  2,  6,  7,  1,  7,  8,  2,  7,  9,  1,  9,  10, 2,  10,
      11, 1,  11, 8,  1,  11, 12, 2,  12, 13, 1,  13, 14, 1,  13, 15, 2,
      15, 16, 1,  16, 17, 2,  17, 14, 1,  17, 18, 1,  18, 19, 2,  19, 20,
      1,  19, 21, 1,  21, 22, 2,  22, 23, 1,  23, 20, 2}},
    {"[#7D2][#7D2^1][#7D1]", {0, 1, 2, 1, 2, 2}},
    {"[#8D1][#7D3^2]([#8D1])*", {0, 1, 2, 1, 2, 2, 1, 3, 1}},
    {"[#16D4]([#8D1])([#8D1])([*!#8])([*!#8])",
     {0, 1, 2, 0, 2, 2, 0, 3, 1, 0, 4, 1}},
    {"[#16D4]([#8D1])([#8D1])([#8-,#8D1])([#8-,#8D1])",
     {0, 1, 2, 0, 2, 2, 0, 3, 1, 0, 4, 1}},
    {"[#16D4]([#16D1])([#8D1])([#8-,#8])([#8-,#8])",
     {0, 1, 2, 0, 2, 2, 0, 3, 1, 0, 4, 1}},
    {"[#16D3]([#8D1])([*!#8])([*!#8])", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    {"[#16D3]([#8D1])([#8D1-])([#8D1-])", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    {"[#16D3]([#8D1])([#8])([#8])", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    {"[#16D2]([#8D1])([#16D1])", {0, 1, 2, 0, 2, 2}},
    {"[#16D2]([#8D1])([*!#8])", {0, 1, 2, 0, 2, 1}},
    {"[#16D2]([#8D1])([#8D1])", {0, 1, 2, 0, 2, 2}},
    {"[#15D3]([#8D1])([#8D1])([#8D2])", {0, 1, 2, 0, 2, 2, 0, 3, 1}},
    {"[#7D2]([#8D1])([#1])", {0, 1, 2, 0, 2, 1}},
    {"[#15D4]([#8D1])(*)(*)(*)", {0, 1, 2, 0, 2, 1, 0, 3, 1, 0, 4, 1}},
    {"[#6D3^2]([#8D1])([#8])*", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    {"[#8D1][#6D2^1][#8D1]", {0, 1, 2, 1, 2, 2}},
    {"[#6D3^2]([#8D1;!-])([#7])*", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    {"[#34D3^2]([#8D1])([#8])*", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    {"[#6D3^2]([#8D1])([#16])*", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    {"[#6D3^2]([#16D1])([#16])*", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    {"[CD3^2]([#16D1])([N])*", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    // Do not seem to work correctly with RDKit:
    // {"[#6^2][#6D2^1][#6^2]", {0, 1, 2, 1, 2, 2}},
    // {"[#6^2][#6D2^1][#8D1]", {0, 1, 2, 1, 2, 2}},
    {"[#6D3^2;!R]([#7D1H0;!R])([#7;!R])*", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    {"[#6D3^2;!R]([#7D2H1;!R])([#7;!R])*", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    {"[#6D3^2;!R]([#7D3H2;!R])([#7;!R])*", {0, 1, 2, 0, 2, 1, 0, 3, 1}},
    {"[#6D3^2;!R]([#1,#6])([#1,#6])[#7D3^2;!R]([#1])[#6]",
     {0, 1, 1, 0, 2, 1, 0, 3, 2, 3, 4, 1, 3, 5, 1}},
};

SubStrMatches performSubstructMatch(RDKit::ROMol &mol, RDKit::ROMol &pattern,
                                    SubstructMatchParameters &params);

void RDKitSmartsMatch(RDKit::ROMol &mol, SubstructMatchParameters &params);

void OBBondTypeAssignment(RDKit::ROMol &mol);

double AverageBondAngle(const RDKit::Atom *atom);

/*!  This function calculates the torsion angle of three vectors, represented
  by four points A--B--C--D, i.e. B and C are vertexes, but none of A--B,
  B--C, and C--D are colinear.  A "torsion angle" is the amount of "twist"
  or torsion needed around the B--C axis to bring A--B into the same plane
  as B--C--D.  The torsion is measured by "looking down" the vector B--C so
  that B is superimposed on C, then noting how far you'd have to rotate
  A--B to superimpose A over D.  Angles are + in the anti-clockwise
  direction.  The operation is symmetrical in that if you reverse the image
  (look from C to B and rotate D over A), you get the same answer.
*/
double CalcTorsionAngle(const RDGeom::Point3D &a, const RDGeom::Point3D &b,
                        const RDGeom::Point3D &c, const RDGeom::Point3D &d);

void PerceiveBondOrders(RDKit::RWMol &mol);
