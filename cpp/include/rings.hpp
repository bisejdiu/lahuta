#ifndef LAHUTA_RINGS_HPP
#define LAHUTA_RINGS_HPP

#include <vector>
#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/RWMol.h>

namespace lahuta {

struct RingData {
    std::vector<int> atom_ids;
    RDGeom::Point3D center;
    RDGeom::Point3D norm1;
    RDGeom::Point3D norm2;

    RingData() = default;
    RingData(RDGeom::Point3D center, RDGeom::Point3D norm1, RDGeom::Point3D norm2,
             std::vector<int> atom_ids)
        : center(center), norm1(norm1), norm2(norm2), atom_ids(atom_ids) {}
};

struct RingDataVec {
  std::vector<RingData> rings;

  RingDataVec() = default;
  RingDataVec(RingDataVec &&other) noexcept : rings(std::move(other.rings)) {}
  RingDataVec &operator=(RingDataVec &&other) noexcept {
    if (this != &other) {
      rings = std::move(other.rings);
    }
    return *this;
  }

  // RingDataVec(std::vector<RingData> rings) : rings(std::move(rings)) {}
  // // copy constructor
  RingDataVec(const RingDataVec &other) : rings(other.rings) {}
  // copy = operator
  RingDataVec &operator=(const RingDataVec &other) {
    if (this != &other) {
      rings = other.rings;
    }
    return *this;
  }
  
  std::vector<std::vector<double>> centers() const {
    std::vector<std::vector<double>> centers;
    for (const auto &ring : rings) {
      centers.push_back({ring.center.x, ring.center.y, ring.center.z});
    }
    return centers;
  }

  // FIX: probably not needed? 
  RDGeom::POINT3D_VECT centers_rkdit() const {
    RDGeom::POINT3D_VECT centers;
    for (const auto &ring : rings) {
      centers.push_back(ring.center);
    }
    return centers;
  }

  std::vector<std::vector<double>> norm1() const {
    std::vector<std::vector<double>> normals;
    for (const auto &ring : rings) {
      normals.push_back({ring.norm1.x, ring.norm1.y, ring.norm1.z});
    }
    return normals;
  }

  std::vector<std::vector<double>> norm2() const {
    std::vector<std::vector<double>> normals;
    for (const auto &ring : rings) {
      normals.push_back({ring.norm2.x, ring.norm2.y, ring.norm2.z});
    }
    return normals;
  }
};

class Rings {
public:
  struct ResId {
    std::string chain;
    int res_num;
    std::string res_name;

    bool operator==(const ResId &other) const {
      return chain == other.chain && res_num == other.res_num &&
             res_name == other.res_name;
    }
  };

  struct ResIdHash {
    std::size_t operator()(const ResId &id) const {
      return std::hash<std::string>()(id.chain) ^
             (std::hash<int>()(id.res_num) << 1) ^
             (std::hash<std::string>()(id.res_name) << 2);
    }
  };

  // FIX: benchmark performance for many ring systems? 
  using RingMap = std::unordered_map<ResId, RingData, ResIdHash>;
  using AtomOrderMap = std::unordered_map<std::string, std::vector<std::string>>;

  Rings() { init_atom_order_table(); }

  void add_ring_atom(const RDKit::Atom *atom) {
    auto *info =
        static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
    if (!info)
      return;

    ResId res_id{info->getChainId(), info->getResidueNumber(),
                 info->getResidueName()};
    const bool is_trp = info->getResidueName() == "TRP";
    if (is_trp) {
      res_id.res_name = "TRP5";
      rings_[res_id].atom_ids.push_back(atom->getIdx());
      res_id.res_name = "TRP6";
    }
    rings_[res_id].atom_ids.push_back(atom->getIdx());
  }

  void process_rings(const RDKit::ROMol &mol) {
    const RDKit::Conformer &conf = mol.getConformer();
    for (auto &[res_id, ring_data] : rings_) {
      if (!ring_data.atom_ids.empty()) {
        order_atoms(res_id.res_name, ring_data.atom_ids, mol);
        find_center_and_normal(conf, ring_data.atom_ids, ring_data.center,
                               ring_data.norm1, ring_data.norm2);
      }
    }
  }

  const RingDataVec get_rings_vector() {
    RingDataVec ring_data;
    for (auto &[res_id, data] : rings_) {
      ring_data.rings.push_back(data);
    }
    return ring_data;
  }

  RingData get_ring(const ResId &res_id) { return rings_[res_id]; }
  const RingMap &get_rings() const { return rings_; }

private:
  RingMap rings_;
  AtomOrderMap atom_orders_;

  void init_atom_order_table() {
    atom_orders_["PRO"] = {"CG", "CB", "CA", "N", "CD"};
    atom_orders_["TYR"] = {"CE2", "CZ", "CE1", "CD1", "CG", "CD2"};
    atom_orders_["HIS"] = {"CD2", "NE2", "CE1", "ND1", "CG"};
    atom_orders_["PHE"] = {"CE1", "CZ", "CE2", "CD2", "CG", "CD1"};
    atom_orders_["TRP5"] = {"NE1", "CE2", "CD2", "CG", "CD1"};
    atom_orders_["TRP6"] = {"CZ2", "CH2", "CZ3", "CE3", "CD2", "CE2"};
  }

  void order_atoms(const std::string &res_name, std::vector<int> &atom_ids,
                   const RDKit::ROMol &mol) {
    if (atom_orders_.find(res_name) == atom_orders_.end())
      return;

    const auto &expected_order = atom_orders_[res_name];
    std::vector<int> ordered_ids;
    ordered_ids.reserve(expected_order.size());

    for (const auto &exptected_name : expected_order) {
      auto it = std::find_if(atom_ids.begin(), atom_ids.end(), [&](int idx) {
        const RDKit::Atom *atom = mol.getAtomWithIdx(idx);
        auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(
            atom->getMonomerInfo());
        return info && info->getName() == exptected_name;
      });
      if (it != atom_ids.end()) {
        ordered_ids.push_back(*it);
      }
    }
    atom_ids = std::move(ordered_ids);
    // NOTE: Data can have missing atoms, which is currently not considered.
    // This probably only matters for PRO since it's not aromatic and hence not
    // planar.
  }

public:
  // FIX: is norm2 needed? 
  static void find_center_and_normal(const RDKit::Conformer &conf,
                                     const std::vector<int> &ring,
                                     RDGeom::Point3D &center,
                                     RDGeom::Point3D &norm1,
                                     RDGeom::Point3D &norm2) {

    center = std::accumulate(ring.begin(), ring.end(), center,
                             [&conf](const RDGeom::Point3D &sum, int idx) {
                               return sum + conf.getAtomPos(idx);
                             });
    center /= static_cast<double>(ring.size());

    for (size_t i = 0; i < ring.size(); ++i) {
      RDGeom::Point3D v1 = conf.getAtomPos(ring[i]) - center;
      RDGeom::Point3D v2 =
          conf.getAtomPos(ring[(i + 1) % ring.size()]) - center;
      norm1 += v1.crossProduct(v2);
    }
    norm1 /= static_cast<double>(ring.size());
    norm1.normalize();
    norm2 = -norm1;
  }
};

} // namespace lahuta

#endif // LAHUTA_RINGS_HPP
