#ifndef LAHUTA_ENTITIES_HPP
#define LAHUTA_ENTITIES_HPP

#include "Geometry/point.h"
#include "GraphMol/RWMol.h"
#include "atom_types.hpp"

// clang-format off
namespace lahuta {

class Luni;
using FeatureTypeCheckFunc = std::function<bool(const AtomType &, const AtomType &)>;

enum class RingGroup { NONE, AROMATIC, ALIPHATIC };
enum class RingType { None };
enum class FeatureGroup {
  None = 0,
  QuaternaryAmine = 1,
  TertiaryAmine = 2,
  Sulfonium = 3,
  SulfonicAcid = 4,
  Sulfate = 5,
  Phosphate = 6,
  Halocarbon = 7,
  Guanidine = 8,
  Acetamidine = 9,
  Carboxylate = 10
};

struct Entity {
  using GetIdFn = size_t (*)(const void*);
  using GetCenterFn = const RDGeom::Point3D& (*)(const void*);

  template <typename T>
  Entity(const T& t)
    : get_id_fn    ([](const void* obj) -> size_t                 { return static_cast<const T*>(obj)->get_id(); }),
      get_center_fn([](const void* obj) -> const RDGeom::Point3D& { return static_cast<const T*>(obj)->get_center(); }),
      obj(&t) {}

  size_t get_id() const { return get_id_fn(obj); }
  const RDGeom::Point3D& get_center() const { return get_center_fn(obj); }
private:
  GetCenterFn get_center_fn;
  GetIdFn get_id_fn;
  const void* obj;
};

// FIX: initialization argument order is inconsistent
struct AtomEntity {
  AtomType type;
  const RDKit::Atom *atom;
  std::vector<const RDKit::Atom*> atoms;
  const RDGeom::Point3D *center;

  size_t get_id() const { return id; }
  const RDGeom::Point3D& get_center() const { return *center; }
  const RDKit::Atom *get_data() const { return atom; }
  bool has_atom(const RDKit::Atom *atom_) const { return std::find(atoms.begin(), atoms.end(), atom_) != atoms.end(); }


  explicit AtomEntity(
    AtomType type_,
    std::vector<const RDKit::Atom*> atoms_,
    const RDGeom::Point3D *pos_, size_t id_
  )
    :
    type(type_),
    atoms(atoms_),
    id(id_), // FIX: check if id (get_id) is used or even needed
    center(pos_),
    atom(atoms_.front()) {
  }

private:
  size_t id;
};

struct RingEntity {
  std::vector<const RDKit::Atom *> atoms;
  RDGeom::Point3D center;
  RDGeom::Point3D norm; // FIX: should be renamed to normal rather than norm

  size_t get_id() const { return id; }
  const RDGeom::Point3D& get_center() const { return center; }
  std::vector<const RDKit::Atom *> get_data() const { return atoms; }
  bool has_atom(const RDKit::Atom *atom) const { return std::find(atoms.begin(), atoms.end(), atom) != atoms.end(); }
  // FIX: add has_atom overload that accepts atom index

public:
  explicit RingEntity(
    RDGeom::Point3D center_,
    RDGeom::Point3D norm_,
    std::vector<const RDKit::Atom *> atoms_,
    size_t id_
  )
    : center(center_), norm(norm_), atoms(atoms_), id(id_) {}

private:
  size_t id;
};

// FIX: GroupEntity does not handle the `id` field the same as AtomEntity and RingEntity
struct GroupEntity {
  AtomType type;
  FeatureGroup group; // FIX: I think we should store the group (or post-compute it)
  std::vector<const RDKit::Atom *> atoms;
  RDGeom::Point3D center;

  size_t get_id() const { return id; }
  const RDGeom::Point3D& get_center() const { return center; }
  std::vector<const RDKit::Atom *> get_data() const { return atoms; }
  bool has_atom(const RDKit::Atom *atom) const { return std::find(atoms.begin(), atoms.end(), atom) != atoms.end(); }

  void set_id(size_t id_) { id = id_; }

  GroupEntity(
    AtomType type_,
    FeatureGroup group_, 
    std::vector<const RDKit::Atom *> members_, 
    RDGeom::Point3D center_
  )
    : type(type_), group(group_), atoms(std::move(members_)), center(center_) {}

private:
  size_t id;
};

template <typename EntityType>
struct EntityCollection {
    using value_type = EntityType;

    EntityCollection() = default;
    EntityCollection(std::vector<EntityType> entities) : data(std::move(entities)) {}
    virtual ~EntityCollection() = default;

    void add_data(const EntityType entity) { data.push_back(entity); }
    const std::vector<EntityType>& get_data() const { return data; }
    /*virtual const EntityVecBase& filter_data(const Luni* luni, AtomType type, FeatureTypeCheckFunc check_func) const = 0;*/

    EntityType& operator[](size_t index) { return data[index]; }
    const EntityType& operator[](size_t index) const { return data[index]; }

    int size() const { return static_cast<int>(data.size()); }
    void reserve(size_t size) { data.reserve(size); }

    // FIX: positions returns a copy of the points but performance impact is negligible 
    const RDGeom::POINT3D_VECT positions() const;
    const std::vector<size_t> atom_ids() const;
    const std::vector<const RDKit::Atom *> atoms() const;

    typename std::vector<EntityType>::const_iterator begin() const { return data.begin(); }
    typename std::vector<EntityType>::const_iterator end() const { return data.end(); }
    typename std::vector<EntityType>::const_iterator cbegin() const { return data.cbegin(); }
    typename std::vector<EntityType>::const_iterator cend() const { return data.cend(); }


    std::vector<EntityType> data;
};

class AtomEntityCollection : public EntityCollection<AtomEntity> {
public:
    AtomEntityCollection() = default;

    using EntityCollection::add_data;
    void add_data(const RDKit::RWMol& mol, const RDKit::Atom* atom, AtomType type);
    static const AtomEntityCollection filter(const Luni* luni, AtomType type, FeatureTypeCheckFunc check_func = AtomTypeFlags::has_any);

    const RDGeom::POINT3D_VECT positions() const;
};

class RingEntityCollection : public EntityCollection<RingEntity> {
public:
    RingEntityCollection() = default;

    using EntityCollection::add_data;
    void add_data(const RDKit::RWMol& mol, const std::vector<int>& ring, int id);

    // FIX: move to ring props?
    double compute_angle(const RingEntity& rd, const std::vector<double>& point) const;
    std::vector<double> compute_angles(const std::vector<int>& ring_indices, const std::vector<std::vector<double>>& points) const;
};

class GroupEntityCollection : public EntityCollection<GroupEntity> {
public:
    GroupEntityCollection() = default;
    GroupEntityCollection(std::vector<GroupEntity> features) : EntityCollection(std::move(features)) {} 

    using EntityCollection::add_data;
    void add_data(AtomType type, FeatureGroup group, const std::vector<const RDKit::Atom*>& members);

    static const GroupEntityCollection filter(const Luni* luni, AtomType type, FeatureTypeCheckFunc check_func = AtomTypeFlags::has_any);
};

/*std::vector<const RDKit::Atom *> get_atom_types(const Luni *luni, AtomType type);*/

// clang-format off
template <typename EntityType> const std::vector<size_t>
EntityCollection<EntityType>::atom_ids() const {
  std::vector<size_t> ids;
  ids.reserve(data.size());
  for (const auto &ring : data) {
    for (const auto *atom : ring.atoms) {
      ids.push_back(atom->getIdx());
    }
  }
  return ids;
}

template <typename EntityType> const RDGeom::POINT3D_VECT
EntityCollection<EntityType>::positions() const {
  RDGeom::POINT3D_VECT pos;
  pos.reserve(data.size());
  for (const auto &feature : data) {
    pos.push_back(feature.center);
  }
  return pos;
}

template <typename EntityType>
const std::vector<const RDKit::Atom *>
EntityCollection<EntityType>::atoms() const {
  std::vector<const RDKit::Atom *> atoms_vec;
  atoms_vec.reserve(data.size());
  for (const auto &d : data) {
    atoms_vec.insert(atoms_vec.end(), d.atoms.begin(), d.atoms.end());
  }
  return atoms_vec;
}

} // namespace lahuta

#endif // LAHUTA_ENTITIES_HPP
