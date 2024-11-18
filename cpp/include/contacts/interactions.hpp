#ifndef LAHUTA_INTERACTIONS_HPP
#define LAHUTA_INTERACTIONS_HPP

#include "nn.hpp"

namespace lahuta {

class Luni;

class InteractionBase {
  Contacts contacts;

public:
  virtual void handle(const GroupEntity &feature1, const GroupEntity &feature2) = 0;
  virtual ~InteractionBase() = default;
  void add_interaction(const Contact &interaction) { contacts.add(interaction); }
  Contacts release() { return std::move(contacts); }
};

inline double compute_dist_sq(const RDGeom::Point3D &p1, const RDGeom::Point3D &p2) {
  double p1_c[3] = {p1.x, p1.y, p1.z};
  double p2_c[3] = {p2.x, p2.y, p2.z};
  return FastNS::dist_sq(p1_c, p2_c);
}

/*class IonicInteraction : public InteractionBase {*/
/*public:*/
/*  void handle(const Feature &feature1, const Feature &feature2) override {*/
/*    EntityID entity1 = make_entity_id(EntityType::Group, feature1.get_id());*/
/*    EntityID entity2 = make_entity_id(EntityType::Group, feature2.get_id());*/
/**/
/*    auto center_dist_sq = compute_dist_sq(feature1.center, feature2.center);*/
/*    add_interaction(Contact(entity1, entity2, center_dist_sq, InteractionType::Ionic));*/
/*  }*/
/*};*/

struct InteractionOptions {
  float distance_cutoff;
};

class Interactions {
public:
  Interactions(const Luni &luni, InteractionOptions opts) : luni_(luni), opts_(opts) {}

  [[nodiscard]] Contacts hbond();
  [[nodiscard]] Contacts weak_hbond();
  [[nodiscard]] Contacts hydrophobic();
  [[nodiscard]] Contacts halogen();
  [[nodiscard]] Contacts ionic();
  [[nodiscard]] Contacts metalic();
  [[nodiscard]] Contacts cationpi();
  [[nodiscard]] Contacts pistacking();

private:
  const Luni &luni_;
  /*const std::vector<Feature> *group_features_;*/
  InteractionOptions opts_{5.0};
};

} // namespace lahuta

#endif // LAHUTA_INTERACTIONS_HPP
