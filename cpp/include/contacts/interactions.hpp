#ifndef LAHUTA_INTERACTIONS_HPP
#define LAHUTA_INTERACTIONS_HPP

#include "nn.hpp"

namespace lahuta {

class Luni;

inline void log_feature_atoms(const std::string &feature_type, const Feature &feature) {
  std::cout << "i. " << feature_type << " Charge Feature Atoms:" << std::endl;
  std::cout << "Number of atoms: " << feature.members.size() << std::endl;
  for (const auto *atom : feature.members) {
    unsigned int atom_index = atom->getIdx();
    auto *res_info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
    std::string atom_name = res_info->getName();
    std::string residue_name = res_info->getResidueName();
    std::string chain_id = res_info->getChainId();
    auto residue_number = res_info->getResidueNumber();

    std::cout << "i.  Atom Index: " << atom_index << ", Atom Name: " << atom_name
              << ", Residue: " << residue_name << " " << chain_id << residue_number << std::endl;
  }
}

class InteractionBase {
  Contacts contacts;

public:
  virtual void handle(const Feature &feature1, const Feature &feature2) = 0;
  virtual ~InteractionBase() = default;
  void add_interaction(const Contact &interaction) { contacts.add(interaction); }
  Contacts release() { return std::move(contacts); }
};

class IonicInteraction : public InteractionBase {
public:
  void handle(const Feature &feature1, const Feature &feature2) override {
    EntityID entity1 = make_entity_id(EntityType::Group, feature1.get_id());
    EntityID entity2 = make_entity_id(EntityType::Group, feature2.get_id());

    auto center_dist_sq = FastNS::dist_sq(feature1.center, feature2.center);
    add_interaction(Contact(entity1, entity2, center_dist_sq, InteractionType::Ionic));
  }
};

struct InteractionOptions {
  float distance_cutoff;
};

class Interactions {
public:
  // TODO: InteractionOptions should have default values, and each of the finder functions
  // should take another custom options struct
  Interactions(Luni *luni, const std::vector<Feature> &group_features, InteractionOptions opts)
      : luni_(luni), group_features_(group_features), opts_(opts) {}

  [[nodiscard]] Contacts find_ionic_interactions();
  [[nodiscard]] Contacts find_hbond_interactions();

private:
  Luni *luni_;
  const std::vector<Feature> &group_features_;
  InteractionOptions opts_;
};

} // namespace lahuta

#endif // LAHUTA_INTERACTIONS_HPP
