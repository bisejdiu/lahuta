#ifndef LAHUTA_DISTANCES_HPP
#define LAHUTA_DISTANCES_HPP

#include "contacts/groups.hpp"
#include "contacts/interactions.hpp"
#include "kdtree/KDTree.hpp"
#include <GraphMol/Atom.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/RWMol.h>
#include <set>
#include <vector>

namespace lahuta {

// check a similar definition in contacts/contacts.hpp
struct _AtomData_ {
  const RDKit::Atom *atom;
  point_t position; // atom position
  size_t feature_id;
  bool is_type_a;
};

class KDTreeSearch {
  struct Pairing {
    std::pair<size_t, size_t> indices;
    double distance_sq;
  };

public:
  KDTreeSearch(const std::vector<_AtomData_> &all_atoms) {
    for (const auto &atom_data : all_atoms) {
      positions_.push_back(atom_data.position);
      atom_data_.push_back(atom_data);
    }
    tree_ = std::make_unique<KDTree>(positions_);
  }

  void find_interactions(
      double distance_threshold, InteractionBase &handler, const std::vector<const Feature *> &feature_map) {
    std::set<std::pair<size_t, size_t>> interactions;
    for (size_t i = 0; i < atom_data_.size(); ++i) {
      const auto &atom_data = atom_data_[i];
      if (!atom_data.is_type_a)
        continue;

      indexArr neighbors = tree_->neighborhood_indices(atom_data.position, distance_threshold);
      for (size_t idx : neighbors) {
        if (idx == i)
          continue;

        const auto &neighbor_atom = atom_data_[idx];
        if (!neighbor_atom.is_type_a) {
          interactions.emplace(atom_data.feature_id, neighbor_atom.feature_id);
        }
      }
    }

    for (const auto &interaction : interactions) {
      const auto *feature1 = feature_map[interaction.first];
      const auto *feature2 = feature_map[interaction.second];
      handler.handle(*feature1, *feature2);
    }
  }

private:
  std::vector<point_t> positions_;
  std::vector<_AtomData_> atom_data_;
  std::unique_ptr<KDTree> tree_;
};

class BruteForceSearch {
public:
  BruteForceSearch(
      const std::vector<const Feature *> &features_type_a,
      const std::vector<const Feature *> &features_type_b)
      : features_type_a_(features_type_a), features_type_b_(features_type_b) {}

  void find_interactions(double distance_threshold, InteractionBase &handler) {
    double dist_sq = distance_threshold * distance_threshold;
    for (const auto *feature_a : features_type_a_) {
      for (const auto *feature_b : features_type_b_) {
        if (are_features_within_dist_sq(*feature_a, *feature_b, dist_sq)) {
          handler.handle(*feature_a, *feature_b);
        }
      }
    }
  }

private:
  bool are_features_within_dist_sq(const Feature &info_a, const Feature &info_b, double dist_sq) {
    size_t num_checks = info_a.members.size() * info_b.members.size();
    auto &conformer = info_a.members.front()->getOwningMol().getConformer();

    static auto is_within_dist_sq = [&](const RDKit::Atom *atom_a, const RDKit::Atom *atom_b) {
      auto pos_a = conformer.getAtomPos(atom_a->getIdx());
      auto pos_b = conformer.getAtomPos(atom_b->getIdx());
      return (pos_a - pos_b).lengthSq() < dist_sq;
    };

    if (num_checks == 1) {
      return is_within_dist_sq(info_a.members.front(), info_b.members.front());
    }

    if (num_checks == 2) {
      if (is_within_dist_sq(info_a.members.front(), info_b.members.front())) {
        return true;
      }
      const auto *atom_a2 = (info_a.members.size() > 1) ? info_a.members.back() : info_a.members.front();
      const auto *atom_b2 = (info_b.members.size() > 1) ? info_b.members.back() : info_b.members.front();
      return is_within_dist_sq(atom_a2, atom_b2);
    }

    for (const auto *atom_a : info_a.members) {
      for (const auto *atom_b : info_b.members) {
        if (is_within_dist_sq(atom_a, atom_b)) {
          return true;
        }
      }
    }

    return false;
  }

private:
  const std::vector<const Feature *> &features_type_a_;
  const std::vector<const Feature *> &features_type_b_;
};

class SimplePairFeatures {
public:
  SimplePairFeatures(const FeatureVec &group_features) : group_features_(group_features.features) {}

  void process_features(AtomType type1, AtomType type2) {
    features_a = get_features(group_features_, type1);
    features_b = get_features(group_features_, type2);
  }

  Contacts find_interactions(InteractionType interaction_type, double distance_threshold) {
    auto handler = create_handler(interaction_type);
    compute_interactions(distance_threshold, *handler);
    return handler->release();
  }

private:
  void compute_interactions(double distance_threshold, InteractionBase &handler) {
    BruteForceSearch brute_force_finder(features_a, features_b);
    brute_force_finder.find_interactions(distance_threshold, handler);
  }

  static std::unique_ptr<InteractionBase> create_handler(InteractionType type) {
    switch (type) {
      case InteractionType::Ionic:
        return std::make_unique<IonicInteraction>();
      default:
        throw std::invalid_argument("Unsupported interaction type");
    }
  }

  const std::vector<Feature> &group_features_;
  std::vector<const Feature *> features_a;
  std::vector<const Feature *> features_b;
};

class PairFeatures {
public:
  static const size_t MAX_DISTANCE_THRESHOLD = 1;

  PairFeatures(const std::vector<Feature> &group_features) : group_features_(group_features) {}

  void process_features(AtomType type1, AtomType type2) {
    auto feature_a = get_features(group_features_, type1);
    auto feature_b = get_features(group_features_, type2);
    if (should_brute_force(feature_a, feature_b)) {
      features_type_a_ = feature_a;
      features_type_b_ = feature_b;
      return;
    }
    prepare_atom_data(feature_a, true);
    prepare_atom_data(feature_b, false);
  }

  Contacts find_interactions(InteractionType interaction_type, double distance_threshold) {
    auto handler = create_handler(interaction_type);
    compute_interactions(distance_threshold, *handler);
    return handler->release();
  }

  bool should_brute_force(
      const std::vector<const Feature *> &features_a, const std::vector<const Feature *> &features_b) {
    /*return false;*/
    size_t num_distance_computations = features_a.size() * features_b.size();
    return num_distance_computations <= MAX_DISTANCE_THRESHOLD;
  }

private:
  void prepare_atom_data(const std::vector<const Feature *> &features, bool is_type_a) {
    for (const auto *feature : features) {
      feature_map_.push_back(feature);
      add_atoms_from_feature(feature, is_type_a, feature_map_.size() - 1);
      if (is_type_a) {
        features_type_a_.push_back(feature);
      } else {
        features_type_b_.push_back(feature);
      }
    }
  }

  void add_atoms_from_feature(const Feature *feature, bool is_type_a, size_t feature_id) {
    for (const RDKit::Atom *atom : feature->members) {
      auto position = get_atom_position(atom);
      all_atoms_.push_back({atom, position, feature_id, is_type_a});
    }
  }

  void compute_interactions(double distance_threshold, InteractionBase &handler) {
    if (!should_brute_force(features_type_a_, features_type_b_)) {
      KDTreeSearch kd_tree_finder(all_atoms_);
      kd_tree_finder.find_interactions(distance_threshold, handler, feature_map_);
    } else {
      BruteForceSearch brute_force_finder(features_type_a_, features_type_b_);
      brute_force_finder.find_interactions(distance_threshold, handler);
    }
  }

  static std::unique_ptr<InteractionBase> create_handler(InteractionType type) {
    switch (type) {
      case InteractionType::Ionic:
        return std::make_unique<IonicInteraction>();
      default:
        throw std::invalid_argument("Unsupported interaction type");
    }
  }

  static point_t get_atom_position(const RDKit::Atom *atom) {
    const auto &conformer = atom->getOwningMol().getConformer();
    auto pos = conformer.getAtomPos(atom->getIdx());
    return {pos.x, pos.y, pos.z};
  }

  const std::vector<Feature> &group_features_;
  std::vector<_AtomData_> all_atoms_;
  std::vector<const Feature *> feature_map_;
  std::vector<const Feature *> features_type_a_;
  std::vector<const Feature *> features_type_b_;
};

} // namespace lahuta

#endif // LAHUTA_DISTANCES_HPP
