#ifndef LAHUTA_CONTACTS_ARPEGGIO_HPP
#define LAHUTA_CONTACTS_ARPEGGIO_HPP

#include "contacts/utils.hpp"
#include "entities/find_contacts.hpp"
#include "elements.hpp"
#include <logging.hpp>
#include "contacts/pistacking.hpp"
#include "contacts/arpeggio_helpers.hpp"

namespace lahuta {

struct IonicParams {
  float distance_max = 4.0;
};
inline ContactSet find_ionic(const Topology& topology, const IonicParams& params = IonicParams{}) {
  return find_contacts(
    topology,
    [](const AtomRec &r) { return (r.type & AtomType::PositiveCharge) == AtomType::PositiveCharge; },
    [](const AtomRec &r) { return (r.type & AtomType::NegativeCharge) == AtomType::NegativeCharge; },
    {params.distance_max, 0.5, 0.5},
    [&topology, &params](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist) -> InteractionType {
      const auto &mol = topology.molecule();
      const auto *atom_a = mol.getAtomWithIdx(topology.atom(rec_idx_a).idx);
      const auto *atom_b = mol.getAtomWithIdx(topology.atom(rec_idx_b).idx);
      if (are_residueids_close(mol, *atom_a, *atom_b, 1)) return InteractionType::None;
      return InteractionType::Ionic;
    }
  );
}

struct AromaticParams {
  float distance_max = 4.0;
};

inline ContactSet find_aromatic(const Topology &topology, const AromaticParams& params = AromaticParams{}) {
  return find_contacts(
    topology,
    [](const AtomRec& rec) { return (rec.type & AtomType::Aromatic) == AtomType::Aromatic; },
    {params.distance_max, 0.1, 1},
    [&topology, &params](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist) -> InteractionType {

      const auto& ring_rec_a = topology.atom(rec_idx_a);
      const auto& ring_rec_b = topology.atom(rec_idx_b);

      const auto& mol = topology.molecule();
      const auto* r_atom_a = mol.getAtomWithIdx(ring_rec_a.idx);
      const auto* r_atom_b = mol.getAtomWithIdx(ring_rec_b.idx);

      if (are_residueids_close(mol, *r_atom_a, *r_atom_b, 1)) return InteractionType::None;
      return InteractionType::PiStackingP;
    }
  );
}

struct CarbonylParams {
  float distance_max = 3.6;
};

inline ContactSet find_carbonyl(const Topology &topology, const CarbonylParams& params = CarbonylParams{}) {
  return find_contacts(
    topology,
    [](const AtomRec& rec) { return (rec.type & AtomType::CarbonylCarbon) == AtomType::CarbonylCarbon; },
    [](const AtomRec& rec) { return (rec.type & AtomType::CarbonylOxygen) == AtomType::CarbonylOxygen; },
    {params.distance_max, 0.1, 1},
    [&topology, &params](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist) -> InteractionType {

      const auto& ring_rec_a = topology.atom(rec_idx_a);
      const auto& ring_rec_b = topology.atom(rec_idx_b);

      const auto& mol = topology.molecule();
      const auto* r_atom_a = mol.getAtomWithIdx(ring_rec_a.idx);
      const auto* r_atom_b = mol.getAtomWithIdx(ring_rec_b.idx);

      if (are_residueids_close(mol, *r_atom_a, *r_atom_b, 1)) return InteractionType::None;
      return InteractionType::PiStackingP;
    }
  );
}


struct HydrophobicParams {
  float distance_max = 4.5;
};

inline ContactSet find_hydrophobic(const Topology& topology, const HydrophobicParams& params = HydrophobicParams{}) {
  return find_contacts(
      topology,
      [](const AtomRec &rec) { return (rec.type & AtomType::Hydrophobic) == AtomType::Hydrophobic; },
      {params.distance_max, 0.4},
      [&topology, &params](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist) -> InteractionType {
        const auto &mol = topology.molecule();

        const auto atom_a_rec = topology.atom(rec_idx_a);
        const auto atom_b_rec = topology.atom(rec_idx_b);
        const auto *atom_a = mol.getAtomWithIdx(atom_a_rec.idx);
        const auto *atom_b = mol.getAtomWithIdx(atom_b_rec.idx);

        if (are_residueids_close(mol, *atom_a, *atom_b, 1)) return InteractionType::None;
        return InteractionType::Hydrophobic;
      }
  );
}

struct VanDerWaalsParams {
  float distance_max = 4.5;
};

inline ContactSet find_vdw(const Topology &topology, const VanDerWaalsParams& params = VanDerWaalsParams{}) {
  return find_contacts(
    topology,
    [](const AtomRec& rec) { return true; },
    [](const AtomRec& rec) { return true; },
    {params.distance_max, 0.1, 1},
    [&topology, &params](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist_sq) -> InteractionType {

      static const auto vdw_comp_factor = 0.1;
      static const bool remove_clashes = true;

      const auto& ring_rec_a = topology.atom(rec_idx_a);
      const auto& ring_rec_b = topology.atom(rec_idx_b);

      const auto& mol = topology.molecule();
      const auto* r_atom_a = mol.getAtomWithIdx(ring_rec_a.idx);
      const auto* r_atom_b = mol.getAtomWithIdx(ring_rec_b.idx);

      auto vdw_a = vdw_radius(static_cast<Element>(r_atom_a->getAtomicNum()));
      auto vdw_b = vdw_radius(static_cast<Element>(r_atom_b->getAtomicNum()));

      float sum_vdw = vdw_a + vdw_b;
      float max_sq = (sum_vdw + vdw_comp_factor) * (sum_vdw + vdw_comp_factor);

      if (dist_sq > max_sq || dist_sq < (sum_vdw * sum_vdw))  return InteractionType::None;
      if (are_residueids_close(mol, *r_atom_a, *r_atom_b, 1)) return InteractionType::None;

      return InteractionType::PiStackingP;
    }
  );
}

struct HbondParams {
  float distance_max = 4.5;
};

inline ContactSet find_hbond(const Topology& topology, const HbondParams& params = HbondParams{}) {
    return find_contacts(
        topology,
        [](const AtomRec &rec) { return (rec.type & AtomType::HbondDonor) == AtomType::HbondDonor; },
        [](const AtomRec &rec) { return (rec.type & AtomType::HbondAcceptor) == AtomType::HbondAcceptor; },
        {params.distance_max, 0.7},
        [&topology, &params](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist_sq) -> InteractionType {
          const auto& mol = topology.molecule();
          const auto& conformer = topology.conformer();

          const auto  donor_rec     = topology.atom(rec_idx_a);
          const auto  acceptor_rec  = topology.atom(rec_idx_b);
          const auto* donor_atom    = mol.getAtomWithIdx(donor_rec.idx);
          const auto* acceptor_atom = mol.getAtomWithIdx(acceptor_rec.idx);

          if (are_residueids_close(mol, *donor_atom, *acceptor_atom, 1)) return InteractionType::None;

          if (!passes_hbond_distance_filter(mol, conformer, donor_atom, acceptor_atom, 0.1)) return InteractionType::None;
          if (!passes_hbond_angle_filter(mol, conformer, donor_atom, acceptor_atom, 1.57, false)) return InteractionType::None; // 90 degrees

          return InteractionType::HydrogenBond;
        }
    );
}

struct WeakHbondParams {
  float distance_max = 4.5;
};

inline ContactSet find_weak_hbond(const Topology& topology, const WeakHbondParams& params = WeakHbondParams{}) {
    return find_contacts(
        topology,
        [](const AtomRec &rec) { return (rec.type & AtomType::WeakHbondDonor) == AtomType::WeakHbondDonor; },
        [](const AtomRec &rec) { return (rec.type & AtomType::HbondAcceptor) == AtomType::HbondAcceptor; },
        {params.distance_max, 0.5},
        [&topology, &params](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist_sq) -> InteractionType {
          const auto& mol = topology.molecule();
          const auto& conformer = topology.conformer();

          const auto  donor_rec    = topology.atom(rec_idx_a);
          const auto  acceptor_rec = topology.atom(rec_idx_b);
          const auto* donor_atom    = mol.getAtomWithIdx(donor_rec.idx);
          const auto* acceptor_atom = mol.getAtomWithIdx(acceptor_rec.idx);

          if (are_residueids_close(mol, *donor_atom, *acceptor_atom, 1)) return InteractionType::None;

          if (!passes_hbond_distance_filter(mol, conformer, donor_atom, acceptor_atom, 0.1))     return InteractionType::None;
          if (!passes_hbond_angle_filter(mol, conformer, donor_atom, acceptor_atom, 2.27, true)) return InteractionType::None; // 130 degrees

          return InteractionType::WeakHydrogenBond;
        }
    );
}

struct PolarHbondParams {
  float distance_max = 3.5;
};

inline ContactSet find_polar_hbond(const Topology &topology, const PolarHbondParams& params = PolarHbondParams{}) {
  return find_contacts(
    topology,
    [](const AtomRec& rec) { return (rec.type & AtomType::HbondDonor) == AtomType::HbondDonor; },
    [](const AtomRec& rec) { return (rec.type & AtomType::HbondAcceptor) == AtomType::HbondAcceptor; },
    {params.distance_max, 0.5},
    [&topology, &params](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist) -> InteractionType {

      const auto& ring_rec_a = topology.atom(rec_idx_a);
      const auto& ring_rec_b = topology.atom(rec_idx_b);

      const auto& mol = topology.molecule();
      const auto* r_atom_a = mol.getAtomWithIdx(ring_rec_a.idx);
      const auto* r_atom_b = mol.getAtomWithIdx(ring_rec_b.idx);

      if (are_residueids_close(mol, *r_atom_a, *r_atom_b, 1)) return InteractionType::None;
      return InteractionType::PiStackingP;
    }
  );
}

struct WeakPolarHbondParams {
  float distance_max = 3.5;
};

inline ContactSet find_weak_polar_hbond(const Topology &topology, const WeakPolarHbondParams& params = WeakPolarHbondParams{}) {
  return find_contacts(
    topology,
    [](const AtomRec& rec) { return (rec.type & AtomType::WeakHbondDonor) == AtomType::WeakHbondDonor; },
    [](const AtomRec& rec) { return (rec.type & AtomType::HbondAcceptor)  == AtomType::HbondAcceptor; },
    {params.distance_max, 0.5},
    [&topology, &params](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist) -> InteractionType {

      const auto& ring_rec_a = topology.atom(rec_idx_a);
      const auto& ring_rec_b = topology.atom(rec_idx_b);

      const auto& mol = topology.molecule();
      const auto* r_atom_a = mol.getAtomWithIdx(ring_rec_a.idx);
      const auto* r_atom_b = mol.getAtomWithIdx(ring_rec_b.idx);

      if (are_residueids_close(mol, *r_atom_a, *r_atom_b, 1)) return InteractionType::None;

      return InteractionType::PiStackingP;
    }
  );
}

struct DonorPiParams {
  float distance_max = 4.5;
};

inline ContactSet donor_pi(const Topology& topology, const DonorPiParams& params = DonorPiParams{}) {
  return find_contacts(
      topology,
      [](const AtomRec &rec) { return (rec.type & AtomType::HbondDonor) == AtomType::HbondDonor; },
      [](const RingRec &rec) { return rec.aromatic; },
      {params.distance_max, 0.4},
      [&topology, &params](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist_sq) -> InteractionType {
        const auto &mol = topology.molecule();

        const auto atom_a_rec = topology.atom(rec_idx_a);
        const auto atom_b_rec = topology.ring(rec_idx_b);
        const auto *atom_a = mol.getAtomWithIdx(atom_a_rec.idx);
        const auto *atom_b = mol.getAtomWithIdx(atom_b_rec.atoms[0]);

        auto angle = compute_angle(atom_b_rec, topology.conformer().getAtomPos(atom_a->getIdx()));

        const auto angle_cutoff = 30.0;
        if (atom_a->getIsAromatic()) return InteractionType::None;
        if (!passes_angle_filter(angle, angle_cutoff)) return InteractionType::None;
        if (are_residueids_close(mol, *atom_a, *atom_b, 1)) return InteractionType::None;

        return InteractionType::Hydrophobic;
      }
  );
}

struct SulphurPiParams {
  float distance_max = 6.0;
};

inline ContactSet sulphur_pi(const Topology& topology, const SulphurPiParams& params = SulphurPiParams{}) {
  return find_contacts(
      topology,
      [](const AtomRec &rec) { return (rec.atom.getAtomicNum() == Element::S); },
      [](const RingRec &rec) { return rec.aromatic; },
      {params.distance_max, 0.4},
      [&topology, &params](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist_sq) -> InteractionType {
        const auto &mol = topology.molecule();

        const auto atom_a_rec = topology.atom(rec_idx_a);
        const auto atom_b_rec = topology.ring(rec_idx_b);
        const auto *atom_a = mol.getAtomWithIdx(atom_a_rec.idx);
        const auto *atom_b = mol.getAtomWithIdx(atom_b_rec.atoms[0]);

        auto res_a_info = static_cast<const RDKit::AtomPDBResidueInfo*>(atom_a->getMonomerInfo());

        // this can be added as part of the predicate
        if (!res_a_info || res_a_info->getResidueName() != "MET") return InteractionType::None;
        if (are_residueids_close(mol, *atom_a, *atom_b, 1)) return InteractionType::None;

        return InteractionType::Hydrophobic;
      }
  );
}

struct CarbonPiParams {
  float distance_max = 4.5;
};

inline ContactSet carbon_pi(const Topology& topology, const CarbonPiParams& params = CarbonPiParams{}) {
  return find_contacts(
      topology,
      [](const AtomRec &rec) {
        bool is_C_atom = rec.atom.getAtomicNum() == Element::C;
        bool is_whd = (rec.type & AtomType::WeakHbondDonor) == AtomType::WeakHbondDonor;
        return is_C_atom && is_whd;
      },
      [](const RingRec &rec) { return rec.aromatic; },
      {params.distance_max, 0.4},
      [&topology, &params](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist_sq) -> InteractionType {
        const auto &mol = topology.molecule();

        const auto atom_a_rec = topology.atom(rec_idx_a);
        const auto atom_b_rec = topology.ring(rec_idx_b);
        const auto *atom_a = mol.getAtomWithIdx(atom_a_rec.idx);
        const auto *atom_b = mol.getAtomWithIdx(atom_b_rec.atoms[0]);

        auto angle = compute_angle(atom_b_rec, topology.conformer().getAtomPos(atom_a->getIdx()));

        const auto angle_cutoff = 30.0;
        if (atom_a->getIsAromatic()) return InteractionType::None;
        if (!passes_angle_filter(angle, angle_cutoff)) return InteractionType::None;

        // FIX: need to make sure we remove metals here

        /*if (are_residueids_close(mol, *atom_a, *atom_b, 1)) return InteractionType::None;*/
        return InteractionType::Hydrophobic;
      }
  );
}

struct CationPiParams {
  float distance_max = 4.5;
};

inline ContactSet cation_pi(const Topology& topology, const CationPiParams& params = CationPiParams{}) {
  return find_contacts(
      topology,
      [](const AtomRec &rec) { return (rec.type & AtomType::PositiveCharge) == AtomType::PositiveCharge; },
      [](const RingRec &rec) { return rec.aromatic; },
      {params.distance_max, 0.4},
      [&topology, &params](std::uint32_t rec_idx_a, std::uint32_t rec_idx_b, float dist_sq) -> InteractionType {
        const auto &mol = topology.molecule();

        const auto atom_a_rec = topology.atom(rec_idx_a);
        const auto atom_b_rec = topology.ring(rec_idx_b);
        const auto *atom_a = mol.getAtomWithIdx(atom_a_rec.idx);
        const auto *atom_b = mol.getAtomWithIdx(atom_b_rec.atoms[0]);

        auto angle = compute_angle(atom_b_rec, topology.conformer().getAtomPos(atom_a->getIdx()));

        const auto angle_cutoff = 30.0;
        if (atom_a->getIsAromatic()) return InteractionType::None;
        if (!passes_angle_filter(angle, angle_cutoff)) return InteractionType::None;
        if (are_residueids_close(mol, *atom_a, *atom_b, 1)) return InteractionType::None;

        return InteractionType::Hydrophobic;
      }
  );
}

} // namespace lahuta

#endif // LAHUTA_CONTACTS_ARPEGGIO_HPP
