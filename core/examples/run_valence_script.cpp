#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

#include <GraphMol/MonomerInfo.h>
#include <GraphMol/RWMol.h>

#include "chemistry/valence_model.hpp"
#include "lahuta.hpp"
#include "serialization/json.hpp"

using namespace lahuta;

namespace {

enum class AtomGeometry : int {
  Spherical           = 0,
  Terminal            = 1,
  Linear              = 2,
  Trigonal            = 3,
  Tetrahedral         = 4,
  TrigonalBiPyramidal = 5,
  Octahedral          = 6,
  SquarePlanar        = 7,
  Unknown             = 8
};

// clang-format off
const char *geometry_label(AtomGeometry geometry) {
  switch (geometry) {
    case AtomGeometry::Spherical:          return "Spherical";
    case AtomGeometry::Terminal:           return "Terminal";
    case AtomGeometry::Linear:             return "Linear";
    case AtomGeometry::Trigonal:           return "Trigonal";
    case AtomGeometry::Tetrahedral:        return "Tetrahedral";
    case AtomGeometry::TrigonalBiPyramidal:return "Trigonal Bi-Pyramidal";
    case AtomGeometry::Octahedral:         return "Octahedral";
    case AtomGeometry::SquarePlanar:       return "Square Planar";
    case AtomGeometry::Unknown:            return "Unknown";
  }
  return "Unknown";
}
// clang-format on

AtomGeometry geometry_from_hybridization(const RDKit::Atom &atom) {
  const auto hyb    = atom.getHybridization();
  const auto degree = atom.getDegree();
  switch (hyb) {
    case RDKit::Atom::HybridizationType::S:
      if (degree == 0) return AtomGeometry::Spherical;
      if (degree == 1) return AtomGeometry::Terminal;
      return AtomGeometry::Unknown;
    case RDKit::Atom::HybridizationType::SP:
      return AtomGeometry::Linear;
    case RDKit::Atom::HybridizationType::SP2:
      return AtomGeometry::Trigonal;
    case RDKit::Atom::HybridizationType::SP3:
      return AtomGeometry::Tetrahedral;
    case RDKit::Atom::HybridizationType::SP3D:
      return AtomGeometry::TrigonalBiPyramidal;
    case RDKit::Atom::HybridizationType::SP3D2:
      return AtomGeometry::Octahedral;
    default:
      return AtomGeometry::Unknown;
  }
}

int count_bonded_hydrogens(const RDKit::ROMol &mol, const RDKit::Atom &atom) {
  int count = 0;
  for (const auto &bond : mol.atomBonds(&atom)) {
    const auto *other = bond->getOtherAtom(&atom);
    if (other && other->getAtomicNum() == 1) ++count;
  }
  return count;
}

struct Args {
  std::string input_path;
  std::optional<std::string> output_path;
  std::string format = "json";
};

void print_usage(const char *prog) {
  std::cout << "Usage: " << prog << " <structure.cif|.cif.gz|.pdb|.pdb.gz|.mmcif> [--o <output.json>]\n"
            << "Options:\n"
            << "  --o, -o         Output file path (default: stdout)\n"
            << "  --format        Output format (only 'json' is supported)\n";
}

enum class ParseResult { Ok, Help, Error };

ParseResult parse_args(int argc, char **argv, Args &out) {
  if (argc < 2) {
    print_usage(argv[0]);
    return ParseResult::Error;
  }

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--help" || arg == "-h") {
      print_usage(argv[0]);
      return ParseResult::Help;
    }
    if (arg == "--o" || arg == "-o" || arg == "--output") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << arg << "\n";
        return ParseResult::Error;
      }
      out.output_path = std::string(argv[++i]);
      continue;
    }
    if (arg == "--format") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for --format\n";
        return ParseResult::Error;
      }
      out.format = std::string(argv[++i]);
      continue;
    }
    if (!arg.empty() && arg[0] == '-') {
      std::cerr << "Unknown option: " << arg << "\n";
      return ParseResult::Error;
    }
    if (out.input_path.empty()) {
      out.input_path = std::move(arg);
    } else {
      std::cerr << "Unexpected extra argument: " << arg << "\n";
      return ParseResult::Error;
    }
  }

  if (out.input_path.empty()) {
    std::cerr << "Missing input structure file\n";
    return ParseResult::Error;
  }
  if (out.format != "json") {
    std::cerr << "Unsupported format: " << out.format << " (only 'json' supported)\n";
    return ParseResult::Error;
  }
  return ParseResult::Ok;
}

} // namespace

int main(int argc, char **argv) {
  Args args;
  const auto parse_result = parse_args(argc, argv, args);
  if (parse_result == ParseResult::Help) return 0;
  if (parse_result == ParseResult::Error) return 1;

  if (!std::filesystem::exists(args.input_path)) {
    std::cerr << "Input file not found: " << args.input_path << "\n";
    return 1;
  }

  Luni luni(args.input_path);
  if (!luni.build_topology()) {
    std::cerr << "Failed to build topology from " << args.input_path << "\n";
    return 1;
  }

  const auto &mol     = luni.get_molecule();
  const int model_num = 1;
  const int unit_id   = 0;

  std::vector<int> formal_charge_in(mol.getNumAtoms(), 0);
  for (const auto atom : mol.atoms()) {
    formal_charge_in[atom->getIdx()] = atom->getFormalCharge();
  }

  ValenceModel valence_model;
  valence_model.apply(mol);

  std::ofstream out_file;
  std::ostream *out = &std::cout;
  if (args.output_path) {
    out_file.open(*args.output_path, std::ios::binary | std::ios::trunc);
    if (!out_file) {
      std::cerr << "Failed to open output file: " << *args.output_path << "\n";
      return 1;
    }
    out = &out_file;
  }

  *out << "[";
  bool first = true;
  for (const auto atom : mol.atoms()) {
    const auto idx   = atom->getIdx();
    const auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());

    const std::string chain_id  = info ? info->getChainId() : "";
    const std::string res_name  = info ? info->getResidueName() : "";
    const int res_seq           = info ? info->getResidueNumber() : 0;
    const std::string ins_code  = info ? info->getInsertionCode() : "";
    const std::string atom_name = info ? info->getName() : atom->getSymbol();
    const std::string alt_id    = info ? info->getAltLoc() : "";
    int atom_id                 = info ? info->getSerialNumber() : static_cast<int>(idx);
    if (atom_id <= 0) atom_id = static_cast<int>(idx);

    const int formal_charge     = formal_charge_in[idx];
    const int vm_charge         = atom->getFormalCharge();
    const int implicit_h        = atom->getNumCompImplicitHs();
    const int explicit_h        = count_bonded_hydrogens(mol, *atom);
    const int total_h           = implicit_h + explicit_h;
    const AtomGeometry geometry = geometry_from_hybridization(*atom);

    JsonBuilder json(256);
    json.key("model_num")
        .value(model_num)
        .key("unit_id")
        .value(unit_id)
        .key("chain_id")
        .value(chain_id)
        .key("res_name")
        .value(res_name)
        .key("res_seq")
        .value(res_seq)
        .key("ins_code")
        .value(ins_code)
        .key("atom_name")
        .value(atom_name)
        .key("alt_id")
        .value(alt_id)
        .key("element")
        .value(atom->getSymbol())
        .key("atom_id")
        .value(atom_id)
        .key("source_index")
        .value(idx)
        .key("formal_charge")
        .value(formal_charge)
        .key("vm_charge")
        .value(vm_charge)
        .key("implicit_h")
        .value(implicit_h)
        .key("total_h")
        .value(total_h)
        .key("ideal_geometry")
        .value(geometry_label(geometry))
        .key("ideal_geometry_code")
        .value(static_cast<int>(geometry));

    if (!first) *out << ",";
    *out << json.str();
    first = false;
  }
  *out << "]\n";

  return 0;
}
