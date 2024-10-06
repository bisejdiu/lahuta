#include "lahuta.hpp"
#include "parser.hpp"

namespace lahuta {

namespace {
// Check if a substring at 'pos' matches an operator
bool is_operator(const std::string &word, size_t pos, const std::string &op) {
  return word.compare(pos, op.length(), op) == 0;
}

// Split a word into tokens (operators and identifiers)
std::vector<std::string> split_word(const std::string &word) {
  std::vector<std::string> tokens;
  size_t pos = 0;

  while (pos < word.length()) {
    if (word.compare(pos, 3, "and") == 0) {
      tokens.push_back("and");
      pos += 3;
    } else if (word.compare(pos, 2, "or") == 0) {
      tokens.push_back("or");
      pos += 2;
    } else if (word.compare(pos, 3, "not") == 0) {
      tokens.push_back("not");
      pos += 3;
    } else if (word.compare(pos, 5, "resid") == 0) {
      tokens.push_back("resid");
      pos += 5;
    } else if (word.compare(pos, 7, "resname") == 0) {
      tokens.push_back("resname");
      pos += 7;
    } else if (std::isalpha(word[pos])) {
      // Extract letters
      size_t start = pos;
      while (pos < word.length() && std::isalpha(word[pos])) {
        ++pos;
      }
      tokens.push_back(word.substr(start, pos - start));
    } else if (std::isdigit(word[pos])) {
      // Extract digits
      size_t start = pos;
      while (pos < word.length() && std::isdigit(word[pos])) {
        ++pos;
      }
      tokens.push_back(word.substr(start, pos - start));
    } else if (word[pos] == '-') {
      tokens.push_back("-");
      ++pos;
    } else {
      // Unknown character, skip or handle error
      ++pos;
    }
  }

  return tokens;
}
} // namespace

const std::vector<std::string> Luni::symbols() const {
  return atom_attrs<std::string>(
      [](const RDKit::Atom *atom) { return atom->getSymbol(); });
}

const std::vector<std::string> Luni::names() const {
  return atom_attrs<std::string>(
      [](const RDKit::Atom *atom) -> const std::string & {
        auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(
            atom->getMonomerInfo());
        return info->getName();
      });
}

const std::vector<int> Luni::indices() const {
  return atom_attrs<int>(
      [](const RDKit::Atom *atom) { return atom->getIdx(); });
}

const std::vector<int> Luni::atomic_numbers() const {
  return atom_attrs<int>(
      [](const RDKit::Atom *atom) { return atom->getAtomicNum(); });
}

const std::vector<std::string> Luni::elements() const {
  const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
  return atom_attrs<std::string>([&tbl](const RDKit::Atom *atom) {
    return tbl->getElementSymbol(atom->getAtomicNum());
  });
}

const std::vector<std::string> Luni::resnames() const {
  return atom_attrs<std::string>(
      [](const RDKit::Atom *atom) -> const std::string & {
        auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(
            atom->getMonomerInfo());
        return info->getResidueName();
      });
}

const std::vector<int> Luni::resids() const {
  return atom_attrs<int>([](const RDKit::Atom *atom) -> int {
    auto *info =
        static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
    return info->getResidueNumber();
  });
}

const std::vector<int> Luni::resindices() const {
  return atom_attrs<int>([](const RDKit::Atom *atom) -> int {
    auto *info =
        static_cast<const RDKit::AtomPDBResidueInfo *>(atom->getMonomerInfo());
    return info->getSegmentNumber();
  });
}

const std::vector<std::string> Luni::chainlabels() const {
  return atom_attrs<std::string>(
      [](const RDKit::Atom *atom) -> const std::string & {
        auto *info = static_cast<const RDKit::AtomPDBResidueInfo *>(
            atom->getMonomerInfo());
        return info->getChainId();
      });
}

template <typename T>
std::vector<T>
Luni::atom_attrs(std::function<T(const RDKit::Atom *)> func) const {
  std::vector<T> attrs;
  attrs.reserve(mol->getNumAtoms());
  for (const auto atom : mol->atoms()) {
    attrs.push_back(func(atom));
  }
  return attrs;
}

template <typename T>
std::vector<std::reference_wrapper<const T>>
Luni::atom_attrs_ref(std::function<const T &(const RDKit::Atom *)> func) const {
  std::vector<std::reference_wrapper<const T>> attributes;
  attributes.reserve(mol->getNumAtoms());
  for (const auto atom : mol->atoms()) {
    attributes.push_back(std::cref(func(atom)));
  }
  return attributes;
}

std::vector<RDKit::MatchVectType>
Luni::match_smarts_string(std::string sm, std::string atype,
                          bool log_values) const {

  if (!mol->getRingInfo()->isInitialized()) {
    RDKit::MolOps::symmetrizeSSSR(*mol);
  }
  std::vector<RDKit::MatchVectType> match_list;
  auto sm_mol = RDKit::SmartsToMol(sm);
  mol->updatePropertyCache(false);
  RDKit::SubstructMatch(*mol, *sm_mol, match_list);
  return match_list;
};

NSResults Luni::find_neighbors_opt(double cutoff) {

  if (cutoff == _cutoff) {
    return bonded_nps;
  } else if (cutoff < _cutoff) {
    return bonded_nps.filter(cutoff);
  }

  grid.update_cutoff(cutoff);
  auto ns = grid.self_search();
  return ns;
}

Luni Luni::filter_luni(const std::vector<int> &atom_indices) const {

  auto new_mol = filter_with_bonds(*mol, atom_indices);
  new_mol.updatePropertyCache(false);

  Luni new_luni;
  new_luni.mol = std::make_shared<RDKit::RWMol>(new_mol);
  new_luni.bonded_nps = bonded_nps.filter(atom_indices);
  new_luni.topology.assign_atom_types(new_mol);

  return new_luni;
}

std::vector<std::string> Luni::tokenize_simple(const std::string &str) {
  std::vector<std::string> tokens;
  std::string token;

  std::istringstream iss(str);
  while (iss >> token) {
    tokens.push_back(token);
  }
  return tokens;
}

std::vector<std::string> Luni::tokenize(const std::string &str) {
  std::vector<std::string> tokens;
  size_t i = 0;

  while (i < str.length()) {
    char c = str[i];

    // Skip whitespace
    if (std::isspace(c)) {
      ++i;
      continue;
    }

    // Parentheses
    if (c == '(' || c == ')') {
      tokens.push_back(std::string(1, c));
      ++i;
      continue;
    }

    // Dash or negative sign
    if (c == '-') {
      bool is_negative = false;
      if (tokens.empty() || tokens.back() == "(" || tokens.back() == "and" ||
          tokens.back() == "or" || tokens.back() == "not" ||
          tokens.back() == "resid" || tokens.back() == "resname") {
        is_negative = true;
      }

      if (is_negative && i + 1 < str.length() && std::isdigit(str[i + 1])) {
        // Negative number
        size_t start = i;
        ++i; // Move past '-'
        while (i < str.length() && std::isdigit(str[i])) {
          ++i;
        }
        tokens.push_back(str.substr(start, i - start));
      } else {
        // Dash operator
        tokens.push_back("-");
        ++i;
      }
      continue;
    }

    // Identifiers and operators
    if (std::isalpha(c)) {
      size_t start = i;
      while (i < str.length() &&
             (std::isalpha(str[i]) || std::isdigit(str[i]) || str[i] == '-')) {
        ++i;
      }
      std::string word = str.substr(start, i - start);

      // Split the word into tokens
      std::vector<std::string> word_tokens = split_word(word);

      // Append the tokens
      tokens.insert(tokens.end(), word_tokens.begin(), word_tokens.end());
      continue;
    }

    // Numbers
    if (std::isdigit(c)) {
      size_t start = i;
      while (i < str.length() && std::isdigit(str[i])) {
        ++i;
      }
      tokens.push_back(str.substr(start, i - start));
      continue;
    }

    // Handle any other character sequences
    size_t start = i;
    while (i < str.length() && !std::isspace(str[i]) && str[i] != '(' &&
           str[i] != ')' && str[i] != '-') {
      ++i;
    }
    tokens.push_back(str.substr(start, i - start));
  }

  return tokens;
}

bool Luni::parse_expression(const std::string &selection) {
  try {
    std::vector<std::string> tokens = Luni::tokenize(selection);
    lahuta::Parser parser(tokens);

    // Parse the expression
    lahuta::NodePtr root = parser.parse_expression();

    lahuta::FilterVisitor visitor(*this);
    root->accept(visitor);

    filtered_indices = visitor.get_result();
    return true;
  } catch (const std::exception &e) {
    // Handle parsing errors
    filtered_indices.clear();
    return false;
  }
}

Luni Luni::filter() const {
  if (filtered_indices.empty()) {
    std::cerr << "Selection not parsed or empty" << std::endl;
  }
  return filter_luni(filtered_indices);
}

std::vector<int> Luni::parse_and_filter(const std::string &selection) const {
  std::vector<std::string> tokens = Luni::tokenize(selection);
  lahuta::Parser parser(tokens);

  // Parse the expression
  lahuta::NodePtr root = parser.parse_expression();

  lahuta::FilterVisitor visitor(*this);
  root->accept(visitor);

  const std::vector<int> &filtered_indices = visitor.get_result();
  return filtered_indices;
}

NSResults Luni::remove_adjascent_residueid_pairs(NSResults &results,
                                                 int res_diff) {
  Pairs filtered;
  Distances dists;
  for (size_t i = 0; i < results.get_pairs().size(); ++i) {
    auto *fatom = get_molecule().getAtomWithIdx(results.get_pairs()[i].first);
    auto *satom = get_molecule().getAtomWithIdx(results.get_pairs()[i].second);

    auto *finfo =
        static_cast<const RDKit::AtomPDBResidueInfo *>(fatom->getMonomerInfo());
    auto *sinfo =
        static_cast<const RDKit::AtomPDBResidueInfo *>(satom->getMonomerInfo());

    if (fatom->getAtomicNum() == 1 || satom->getAtomicNum() == 1)
      continue; // skip H atoms (hydrogens)

    auto f_resid = finfo->getResidueNumber();
    auto s_resid = sinfo->getResidueNumber();

    if (std::abs(f_resid - s_resid) > res_diff) {
      filtered.push_back(results.get_pairs()[i]);
      dists.push_back(results.get_distances()[i]);
    }
  }
  return NSResults(filtered, dists);
}

std::vector<int> Luni::factorize(const std::vector<std::string>& labels) {
    std::vector<int> ids(labels.size());

    // Hash map from labels to ids
    std::unordered_map<std::string_view, int> label_to_id;
    label_to_id.reserve(labels.size());

    int current_id = 0;
    for (size_t i = 0; i < labels.size(); ++i) {
        std::string_view label = labels[i];
        auto it = label_to_id.find(label);
        if (it == label_to_id.end()) {
            label_to_id[label] = current_id;
            ids[i] = current_id;
            ++current_id;
        } else {
            ids[i] = it->second;
        }
    }

    return ids;
}

int Luni::count_unique(const std::vector<int>& vec) {
    std::unordered_set<int> unique_elements(vec.begin(), vec.end());
    return unique_elements.size();
}

int Luni::count_unique(const std::vector<std::string>& vec) {
    std::unordered_set<std::string_view> unique_elements;
    unique_elements.reserve(vec.size()); 

    for (const auto& str : vec) {
        unique_elements.insert(std::string_view(str));
    }

    return unique_elements.size();
}

std::vector<std::string> Luni::find_elements(const std::vector<int>& atomic_numbers) {
    const RDKit::PeriodicTable *tbl = RDKit::PeriodicTable::getTable();
    std::vector<std::string> elements;
    elements.reserve(atomic_numbers.size());

    for (int atomic_number : atomic_numbers) {
        elements.push_back(tbl->getElementSymbol(atomic_number));
    }

    return elements;
}

} // namespace lahuta
