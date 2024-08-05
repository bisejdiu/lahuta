#include <array>
#include <string_view>

// clang-format off

namespace lahuta {
namespace rules {

struct Rule {
    std::string_view s2;
    std::string_view s3;
    int result;
};

constexpr std::array<Rule, 3> def = {{
    {"C", "O", 2},
}};

constexpr std::array<Rule, 1> default_amino_acid = {{
    {"C", "O", 2},
}};

constexpr std::array<Rule, 1> default_base = {{
    {"OP1", "P", 2},
}};

constexpr std::array<Rule, 3> his_rules = {{
    {"CD2", "CG", 2},
    {"CE1", "ND1", 2},
    {"C", "O", 2}
}};

constexpr std::array<Rule, 2> arg_rule = {{
    {"CZ", "NH2", 2},
    {"C", "O", 2}
}};

constexpr std::array<Rule, 4> phe_rule = {{
    {"CE1", "CZ", 2},
    {"CD2", "CE2", 2},
    {"CD1", "CG", 2},
    {"C", "O", 2}
}};

constexpr std::array<Rule, 5> trp_rule = {{
    {"CD1", "CG", 2},
    {"CD2", "CE2", 2},
    {"CE3", "CZ3", 2},
    {"CH2", "CZ2", 2},
    {"C", "O", 2}
}};

constexpr std::array<Rule, 2> asn_rule = {{
    {"CG", "OD1", 2},
    {"C", "O", 2}
}};

constexpr std::array<Rule, 2> gln_rule = {{
    {"CD", "OE1", 2},
    {"C", "O", 2}
}};

constexpr std::array<Rule, 4> tyr_rule = {{
    {"CD1", "CG", 2},
    {"CD2", "CE2", 2},
    {"CE1", "CZ", 2},
    {"C", "O", 2}
}};

constexpr std::array<Rule, 2> asp_rule = {{
    {"CG", "OD1", 2},
    {"C", "O", 2}
}};

constexpr std::array<Rule, 2> glu_rule = {{
    {"CD", "OE1", 2},
    {"C", "O", 2}
}};

constexpr std::array<Rule, 5> g_rule = {{
    {"C8", "N7", 2},
    {"C4", "C5", 2},
    {"C2", "N3", 2},
    {"C6", "O6", 2},
    {"OP1", "P", 2}
}};


constexpr std::array<Rule, 4> c_rule = {{
    {"C4", "N3", 2},
    {"C5", "C6", 2},
    {"C2", "O2", 2},
    {"OP1", "P", 2}
}};


constexpr std::array<Rule, 5> a_rule = {{
    {"C2", "N3", 2},
    {"C6", "N1", 2},
    {"C4", "C5", 2},
    {"C8", "N7", 2},
    {"OP1", "P", 2}
}};

constexpr std::array<Rule, 4> u_rule = {{
    {"C5", "C6", 2},
    {"C2", "O2", 2},
    {"C4", "O4", 2},
    {"OP1", "P", 2}
}};

constexpr std::array<Rule, 5> dg_rule = {{
    {"C8", "N7", 2},
    {"C4", "C5", 2},
    {"C2", "N3", 2},
    {"C6", "O6", 2},
    {"OP1", "P", 2}
}};

constexpr std::array<Rule, 4> dc_rule = {{
    {"C4", "N3", 2},
    {"C5", "C6", 2},
    {"C2", "O2", 2},
    {"OP1", "P", 2}
}};

constexpr std::array<Rule, 5> da_rule = {{
    {"C2", "N3", 2},
    {"C6", "N1", 2},
    {"C4", "C5", 2},
    {"C8", "N7", 2},
    {"OP1", "P", 2}
}};

constexpr std::array<Rule, 4> dt_rule = {{
    {"C5", "C6", 2},
    {"C2", "O2", 2},
    {"C4", "O4", 2},
    {"OP1", "P", 2}
}};


} // namespace rules
} // namespace lahuta
