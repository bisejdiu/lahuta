#!/usr/bin/env python

import re


def parse_enum(file_path):
    """Parse the resTokenType enum from the given file."""
    enum_values = []
    in_enum = False

    with open(file_path, "r") as file:
        for line in file:
            stripped_line = line.strip()
            if not in_enum:
                if stripped_line.startswith("enum class resTokenType"):
                    in_enum = True
            else:
                if stripped_line == "};":
                    break
                # Remove comments and split by comma
                value = re.split(
                    r"\s*//.*", stripped_line)[0].strip().rstrip(",")
                if value and value != "UNKNOWN":
                    enum_values.append(value)

    return enum_values


def generate_gperf_input(enum_values):
    """Generate the gperf input from the list of enum values."""
    gperf_input = """%compare-strncmp
%define class-name resNameKeyword
%define initializer-suffix ,static_cast<resTokenType>(0)
%define lookup-function-name look_up
%define slot-name string_offset
%enum
%language=C++
%pic
%readonly-tables
%struct-type

%{
#include <stddef.h>
#include <string.h>
#include "token.h"

namespace lahuta {
namespace {
%}
struct keywordEntry {
  int string_offset;
  resTokenType type;
};

%%
"""
    for value in enum_values:
        key = value + ","
        gperf_input += f"{key:<8}resTokenType::{value}\n"

    gperf_input += """%%
}

#ifndef LOOKUP_IDENTIFIER_DEFINED
#define LOOKUP_IDENTIFIER_DEFINED

inline resTokenType res_name_table(const char* res_name, std::size_t size) noexcept {
  const keywordEntry *entry = resNameKeyword::look_up(res_name, size);
  if (entry) {
    return entry->type;
  } else {
    return resTokenType::UNKNOWN;
  }
}
}
#endif
"""

    return gperf_input


def main():
    input_file = "token.h"
    output_file = input_file.replace(".h", ".gperf")

    enum_values = parse_enum(input_file)
    gperf_input = generate_gperf_input(enum_values)

    with open(output_file, "w") as file:
        file.write(gperf_input)

    print(f"gperf input file generated successfully as {output_file}")


if __name__ == "__main__":
    main()
