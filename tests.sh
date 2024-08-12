#!/usr/bin/env bash

set -e

usage() {
    echo "Usage: $0 <executable> <data_dir>"
    exit 1
}

if [ $# -ne 2 ]; then
    usage
fi

executable=$1
data_dir=$2
output_file="OUT"
expected_results="$data_dir/results_arom.txt"

trap "rm -f $output_file" EXIT

if [ ! -f "$expected_results" ]; then
    echo "Error: Expected results file '$expected_results' not found."
    exit 1
fi

# Run the tests
for input_file in "$data_dir"/*.cif; do
    echo "Running test on $input_file..."
    if ! ./"$executable" "$input_file" | grep "^1." >> "$output_file"; then
        echo "Error: Executable failed on file '$input_file'."
        exit 1
    fi
done

if diff "$output_file" "$expected_results" > /dev/null; then
    echo "All tests passed."
    rm "$output_file"
else
    echo "Test failed. Differences found:"
    diff "$output_file" "$expected_results"
    exit 1
fi

