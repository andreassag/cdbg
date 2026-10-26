#!/usr/bin/env bash
set -euo pipefail

BIN="./bin/cdbg"
if [[ ! -x "$BIN" ]]; then
    make all
fi

echo "=== Benchmark: cDBG Graph Construction & Succinct Filter Overhead ==="

echo "1. Building Graph & Peeling Xor Filter (k=31, 3 genomes):"
/usr/bin/time -v $BIN build -k 31 -o /tmp/bench.cdbg data/sample1.fa data/sample2.fa data/sample3.fa 2>&1 | grep -E "Elapsed|Maximum resident"

echo -e "\n2. Querying 510 k-mers against graph with succinct filter:"
/usr/bin/time -v $BIN query -i /tmp/bench.cdbg -q data/query.fa > /dev/null 2>&1

echo -e "\n3. Filter Size Efficiency:"
$BIN info -i /tmp/bench.cdbg | grep -E "Unitigs|Filter"

rm -f /tmp/bench.cdbg
echo "Benchmark completed."
