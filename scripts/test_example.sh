#!/usr/bin/env bash
set -euo pipefail

BIN="./bin/cdbg"
if [[ ! -x "$BIN" ]]; then
    make all
fi

echo "=== 1. Building Colored Compacted De Bruijn Graph ==="
$BIN build -k 31 -o pangenome.cdbg data/sample1.fa data/sample2.fa data/sample3.fa

echo -e "\n=== 2. Graph and Filter Info ==="
$BIN info -i pangenome.cdbg

echo -e "\n=== 3. Querying Sequences ==="
$BIN query -i pangenome.cdbg -q data/query.fa

echo -e "\n=== 4. Exporting GFA v1 ==="
$BIN export-gfa -i pangenome.cdbg -o pangenome.gfa
head -n 20 pangenome.gfa

rm -f pangenome.cdbg pangenome.gfa
echo -e "\nExample run completed successfully!"
