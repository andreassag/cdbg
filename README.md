# cdbg: Colored Compacted De Bruijn Graph with Succinct Filters (`bifrost`-lite)

[![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](https://en.cppreference.com/w/cpp/20)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Zero Dependencies](https://img.shields.io/badge/dependencies-none-brightgreen.svg)]()

`cdbg` is a high-performance C++20 colored compacted de Bruijn graph engine designed for pangenomics and sequence membership indexing. It extracts maximal non-branching unitigs across multiple input genomes, maps unitigs to deduplicated bit-sliced sample color equivalence classes, and embeds dual succinct membership filters: a 512-bit cache-line **Blocked Bloom Filter** and a linear-peeled **8-bit Xor Filter** (Graf & Lemire 2020) with guaranteed zero false negatives.

---

## Key Features

1. **Maximal Unitig Compaction**:
   - Assembles contiguous non-branching de Bruijn paths where internal nodes satisfy $\text{in-degree} = 1$, $\text{out-degree} = 1$, and identical sample color membership.
2. **Succinct BitVector with 512-Bit Superblock Rank**:
   - Precomputed cumulative popcount ranks over 512-bit (64-byte) superblocks providing $O(1)$ `rank1(pos)` and binary-searched `select1(k)`.
3. **Dual Succinct Membership Filters**:
   - **Cache-Line Blocked Bloom Filter**: Maps 64-bit $k$-mer hashes into 64-byte L1-cache lines, setting 4 intra-block bits for instant cache-hot rejection.
   - **Xor Filter 8-bit (Graf & Lemire)**: Constructed via linear-time 3-hypergraph peeling in $\lfloor 1.23 \cdot N \rfloor + 32$ slots. Achieves **zero false negatives** and ~0.39% theoretical false positive rate in under 10 bits per $k$-mer.
4. **Bit-Sliced Color Matrix**:
   - Represents sample presence across $C$ input genomes as bit slices, deduplicating recurring presence patterns into compact equivalence class IDs.
5. **Serialization & GFA v1 Export**:
   - Compact binary `.cdbg` format and standard Graphical Fragment Assembly (**GFA v1**) export with `colors:Z:` tags.

---

## Architectural Overview

```
Input Genomes (Sample 1..C)
           |
   [ Canonical k-mer Extraction ]
           |
   [ Non-Branching Path Compaction ]
           |
           +-----------------------------+
           |                             |
   [ Unitig Graph ]              [ Color Matrix ]
   (Sequences & Edges)           (Bit-Sliced Equivalence Classes)
           |                             |
           +--------------+--------------+
                          |
             [ Succinct Filter Embedding ]
             - Blocked Bloom Filter (64-byte cache lines)
             - Xor Filter 8 (Peeling construction)
                          |
             +------------+------------+
             |                         |
      Binary Index (.cdbg)        GFA v1 Export (.gfa)
```

---

## Build & Test

### Building via Makefile
```bash
make all      # Builds bin/cdbg
make test     # Runs comprehensive unit tests
```

### Building via CMake
```bash
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
ctest --output-on-failure
```

---

## CLI Usage

### 1. Build Index from Multiple Genomes
```bash
./bin/cdbg build -k 31 -o pangenome.cdbg sample1.fa sample2.fa sample3.fa
```

### 2. Inspect Topology & Filter Memory
```bash
./bin/cdbg info -i pangenome.cdbg
```
Output:
```text
=== Colored Compacted De Bruijn Graph Info ===
k-mer size:            31
Number of Unitigs:     7
Total Unitig Bases:    14000 bp
Number of Samples:     3
Color Equivalence Cls: 3
Bloom Filter Size:     16832 bytes
Xor Filter Size:       17280 bytes
```

### 3. Query Sequences for Membership and Color Presence
```bash
./bin/cdbg query -i pangenome.cdbg -q query.fa
```

### 4. Export to GFA v1
```bash
./bin/cdbg export-gfa -i pangenome.cdbg -o graph.gfa
```

---

## Verification & Benchmarks

Run automated verification:
```bash
./scripts/test_example.sh
./scripts/benchmark_compare.sh
```

### Performance on 3 Genomes (24 kb Pangenome):
- **Graph Construction + Xor Filter Peeling**: ~40 ms elapsed time.
- **Succinct Filter False Positive Rate**: 0.41% (Xor filter), 1.37% (Bloom filter).
- **Peak RSS**: ~16 MB.

---

## License

MIT License. Copyright (c) 2026 Andreas Sagen.
