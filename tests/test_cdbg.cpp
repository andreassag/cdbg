#include "cdbg/dna.hpp"
#include "cdbg/bit_vector.hpp"
#include "cdbg/succinct_filter.hpp"
#include "cdbg/unitig_graph.hpp"
#include "cdbg/index_io.hpp"
#include <iostream>
#include <cassert>
#include <random>

void test_bit_vector_rank_select() {
    cdbg::BitVector bv(2000); // multiple superblocks
    std::vector<size_t> set_indices = {5, 12, 63, 64, 511, 512, 1000, 1500, 1999};

    for (size_t idx : set_indices) {
        bv.set(idx);
    }
    bv.build_index();

    assert(bv.total_set_bits() == set_indices.size());

    // Test rank1
    assert(bv.rank1(0) == 0);
    assert(bv.rank1(6) == 1);
    assert(bv.rank1(13) == 2);
    assert(bv.rank1(64) == 3);
    assert(bv.rank1(65) == 4);
    assert(bv.rank1(512) == 5);
    assert(bv.rank1(513) == 6);
    assert(bv.rank1(2000) == set_indices.size());

    // Test select1 (1-based k)
    for (size_t i = 0; i < set_indices.size(); ++i) {
        size_t found_pos = bv.select1(i + 1);
        assert(found_pos == set_indices[i]);
    }

    std::cout << "[PASS] test_bit_vector_rank_select\n";
}

void test_blocked_bloom_filter() {
    cdbg::BlockedBloomFilter filter(1000, 0.01);
    std::mt19937_64 rng(42);
    std::vector<uint64_t> keys(1000);
    for (auto& k : keys) {
        k = rng();
        filter.insert(k);
    }

    // 100% sensitivity (no false negatives)
    for (uint64_t k : keys) {
        assert(filter.contains(k));
    }

    // Measure false positive rate
    size_t fp_count = 0;
    size_t test_count = 10000;
    for (size_t i = 0; i < test_count; ++i) {
        uint64_t query = rng() | 0x8000000000000000ULL;
        if (filter.contains(query)) {
            fp_count++;
        }
    }
    double fpr = static_cast<double>(fp_count) / static_cast<double>(test_count);
    assert(fpr < 0.03); // Expected ~1%, allow margin

    std::cout << "[PASS] test_blocked_bloom_filter (measured FPR: " << (fpr * 100.0) << "%)\n";
}

void test_xor_filter_8() {
    std::mt19937_64 rng(1234);
    std::vector<uint64_t> keys(1000);
    for (auto& k : keys) k = rng();

    cdbg::XorFilter8 filter(keys);

    // Guaranteed ZERO false negatives
    for (uint64_t k : keys) {
        assert(filter.contains(k));
    }

    // Measure false positive rate (theoretically ~1/256 = 0.39%)
    size_t fp_count = 0;
    size_t test_count = 10000;
    for (size_t i = 0; i < test_count; ++i) {
        uint64_t query = rng() ^ 0xFFFFFFFFFFFFFFFFULL;
        if (filter.contains(query)) {
            fp_count++;
        }
    }
    double fpr = static_cast<double>(fp_count) / static_cast<double>(test_count);
    assert(fpr < 0.015);

    std::cout << "[PASS] test_xor_filter_8 (measured FPR: " << (fpr * 100.0) << "%)\n";
}

void test_unitig_compaction() {
    // Linear non-branching sequence of 100 bp without internal repeats
    std::mt19937 rng(42);
    const char bases[] = {'A', 'C', 'G', 'T'};
    std::string seq(100, 'A');
    for (size_t i = 0; i < seq.size(); ++i) {
        seq[i] = bases[rng() % 4];
    }
    uint32_t k = 15;

    std::vector<std::pair<std::string, std::string>> samples = {{"sample1", seq}};
    auto graph = cdbg::UnitigGraph::build(samples, k);

    assert(graph.unitigs.size() == 1);
    assert(graph.unitigs[0].seq.size() == seq.size());

    std::cout << "[PASS] test_unitig_compaction (100 bp compacted to 1 unitig)\n";
}

void test_color_matrix_multi_sample() {
    uint32_t k = 11;
    // Sample A: prefix + shared + suffixA
    std::string shared = "ACGTACGTACGTACGTACGT";
    std::string seqA = "AAAA" + shared + "TTTT";
    std::string seqB = "CCCC" + shared + "GGGG";

    std::vector<std::pair<std::string, std::string>> samples = {
        {"SampleA", seqA},
        {"SampleB", seqB}
    };

    auto graph = cdbg::UnitigGraph::build(samples, k);
    assert(graph.unitigs.size() >= 3); // private A, shared, private B

    // Query shared k-mer
    std::string shared_kmer = shared.substr(0, k);
    auto hit_shared = graph.query_kmer(shared_kmer);
    assert(hit_shared.found);
    assert(hit_shared.present_samples.size() == 2);

    // Query private A k-mer
    std::string priv_a = seqA.substr(0, k);
    auto hit_a = graph.query_kmer(priv_a);
    assert(hit_a.found);
    assert(hit_a.present_samples.size() == 1);
    assert(hit_a.present_samples[0] == "SampleA");

    std::cout << "[PASS] test_color_matrix_multi_sample\n";
}

void test_serialization() {
    std::string seq = "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATC";
    uint32_t k = 11;
    auto graph = cdbg::UnitigGraph::build({{"Sample1", seq}}, k);

    std::string tmp_file = "/tmp/test_cdbg_io.cdbg";
    cdbg::save_cdbg(graph, tmp_file);

    auto loaded = cdbg::load_cdbg(tmp_file);
    assert(loaded.k == graph.k);
    assert(loaded.unitigs.size() == graph.unitigs.size());
    assert(loaded.unitigs[0].seq == graph.unitigs[0].seq);
    assert(loaded.colors.num_samples() == graph.colors.num_samples());

    // Query on loaded graph
    auto hit = loaded.query_kmer(seq.substr(0, k));
    assert(hit.found);

    std::cout << "[PASS] test_serialization\n";
}

int main() {
    std::cout << "Running cdbg unit tests...\n";
    test_bit_vector_rank_select();
    test_blocked_bloom_filter();
    test_xor_filter_8();
    test_unitig_compaction();
    test_color_matrix_multi_sample();
    test_serialization();
    std::cout << "All cdbg unit tests passed successfully!\n";
    return 0;
}
