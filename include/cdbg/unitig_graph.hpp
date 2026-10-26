#pragma once

#include "dna.hpp"
#include "fasta.hpp"
#include "color_matrix.hpp"
#include "succinct_filter.hpp"
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <algorithm>

namespace cdbg {

struct Unitig {
    uint32_t id = 0;
    std::string seq;
    uint32_t color_class = 0;
    std::vector<uint32_t> out_edges;
};

struct KmerHit {
    bool found = false;
    uint32_t unitig_id = 0;
    uint32_t offset = 0;
    uint32_t color_class = 0;
    std::vector<std::string> present_samples;
};

class UnitigGraph {
public:
    uint32_t k = 31;
    std::vector<Unitig> unitigs;
    ColorMatrix colors;
    XorFilter8 xor_filter;
    BlockedBloomFilter bloom_filter;

    UnitigGraph() = default;

    static UnitigGraph build(const std::vector<std::pair<std::string, std::string>>& samples, uint32_t k_val) {
        UnitigGraph graph;
        graph.k = k_val;

        std::vector<std::string> sample_names;
        sample_names.reserve(samples.size());
        for (const auto& s : samples) sample_names.push_back(s.first);
        graph.colors = ColorMatrix(sample_names);

        struct KmerNode {
            std::vector<bool> sample_presence;
            std::unordered_set<std::string> next;
            std::unordered_set<std::string> prev;
        };

        std::unordered_map<std::string, KmerNode> kmers;
        size_t num_samples = samples.size();

        for (size_t s_idx = 0; s_idx < num_samples; ++s_idx) {
            const std::string& seq = samples[s_idx].second;
            if (seq.size() < k_val) continue;

            std::string prev_kmer;
            for (size_t i = 0; i + k_val <= seq.size(); ++i) {
                std::string kmer = seq.substr(i, k_val);
                bool valid = true;
                for (char c : kmer) {
                    if (base_to_2bit(c) >= 4) { valid = false; break; }
                }
                if (!valid) {
                    prev_kmer.clear();
                    continue;
                }

                auto& node = kmers[kmer];
                if (node.sample_presence.empty()) {
                    node.sample_presence.resize(num_samples, false);
                }
                node.sample_presence[s_idx] = true;

                if (!prev_kmer.empty()) {
                    kmers[prev_kmer].next.insert(kmer);
                    node.prev.insert(prev_kmer);
                }
                prev_kmer = kmer;
            }
        }

        // Compact into maximal non-branching unitigs
        std::unordered_set<std::string> visited;
        std::vector<uint64_t> all_kmer_hashes;

        for (const auto& [kmer, node] : kmers) {
            if (visited.count(kmer)) continue;

            // A unitig starts if in-degree != 1, or predecessor has out-degree != 1,
            // or predecessor has different colors
            bool is_start = false;
            if (node.prev.size() != 1) {
                is_start = true;
            } else {
                const std::string& p_kmer = *node.prev.begin();
                const auto& p_node = kmers[p_kmer];
                if (p_node.next.size() != 1 || p_node.sample_presence != node.sample_presence) {
                    is_start = true;
                }
            }

            if (!is_start) continue;

            // Trace forward non-branching path
            std::string unitig_seq = kmer;
            visited.insert(kmer);
            std::string curr = kmer;

            while (true) {
                const auto& cur_node = kmers[curr];
                if (cur_node.next.size() != 1) break;

                const std::string& nxt = *cur_node.next.begin();
                if (visited.count(nxt)) break;

                const auto& next_node = kmers[nxt];
                if (next_node.prev.size() != 1 || next_node.sample_presence != node.sample_presence) break;

                unitig_seq.push_back(nxt.back());
                visited.insert(nxt);
                curr = nxt;
            }

            uint32_t cid = graph.colors.get_or_create_class(node.sample_presence);
            Unitig u;
            u.id = static_cast<uint32_t>(graph.unitigs.size());
            u.seq = unitig_seq;
            u.color_class = cid;
            graph.unitigs.push_back(std::move(u));
        }

        // Trace any remaining cycles (e.g. isolated circular components)
        for (const auto& [kmer, node] : kmers) {
            if (visited.count(kmer)) continue;

            std::string unitig_seq = kmer;
            visited.insert(kmer);
            std::string curr = kmer;

            while (true) {
                const auto& cur_node = kmers[curr];
                if (cur_node.next.size() != 1) break;

                const std::string& nxt = *cur_node.next.begin();
                if (visited.count(nxt)) break;

                const auto& next_node = kmers[nxt];
                if (next_node.prev.size() != 1 || next_node.sample_presence != node.sample_presence) break;

                unitig_seq.push_back(nxt.back());
                visited.insert(nxt);
                curr = nxt;
            }

            uint32_t cid = graph.colors.get_or_create_class(node.sample_presence);
            Unitig u;
            u.id = static_cast<uint32_t>(graph.unitigs.size());
            u.seq = unitig_seq;
            u.color_class = cid;
            graph.unitigs.push_back(std::move(u));
        }

        // Construct out_edges by matching (k-1)-mer overlaps
        if (graph.k > 1) {
            std::unordered_map<std::string, std::vector<uint32_t>> prefix_map;
            prefix_map.reserve(graph.unitigs.size());
            for (const auto& u : graph.unitigs) {
                if (u.seq.size() >= graph.k) {
                    std::string pfx = u.seq.substr(0, graph.k - 1);
                    prefix_map[pfx].push_back(u.id);
                }
            }

            for (auto& u : graph.unitigs) {
                if (u.seq.size() >= graph.k) {
                    std::string sfx = u.seq.substr(u.seq.size() - (graph.k - 1));
                    auto it = prefix_map.find(sfx);
                    if (it != prefix_map.end()) {
                        for (uint32_t tgt : it->second) {
                            u.out_edges.push_back(tgt);
                        }
                    }
                }
            }
        }

        // Populate succinct filters with canonical k-mer hashes
        for (const auto& u : graph.unitigs) {
            for (size_t i = 0; i + k_val <= u.seq.size(); ++i) {
                std::string km = canonical_kmer(u.seq.substr(i, k_val));
                uint64_t h = hash64(km);
                all_kmer_hashes.push_back(h);
            }
        }

        std::sort(all_kmer_hashes.begin(), all_kmer_hashes.end());
        all_kmer_hashes.erase(std::unique(all_kmer_hashes.begin(), all_kmer_hashes.end()), all_kmer_hashes.end());

        graph.bloom_filter = BlockedBloomFilter(all_kmer_hashes.size(), 0.01);
        for (uint64_t h : all_kmer_hashes) {
            graph.bloom_filter.insert(h);
        }

        if (!all_kmer_hashes.empty()) {
            graph.xor_filter = XorFilter8(all_kmer_hashes);
        }

        return graph;
    }

    KmerHit query_kmer(std::string_view kmer) const {
        KmerHit hit;
        if (kmer.size() != k) return hit;

        std::string canon = canonical_kmer(kmer);
        uint64_t h = hash64(canon);

        // Filter fast rejection
        if (!bloom_filter.contains(h) || !xor_filter.contains(h)) {
            hit.found = false;
            return hit;
        }

        // Scan unitigs for match
        std::string rc = reverse_complement(kmer);
        for (const auto& u : unitigs) {
            if (u.seq.size() < k) continue;
            for (size_t i = 0; i + k <= u.seq.size(); ++i) {
                std::string_view win(u.seq.data() + i, k);
                if (win == kmer || win == rc) {
                    hit.found = true;
                    hit.unitig_id = u.id;
                    hit.offset = static_cast<uint32_t>(i);
                    hit.color_class = u.color_class;

                    for (size_t s = 0; s < colors.num_samples(); ++s) {
                        if (colors.has_sample(u.color_class, static_cast<uint32_t>(s))) {
                            hit.present_samples.push_back(colors.sample_names()[s]);
                        }
                    }
                    return hit;
                }
            }
        }

        hit.found = false;
        return hit;
    }
};

} // namespace cdbg
